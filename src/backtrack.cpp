#include "linearalifold.h"
#include <algorithm>


using namespace std;

void BeamCKYParser::backtrack(State *state, string &structure, bool best_only, default_random_engine *generator) {
    if (state == NULL || state->type == TYPE_MAX) {
        return;
    }

    if (state->type == TYPE_P) {
        structure[state->idx.first] = '(';
        structure[state->idx.second] = ')';
    }

    // for MFE backtracking
    if (best_only) {
        vector<HEdge> edge = get_incoming_hedges(state, inv_n, true);
        backtrack(edge[0].left, structure);
        backtrack(edge[0].right, structure);
    }
    // for sampling mode backtracking
    else {
        StateKey state_key = StateKey(state->idx.first, state->idx.second, state->type);
        if (state_hedges_cache.find(state_key) == state_hedges_cache.end()) {
            state_hedges_cache[state_key] = get_incoming_hedges(state, inv_ktn, false);
        }
        vector<HEdge> edges = state_hedges_cache[state_key];

        if (edges.size() == 0) return;

        vector<value_type> weights;
        for (auto &edge : edges) {
            value_type val = Fast_Exp(edge.left->alpha + edge.weight + (edge.right ? edge.right->alpha : 0) - state->alpha);
            weights.push_back(val);
        }

        // sample an edge
        discrete_distribution<> distribution(weights.begin(), weights.end());
        int idx = distribution(*generator);

        backtrack(edges[idx].left, structure, false, generator);
        backtrack(edges[idx].right, structure, false, generator);
    }
}


void BeamCKYParser::backward_update(int i, int j, State &state, Type type, float global_threshold, int *pruned, int *saved) {
    state.type = type;
    state.idx = make_pair(i, j);

    vector<HEdge> incoming_hedges = get_incoming_hedges(&state, inv_ktn, false);
    if (incoming_hedges.empty()) return;
    vector<HEdge*> saved_hedges;
    HEdge *best_hedge = NULL;

    float edge_threshold = global_threshold - state.beta;
    double best_inside = VALUE_MIN, saved_inside = VALUE_MIN;
    int local_pruned = 0;

    for (auto &hedge: incoming_hedges) {
        double edge_inside = hedge.left->alpha + hedge.weight + (hedge.right ? hedge.right->alpha : 0);

        if (edge_inside > edge_threshold) {
            Fast_LogPlusEquals(saved_inside, edge_inside);
            saved_hedges.push_back(&hedge);
        } else {
            local_pruned++;
            if (saved_hedges.empty() && edge_inside > best_inside) {
                best_inside = edge_inside;
                best_hedge = &hedge;
            }
        }
    }

    *pruned += local_pruned; // global
    float delta;             // scaling factor to compensate for edge pruning
    if (!saved_hedges.empty()) {
        delta = state.alpha - saved_inside;
    } else {        // all pruned: use best hyperedge
        delta = state.alpha - best_inside;
        saved_hedges.push_back(best_hedge); 
        *pruned -= 1;  // one more edge recovered
    }

    for (auto &hedge : saved_hedges) {
        State *left = hedge->left, *right = hedge->right; 
        if (!right) // edge
            left->logplus(state.beta + hedge->weight + delta);
        else {      // hyperedge
            left->logplus(state.beta + hedge->weight + delta + right->alpha);
            right->logplus(state.beta + hedge->weight + delta + left->alpha);
        }
    }
    *saved += saved_hedges.size();
}


vector<HEdge> BeamCKYParser::get_incoming_hedges(State *state, double multiplier, bool best_only = true) {
    vector<HEdge> incoming_edges;

    if (best_only) {
        incoming_edges.push_back(HEdge());
    }

    auto update = [&incoming_edges, best_only] (float edge_weight, State *left, State *right = NULL) {
        if (best_only) {  // for backtrack
            value_type value = edge_weight + left->alpha;
            if (right != NULL) value += right->alpha;

            if (value > incoming_edges[0].weight) {
                incoming_edges[0].set(left, right, value);
            }
        }  
        else {  // for outside
            incoming_edges.push_back(HEdge(left, right, edge_weight));
        }
    };


    int i = state->idx.first;
    int j = state->idx.second;
    Type type = state->type;    

    switch (type) {
    case TYPE_C: {
        if (j < 1) break;
        
        // C = C + U
        float edge_weight = 0;
        bestC[j - 1].set_attributes(0, j - 1, TYPE_C);
        update(edge_weight, &bestC[j - 1]);

        // C = C + P
        for (auto &item : bestP[j]) {
            int i = item.first;
            // the shape is C(0,j) = C(0,i-1) + P(i,j)
            float edge_weight = 0;
            for (int s = 0; s < MSA.size(); s++) {
                int nuci_1 = (i - 1) > -1 ? s5_fast[i][s] : -1;
                int nuci = SS_fast[i][s];
                int nucj = SS_fast[j][s];
                int nucj1 = (j + 1) < seq_length ? s3_fast[j][s] : -1;
                int len = a2s_fast[j - 1][s] - a2s_fast[i][s];

                edge_weight += -score_external_paired(a2s_fast[i][s], a2s_fast[j][s], nuci_1, nuci, nucj, nucj1, len);
            }
                
            edge_weight *= multiplier;
            if (i - 1 > -1) bestC[i - 1].set_attributes(0, i - 1, TYPE_C);
            item.second.set_attributes(i, j, TYPE_P);   
            update(edge_weight, &item.second, i - 1 > -1 ? &bestC[i - 1] : NULL);
        }

    } break;

    case TYPE_P: {

        // P = H
        auto itr = bestH[j].find(i);
        if (itr != bestH[j].end())  {
            float edge_weight = get_pscore(i, j);
            // edge_weight *= multiplier;
            itr->second.set_attributes(i, j, TYPE_MAX);
            update(edge_weight, &itr->second);
        }

        // P = P + P
        for(int p = i + 1; p < j && smart_gap[p - 1] - smart_gap[i] <= SINGLE_MAX_LEN; p++) {
            int q = get_cache_prev_position(p, j);
            while (q != -1 && p < q && (smart_gap[p - 1] - smart_gap[i]) + (smart_gap[j - 1] - smart_gap[q]) <= SINGLE_MAX_LEN) {
                auto itr = bestP[q].find(p);
                if (itr != bestP[q].end()) {
                    // current shape is i...p (pair) q...j
                    float edge_weight = 0; 
                    for (int s = 0; s < MSA.size(); s++) {
                        int nuci = SS_fast[i][s];
                        int nuci1 = s3_fast[i][s];  // malikap: no need to check boundary
                        int nucj_1 = s5_fast[j][s]; // malikap: no need to check boundary
                        int nucj = SS_fast[j][s];

                        int nucp_1 = s5_fast[p][s];  // malikap: no need to check boundary
                        int nucp = SS_fast[p][s];
                        int nucq = SS_fast[q][s];
                        int nucq1 = s3_fast[q][s];  // malikap: no need to check boundary

                        int pair_type_closing = NUM_TO_PAIR(nuci, nucj);
                        int pair_type_enclosed = NUM_TO_PAIR(nucq, nucp);

                        int five_end_loop_size = a2s_fast[p - 1][s] - a2s_fast[i][s];
                        int three_end_loop_size = a2s_fast[j - 1][s] - a2s_fast[q][s];

                        edge_weight += -score_single_alifold(five_end_loop_size, three_end_loop_size, pair_type_closing, pair_type_enclosed, nuci1, nucj_1, nucp_1, nucq1);
                    }
                    edge_weight *= multiplier;
                    edge_weight += get_pscore(i, j);
                    itr->second.set_attributes(p, q, TYPE_P);
                    update(edge_weight, &itr->second);
                }
                q = get_cache_prev_position(p, q);
            }
        }
        // for(int q = j - 1; q > i && (smart_gap[j - 1] - smart_gap[q]) <= SINGLE_MAX_LEN; q--) {
        //     int p = get_cache_next_position(q, i);
        //     while(p != -1 && p < q && (smart_gap[p - 1] - smart_gap[i]) + (smart_gap[j - 1] - smart_gap[q]) <= SINGLE_MAX_LEN) {
                    // ....
        //         p = get_cache_next_position(q, p);
        //     }
        // }

        // P = Multi
        itr = bestMulti[j].find(i);
        if (itr != bestMulti[j].end()) {
            float edge_weight = 0;
            for (int s = 0; s < MSA.size(); s++) {
                int nuci = SS_fast[i][s];
                int nuci1 = s3_fast[i][s];   // malikap: no need to check boundary
                int nucj_1 = s5_fast[j][s];  // malikap: no need to check boundary
                int nucj = SS_fast[j][s];
                int len = a2s_fast[j - 1][s] - a2s_fast[i][s];

                edge_weight += -score_multi(a2s_fast[i][s], a2s_fast[j][s], nuci, nuci1, nucj_1, nucj, len);
            }
            edge_weight *= multiplier;
            edge_weight += get_pscore(i, j);
            itr->second.set_attributes(i, j, TYPE_MULTI);
            update(edge_weight, &itr->second);
        }
    } break;

    case TYPE_M: {
        // M = M + U
        if (j > 0) {
            auto itr = bestM[j - 1].find(i);
            if (itr != bestM[j - 1].end()) {
                float edge_weight = 0;
                if (!multi_approx) {
                    for (int s = 0; s < MSA.size(); s++)
                        edge_weight += -score_multi_unpaired(a2s_fast[j - 1][s], a2s_fast[j][s]); // malikap: M1 unparied base score computation
                } else {
                    edge_weight = -(score_multi_unpaired(j - 1, j) * MSA.size()); // malikap: approximate number of multi bases
                }

                edge_weight *= multiplier;
                itr->second.set_attributes(i, j - 1, TYPE_M);
                update(edge_weight, &itr->second);
            }
        }

        // M = P
        auto itr = bestP[j].find(i);
        if (itr != bestP[j].end()) {
            float edge_weight = 0;
            for (int s = 0; s < MSA.size(); s++) {
                int nuci_i = (i - 1) > -1 ? s5_fast[i][s] : -1;
                int nuci = SS_fast[i][s];
                int nucj = SS_fast[j][s];
                int nucj1 = (j + 1) < seq_length ? s3_fast[j][s] : -1;

                edge_weight += -score_M1(i, j, j, nuci_i, nuci, nucj, nucj1, seq_length);
            }

            edge_weight *= multiplier;
            itr->second.set_attributes(i, j, TYPE_P);
            update(edge_weight, &itr->second);
        }

        // M = M2
        itr = bestM2[j].find(i);
        if (itr != bestM2[j].end()) {
            float edge_weight = 0;
            itr->second.set_attributes(i, j, TYPE_M2);
            update(edge_weight, &itr->second);
        }
    } break;


    case TYPE_M2: {
        if (sortedP[j].size() == 0) {
            // lhuang: M2 might be short, and b might be big
            // so only explore P that fits M2 = M + P
            for (auto const &item : bestP[j])
                sortedP[j].push_back(-item.first);
            sort(sortedP[j].begin(), sortedP[j].end());
        }

        // M2 = M + P
        for (auto &item : sortedP[j]) {
            int k = -item;
            if (k > i) {
                auto itr = bestM[k - 1].find(i);
                if (itr != bestM[k - 1].end()) {
                    // M2(i,j) = M(i,k-1) + P(k,j)
                    float edge_weight = 0;
                    for (int s = 0; s < MSA.size(); s++) {
                        int nuck_1 = s5_fast[k][s];
                        int nuck = SS_fast[k][s];
                        int nucj = SS_fast[j][s];
                        int nucj1 = (j + 1 < seq_length) ? s3_fast[j][s] : -1;

                        edge_weight += -score_M1(k, j, j, nuck_1, nuck, nucj, nucj1, seq_length);
                    }

                    edge_weight *= multiplier;
                    itr->second.set_attributes(i, k - 1, TYPE_M);
                    bestP[j][k].set_attributes(k, j, TYPE_P);
                    update(edge_weight, &itr->second, &bestP[j][k]);
                }
            }
            else {
                break;
            }
        }
    } break;

    case TYPE_MULTI: {
        // Multi = Multi + extend-right (unpaired) (multi jump right)
        int jprev = get_cache_prev_position(i, j);
        if (jprev <= i) break;

        if ((smart_gap[j - 1] - smart_gap[jprev - 1]) <= SINGLE_MAX_LEN) {
            auto itr = bestMulti[jprev].find(i);
            if (itr != bestMulti[jprev].end()) {
                float edge_weight = 0;
                if (!multi_approx) {
                    for (int s = 0; s < MSA.size(); s++)
                        edge_weight += -score_multi_unpaired(a2s_fast[jprev - 1][s], a2s_fast[j - 1][s]);
                } else {
                    edge_weight += -score_multi_unpaired(jprev - 1, j - 1);
                }
                edge_weight *= multiplier;
                itr->second.set_attributes(i, jprev, TYPE_MULTI);
                update(edge_weight, &itr->second);
            }
        }
    

        // Multi = extend-left + M2 + extend-right
        for (int q = j - 1; q >= jprev; q--) {
            for(int p = i + 1; p < q && (smart_gap[p - 1] - smart_gap[i]) <= SINGLE_MAX_LEN; p++) {
                auto itr = bestM2[q].find(p);
                if (itr != bestM2[q].end()) {
                    // the current shape is i..p M2 q..j
                    float edge_weight = 0;
                    if (!multi_approx) {
                        for (int s = 0; s < MSA.size(); s++)
                            edge_weight += -score_multi_unpaired(a2s_fast[i][s], a2s_fast[p - 1][s]) + -score_multi_unpaired(a2s_fast[q][s], a2s_fast[j - 1][s]);
                    } else {
                        edge_weight += -score_multi_unpaired(i, p - 1) + -score_multi_unpaired(q, j - 1); // malikap: unpaired bases approximation
                    }
                    edge_weight *= multiplier;
                    itr->second.set_attributes(p, q, TYPE_M2);
                    update(edge_weight, &itr->second);
                }

            }
        }
    } break;

    default:
        break;
    }

    return incoming_edges;
}