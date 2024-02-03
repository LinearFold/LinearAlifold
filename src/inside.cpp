#include <sys/time.h>
#include <algorithm>

#include "linearalifold.h"
#include "utils/ribo.h"

void BeamCKYParser::run_inside() {
    int n_seq = MSA.size();

    seq_MSA_no_gap.resize(MSA.size());
    if_tetraloops_MSA.resize(MSA.size());
    if_hexaloops_MSA.resize(MSA.size());
    if_triloops_MSA.resize(MSA.size());

#ifdef SPECIAL_HP
    for (int s = 0; s < MSA.size(); s++) {

        seq_MSA_no_gap[s] = "";

        for (auto nuc : MSA[s]) {
            if (nuc != '-') {
                seq_MSA_no_gap[s] += nuc;
            }
        }
        seq_MSA_no_gap[s] += '\0';
    }

    init_tetra_hex_tri(seq_MSA_no_gap, if_tetraloops_MSA, if_hexaloops_MSA, if_triloops_MSA);
#endif

    // start CKY decoding

    if (seq_length > 0)
        bestC[0].set(-score_external_unpaired(0, 0), MANNER_C_eq_C_plus_U);
    if (seq_length > 1)
        bestC[1].set(-score_external_unpaired(0, 1), MANNER_C_eq_C_plus_U);

    if (seq_length > 0)
        bestC[0].alpha = 0.0;
    if (seq_length > 1)
        bestC[1].alpha = 0.0;

    // for runtime statistics
    struct timeval parse_starttime, parse_endtime;
    gettimeofday(&parse_starttime, NULL);


    // from left to right
    for (int j = 0; j < seq_length; ++j) {
        // cout << "j: " << j << endl;
        // beam of H
        {
            // for every state h in H[j]
            //  1. create h(j, jnext): new hairpin loop
            //  2. generate p(i, j)
            //  3. extend h(i, j) to h(i, jnext)

            if (beam > 0 && bestH[j].size() > beam)
                beam_prune(bestH[j]);

            // 1. create h(j, jnext): new hairpin loop
            int jnext = get_cache_next_position(j, j + 3);
            while (jnext != -1 && (smart_gap[jnext] - smart_gap[j] < 4 * smart_gap_threshold)) {
                jnext = get_cache_next_position(j, jnext);
            }

            if (jnext != -1) {
                value_type newscore = 0;
                int tetra_hex_tri = -1;

                for (int s = 0; s < n_seq; s++) {
                    int nucj = SS_fast[j][s];
                    int nucj1 = s3_fast[j][s];
                    int nucjnext_1 = s5_fast[jnext][s];
                    int nucjnext = SS_fast[jnext][s];
                    int len = a2s_fast[jnext - 1][s] - a2s_fast[j][s];

#ifdef SPECIAL_HP
                    if (a2s_fast[j][s] >= 0 && a2s_fast[j][s] < seq_MSA_no_gap[s].size()) {
                        if (len == 4) { // 6:tetra
                            tetra_hex_tri = if_tetraloops_MSA[s][a2s_fast[j][s]];
                        } else if (len == 6) { // 8:hexa
                            tetra_hex_tri = if_hexaloops_MSA[s][a2s_fast[j][s]];
                        } else if (len == 3) { // 5:tri
                            tetra_hex_tri = if_triloops_MSA[s][a2s_fast[j][s]];
                        }
                    }
#endif
                    if (len < 3) {
                        newscore += -600;
                    } else
                        newscore += -score_hairpin(0, len + 1, nucj, nucj1, nucjnext_1, nucjnext, tetra_hex_tri);
                }

                if (partition_mode) {
                    Fast_LogPlusEquals(bestH[jnext][j].alpha, (newscore * inv_ktn));
                } else {
                    update_if_better(bestH[jnext][j], (newscore * inv_n), MANNER_H);
                }
            }

            sort_keys(bestH[j], keys);
            for (auto &item : keys) {
                int i = item.first;
                State &state = item.second;

                // 2. generate p(i, j)
                if (partition_mode) {
                    Fast_LogPlusEquals(bestP[j][i].alpha, state.alpha + get_pscore(i, j));
                } else {
                    update_if_better(bestP[j][i], state.alpha + get_pscore(i, j), MANNER_HAIRPIN);
                }
                // ++nos_P;

                // 3. extend h(i, j) to h(i, jnext)
                int jnext = get_cache_next_position(i, j);
                if (jnext != -1) {
                    value_type newscore = 0;
                    int tetra_hex_tri = -1;

                    for (int s = 0; s < n_seq; s++) {
                        int nuci = SS_fast[i][s];
                        int nuci1 = s3_fast[i][s];
                        int nucjnext_1 = s5_fast[jnext][s];
                        int nucjnext = SS_fast[jnext][s];
                        int len = a2s_fast[jnext - 1][s] - a2s_fast[i][s];

#ifdef SPECIAL_HP
                        if (a2s_fast[i][s] >= 0 && a2s_fast[i][s] < seq_MSA_no_gap[s].size()) {
                            if (len == 4) { // 6:tetra
                                tetra_hex_tri = if_tetraloops_MSA[s][a2s_fast[i][s]];
                            } else if (len == 6) { // 8:hexa
                                tetra_hex_tri = if_hexaloops_MSA[s][a2s_fast[i][s]];
                            } else if (len == 3) { // 5:tri
                                tetra_hex_tri = if_triloops_MSA[s][a2s_fast[i][s]];
                            }
                        }
#endif
                        if (len < 3) {
                            newscore += -600;
                        } else
                            newscore += -score_hairpin(0, len + 1, nuci, nuci1, nucjnext_1, nucjnext, tetra_hex_tri);
                    }

                    if (partition_mode) {
                        Fast_LogPlusEquals(bestH[jnext][i].alpha, (newscore * inv_ktn));
                    } else {
                        update_if_better(bestH[jnext][i], (newscore * inv_n), MANNER_H);
                    }
                }
            }
        }

        if (j == 0)
            continue;

        // beam of Multi
        {
            // for every state in Multi[j]
            //  1. generate P(i, j)
            //  2. extend (i, j) to (i, jnext)

            if (beam > 0 && bestMulti[j].size() > beam)
                beam_prune(bestMulti[j]);

            sort_keys(bestMulti[j], keys);
            for (auto &item : keys) {
                int i = item.first;
                State &state = item.second;

                // 1. generate P(i, j)
                value_type newscore = 0;
                for (int s = 0; s < n_seq; s++) {
                    int nuci = SS_fast[i][s];
                    int nuci1 = s3_fast[i][s];
                    int nucj_1 = s5_fast[j][s];
                    int nucj = SS_fast[j][s];
                    newscore += -score_multi(-1, -1, nuci, nuci1, nucj_1, nucj, -1);
                }
                // newscore += get_pscore(i, j);

                if (partition_mode) {
                    Fast_LogPlusEquals(bestP[j][i].alpha, state.alpha + (newscore * inv_ktn) + get_pscore(i, j));
                } else {
                    update_if_better(bestP[j][i], state.alpha + (newscore * inv_n) + get_pscore(i, j), MANNER_P_eq_MULTI);
                }

                // 2. extend (i, j) to (i, jnext)
                int jnext = get_cache_next_position(i, j);
                if (jnext != -1 && (smart_gap[jnext - 1] - smart_gap[j - 1] <= SINGLE_MAX_LEN)) {
                    int new_l1 = state.trace.paddings.l1;
                    int new_l2 = state.trace.paddings.l2 + jnext - j;

                    value_type newscore = 0;
                    if (!multi_approx) {
                        for (int s = 0; s < n_seq; s++)
                            newscore += -score_multi_unpaired(a2s_fast[j - 1][s], a2s_fast[jnext - 1][s]);
                    } else {
                        newscore += -(score_multi_unpaired(j - 1, jnext - 1) * n_seq); // malikap: LAF approximation
                    }

                    if (partition_mode) {
                        Fast_LogPlusEquals(bestMulti[jnext][i].alpha, state.alpha + (newscore * inv_ktn));
                    } else {
                        update_if_better(bestMulti[jnext][i], state.alpha + (newscore * inv_n), MANNER_MULTI_eq_MULTI_plus_U, new_l1, new_l2);
                    }
                }
            }
        }

        // beam of P
        {
            // for every state in P[j]
            //   1. generate new helix/bulge
            //   2. M = P
            //   3. M2 = M + P
            //   4. C = C + P
            
            if (beam > 0 && bestP[j].size() > beam)
                beam_prune(bestP[j]);

#ifdef is_cube_pruning
            bool use_cube_pruning = beam > MIN_CUBE_PRUNING_SIZE && bestP[j].size() > MIN_CUBE_PRUNING_SIZE;
#else
            bool use_cube_pruning = false;
#endif

            sort_keys(bestP[j], keys);
            for (auto &item : keys) {
                int i = item.first;
                State &state = item.second;

                // 1. P2P generate new helix / single_branch
                if (i > 0 && j < seq_length - 1) {
                    for (int p = i - 1; p >= 0 && (smart_gap[i - 1] - smart_gap[p] <= SINGLE_MAX_LEN); --p) {
                        int q = get_cache_next_position(p, j);
                        while (q != -1 && ((smart_gap[i - 1] - smart_gap[p]) + (smart_gap[q - 1] - smart_gap[j]) <= SINGLE_MAX_LEN)) {
                            // new state is of shape p..i..j..q
                            value_type newscore = 0;
                            for (int s = 0; s < n_seq; s++) {
                                int nucp = SS_fast[p][s];   // nucleotide at position p can be a gap
                                int nucq = SS_fast[q][s];   // nucleotide at position q can be a gap
                                int nucp1 = s3_fast[p][s];  // first non-gap nucleotide after position p in sequence s
                                int nucq_1 = s5_fast[q][s]; // first non-gap nucleotide before position q in sequence s
                                int nuci_1 = s5_fast[i][s]; // first non-gap nucleotide before position i in sequence s
                                int nucj1 = s3_fast[j][s];  // first non-gap nucleotide after position j in sequence s

                                // number of bases (no gaps) between position p and i in sequence s
                                int u1_local = a2s_fast[i - 1][s] - a2s_fast[p][s];

                                // number of bases (no gaps) between position j and q in sequence s
                                int u2_local = a2s_fast[q - 1][s] - a2s_fast[j][s];

                                int type = NUM_TO_PAIR(nucp, nucq); // could be -- or X- or -X
                                int type2 = NUM_TO_PAIR(SS_fast[j][s], SS_fast[i][s]);
                                newscore += -score_single_alifold(u1_local, u2_local, type, type2, nucp1, nucq_1, nuci_1, nucj1);
                            }
                            // newscore += get_pscore(p, q);

                            if (partition_mode) {
                                Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + (newscore * inv_ktn) + get_pscore(p, q));
                            } else if (p == i - 1 && q == j + 1) {
                                update_if_better(bestP[q][p], state.alpha + (newscore * inv_n) + get_pscore(p, q), MANNER_HELIX);
                            } else {
                                update_if_better(bestP[q][p], state.alpha + (newscore * inv_n) + get_pscore(p, q), MANNER_SINGLE, (i - p), (q - j));
                            }

                            // malikap: generate a next q, and generate a new helix / single_branch
                            q = get_cache_next_position(p, q);
                        }
                    }
                }
            

                // 2. M = P
                if (i > 0 && j < seq_length - 1) {
                    value_type newscore = 0;
                    for (int s = 0; s < n_seq; s++) {
                        int nuci_1 = (i - 1) > -1 ? s5_fast[i][s] : -1;
                        int nuci = SS_fast[i][s];
                        int nucj = SS_fast[j][s];
                        int nucj1 = (j + 1) < seq_length ? s3_fast[j][s] : -1;
                        newscore += -score_M1(-1, -1, -1, nuci_1, nuci, nucj, nucj1, -1); // no position information needed
                    }

                    if (partition_mode) {
                        Fast_LogPlusEquals(bestM[j][i].alpha, state.alpha + (newscore * inv_ktn));
                    } else {
                        update_if_better(bestM[j][i], state.alpha + (newscore * inv_n), MANNER_M_eq_P);
                    }
                }
                // 3. M2 = M + P
                if (!use_cube_pruning || partition_mode) {
                    int k = i - 1;
                    if (k > 0 && !bestM[k].empty()) {
                        value_type M1_score = 0;
                        for (int s = 0; s < n_seq; s++) {
                            int nuci_1 = (i - 1) > -1 ? s5_fast[i][s] : -1; // TODO, need to check boundary?
                            int nuci = SS_fast[i][s];
                            int nucj = SS_fast[j][s];
                            int nucj1 = (j + 1) < seq_length ? s3_fast[j][s] : -1; // TODO, need to check boundary? it may diff. from RNAalifold
                            M1_score += -score_M1(-1, -1, -1, nuci_1, nuci, nucj, nucj1, -1);
                        }
                        M1_score *= inv_n;

#ifndef is_candidate_list
                        for (auto &m : bestM[k]) {
                            int newi = m.first;
                            // eq. to first convert P to M1, then M2/M = M + M1
                            // value_type newscore = M1_score + m.second.score;

                            if (partition_mode) {
                                Fast_LogPlusEquals(bestM2[j][newi].alpha, state.alpha + m.second.alpha + (M1_score * inv_kt));
                            } else {
                                update_if_better(bestM2[j][newi], state.alpha + m.second.alpha + M1_score, MANNER_M2_eq_M_plus_P, k);
                            }
                        }
#else 
                        // candidate list
                        auto bestM2_iter = bestM2[j].find(i);
                        if (bestM2_iter == bestM2[j].end() || M1_score > bestM2_iter->second.alpha || partition_mode) {
                            for (auto &m : bestM[k]) {
                                int newi = m.first;
                                // eq. to first convert P to M1, then M2/M = M + M1
                                // value_type newscore = M1_score + m.second.score;
                                if (partition_mode) {
                                    Fast_LogPlusEquals(bestM2[j][newi].alpha, state.alpha + m.second.alpha + (M1_score * inv_kt));
                                } else {
                                    update_if_better(bestM2[j][newi], state.alpha + m.second.alpha + M1_score, MANNER_M2_eq_M_plus_P, k);
                                }
                            }
                        }
#endif
                    }
                }
                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                        State &prefix_C = bestC[k];
                        if (prefix_C.manner != MANNER_NONE || partition_mode) {
                            value_type newscore = 0;
                            for (int s = 0; s < n_seq; s++) {
                                int nuci_1 = (i - 1) > -1 ? s5_fast[i][s] : -1; // external.c line 1165, weird
                                int nuci = SS_fast[i][s];
                                int nucj = SS_fast[j][s];
                                int nucj1 = (j + 1) < seq_length ? s3_fast[j][s] : -1; // external.c line 1165, weird
                                newscore += -score_external_paired(-1, -1, nuci_1, nuci, nucj, nucj1, -1);
                            }

                            if (partition_mode) {
                                Fast_LogPlusEquals(bestC[j].alpha, state.alpha + prefix_C.alpha + (newscore * inv_ktn));
                            } else {
                                update_if_better(bestC[j], state.alpha + prefix_C.alpha + (newscore * inv_n), MANNER_C_eq_C_plus_P, k);
                            }
                        }
                    } 
                    else 
                    {
                        value_type newscore = 0;
                        for (int s = 0; s < MSA.size(); s++) {
                            int nuci = SS_fast[i][s];
                            int nucj = SS_fast[j][s];
                            int nucj1 = (j + 1) < seq_length ? s3_fast[j][s] : -1; // external.c line 1165, weird
                            newscore += -score_external_paired(0, j, -1, nuci, nucj, nucj1, -1);
                        }

                        if (partition_mode) {
                            Fast_LogPlusEquals(bestC[j].alpha, state.alpha + (newscore * inv_ktn));
                        } else {
                            update_if_better(bestC[j], state.alpha + (newscore * inv_n), MANNER_C_eq_C_plus_P, -1);
                        }
                    }
                }
            }

            // lhuang: check j < n-1
            if (use_cube_pruning && j < seq_length - 1 && !partition_mode) {
                // 3. M2 = M + P with cube pruning
                vector<int> valid_Ps;
                vector<value_type> M1_scores;

                sort_keys(bestP[j], keys);
                for (auto &item : keys) {

                    int i = item.first;
                    State &state = item.second;

                    int k = i - 1;

                    // group candidate Ps
                    if (k > 0 && !bestM[k].empty()) {
                        assert(bestM[k].size() == sorted_bestM[k].size());
                        value_type M1_score = 0;

                        for (int s = 0; s < n_seq; s++) {
                            int nucj = SS_fast[j][s];
                            int nucj1 = (j + 1) < seq_length ? s3_fast[j][s] : -1;
                            int nuci = SS_fast[i][s];
                            int nuci_1 = (i - 1 > -1) ? s5_fast[i][s] : -1;
                            M1_score += -score_M1(-1, -1, -1, nuci_1, nuci, nucj, nucj1, -1);
                        }

                        M1_score *= inv_n;
                        M1_score += state.alpha;

                        auto bestM2_iter = bestM2[j].find(i);
#ifndef is_candidate_list
                        valid_Ps.push_back(i);
                        M1_scores.push_back(M1_score);
#else
                        if (bestM2_iter == bestM2[j].end() || M1_score > bestM2_iter->second.alpha) {
                            valid_Ps.push_back(i);
                            M1_scores.push_back(M1_score);
                        }
#endif
                    }
                }

                // build max heap
                // heap is of form (heuristic score, (index of i in valid_Ps, index of M in bestM[i-1]))
                vector<pair<value_type, pair<int, int>>> heap;
                for (int p = 0; p < valid_Ps.size(); ++p) {
                    int i = valid_Ps[p];
                    int k = i - 1;
                    heap.push_back(make_pair(M1_scores[p] + sorted_bestM[k][0].first, make_pair(p, 0)));
                    push_heap(heap.begin(), heap.end());
                }

                // start cube pruning
                // stop after beam size M2 states being filled
                int filled = 0;
                // exit when filled >= beam and current score < prev score
                value_type prev_score = VALUE_MIN;
                value_type current_score = VALUE_MIN;
                while ((filled < beam || current_score == prev_score) && !heap.empty()) {
                    auto &top = heap.front();
                    prev_score = current_score;
                    current_score = top.first;
                    int index_P = top.second.first;
                    int index_M = top.second.second;
                    int i = valid_Ps[top.second.first];
                    int k = i - 1;
                    int newi = sorted_bestM[k][index_M].second;
                    value_type newscore = M1_scores[index_P] + bestM[k][newi].alpha;
                    pop_heap(heap.begin(), heap.end());
                    heap.pop_back();

                    if (bestM2[j][newi].manner == MANNER_NONE) {
                        ++filled;
                        update_if_better(bestM2[j][newi], newscore, MANNER_M2_eq_M_plus_P, k);
                        // ++nos_M2;
                    } else {
                        assert(bestM2[j][newi].alpha > newscore - 1e-8);
                    }

                    ++index_M;
                    while (index_M < sorted_bestM[k].size()) {
                        // candidate_score is a heuristic score
                        value_type candidate_score = M1_scores[index_P] + sorted_bestM[k][index_M].first;
                        int candidate_newi = sorted_bestM[k][index_M].second;
                        if (bestM2[j].find(candidate_newi) == bestM2[j].end()) {
                            heap.push_back(make_pair(candidate_score,
                                                     make_pair(index_P, index_M)));
                            push_heap(heap.begin(), heap.end());
                            break;
                        } else {
                            // based on the property of cube pruning, the new score must be worse
                            // than the state already inserted so we keep iterate through the candidate
                            // list to find the next candidate
                            ++index_M;
                            assert(bestM2[j][candidate_newi].alpha >
                                   M1_scores[index_P] + bestM[k][candidate_newi].alpha - 1e-8);
                        }
                    }
                }
            }
        }

        // beam of M2
        {
            // for every state in M2[j]
            //  1. M = M2
            //  2. multi-loop (by extending M2 on the left)

            if (beam > 0 && bestM2[j].size() > beam)
                beam_prune(bestM2[j]);

            sort_keys(bestM2[j], keys);
            for (auto &item : keys) {
                int i = item.first;
                State &state = item.second;

                // 1. M = M2
                if (partition_mode) {
                    Fast_LogPlusEquals(bestM[j][i].alpha, state.alpha);
                } else {
                    update_if_better(bestM[j][i], state.alpha, MANNER_M_eq_M2);
                }
                // ++nos_M;

                // 2. multi-loop (by extending M2 on the left)
                for (int p = i - 1; p >= 0 && (smart_gap[i - 1] - smart_gap[p] <= SINGLE_MAX_LEN); --p) {
                    int q = get_cache_next_position(p, j);
                    if (q != -1) {
                        // the current shape is p..i M2 j..q
                        value_type newscore = 0;
                        if (!multi_approx) {
                            for (int s = 0; s < n_seq; s++)
                                newscore += -(score_multi_unpaired(a2s_fast[p][s], a2s_fast[i - 1][s]) + score_multi_unpaired(a2s_fast[j][s], a2s_fast[q - 1][s]));
                        } else {
                            newscore += -((score_multi_unpaired(p, i - 1) + score_multi_unpaired(j, q - 1)) * n_seq); // malikap: LinAliFold's approximation
                        }

                        if (partition_mode) {
                            Fast_LogPlusEquals(bestMulti[q][p].alpha, state.alpha + (newscore * inv_ktn));
                        } else {
                            update_if_better(bestMulti[q][p], state.alpha + (newscore * inv_n), MANNER_MULTI, (i - p), (q - j));
                        }
                    }
                }
            }
        }

        // beam of M
        {
            // for every state in M[j]
            //  1. M = M + unpaired

            value_type threshold = VALUE_MIN;
            if (beam > 0 && bestM[j].size() > beam)
                threshold = beam_prune(bestM[j]);
#ifdef is_cube_pruning
            if (!partition_mode)
                sortM(threshold, bestM[j], sorted_bestM[j]);
#endif
            sort_keys(bestM[j], keys);
            for (auto &item : keys) {
                int i = item.first;
                State &state = item.second;

                if (j < seq_length - 1) {
                    value_type newscore = 0;
                    if (!multi_approx) {
                        for (int s = 0; s < n_seq; s++)
                            newscore += -score_multi_unpaired(a2s_fast[j][s], a2s_fast[j + 1][s]);
                    } else {
                        newscore += -(score_multi_unpaired(j, j + 1) * n_seq); // malikap: LinAliFold's approximation
                    }

                    if (partition_mode) {
                        Fast_LogPlusEquals(bestM[j + 1][i].alpha, state.alpha + (newscore * inv_ktn));
                    } else {
                        update_if_better(bestM[j + 1][i], state.alpha + (newscore * inv_n), MANNER_M_eq_M_plus_U);
                    }
                }
            }
        }

        // beam of C
        {
            // C = C + U
            if (j < seq_length - 1) {
                if (partition_mode) {
                    Fast_LogPlusEquals(bestC[j + 1].alpha, bestC[j].alpha);
                } else {
                    update_if_better(bestC[j + 1], bestC[j].alpha, MANNER_C_eq_C_plus_U);
                }
            }
        }
    } // end of for-loo j

    State &viterbi = bestC[seq_length - 1];
    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec - parse_starttime.tv_usec) / 1000000.0;
    inside_time = parse_elapsed_time;
    
    if (partition_mode) {
        fprintf(stdout, "Free Energy of Ensemble: %.2f kcal/mol\n", -viterbi.alpha / 100.0 / MSA.size());
        if (is_verbose)
            fprintf(stdout, "Partition Function Calculation Time: %.2f seconds.\n", parse_elapsed_time);
    } else {
        fprintf(stdout, "Minimum Free Energy: %.2f kcal/mol\n", -viterbi.alpha / 100.0);
        if (is_verbose)
            fprintf(stdout, "MFE Structure Calculation Time: %.2f seconds.\n", parse_elapsed_time);

        char res[seq_length + 1];
        get_parentheses(res, MSA[0]);
        string result = string(res);

        double printscore = (viterbi.alpha / -100.0);
        auto result_pairs = get_pairs(result);
        float pscore_f = 0.;
        for (auto &pair : result_pairs) {
            pscore_f += pscore[pair.first][pair.second - (pair.first + 1)];     // malikap: adjusted 2nd dimesion index
        }
        pscore_f = -pscore_f / 100.;

        printf("\nMFE Structure: \n");
        printf("%s (%.2f = %.2f + %.2f)\n", result.c_str(), printscore, printscore - pscore_f, pscore_f);
        // unsigned long nos_tot = nos_H + nos_P + nos_M2 + nos_Multi + nos_M + nos_C;
        // cout << nos_tot << endl;
    }
}