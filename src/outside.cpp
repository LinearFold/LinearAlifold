#include <sys/time.h>
#include <algorithm>

#include "linearalifold.h"
#include "utils/ribo.h"

using namespace std;


string type2str[5] = {"C", "M", "M2", "P", "MULTI"};


void BeamCKYParser::outside_alifold() {

    bestC[seq_length - 1].beta = 0.0;
    int n_seq = MSA.size();

    struct timeval start_time, end_time;
    gettimeofday(&start_time, NULL);

    // from right to left
    for (int j = seq_length - 1; j > 0; --j) {
        // beam of C
        {
            // C = C + U
            if (j < seq_length - 1)
            {
                Fast_LogPlusEquals(bestC[j].beta, bestC[j + 1].beta);
            }
        }

        // beam of M
        {
            for (auto &item : bestM[j])
            {
                int i = item.first;
                State &state = item.second;
                if (j < seq_length - 1)
                {
                    value_type newscore = 0;
                    if (!multi_approx) {
                        for (int s = 0; s < n_seq; s++)
                            newscore += -score_multi_unpaired(a2s_fast[j][s], a2s_fast[j + 1][s]);
                    }
                    else {
                        newscore += -(score_multi_unpaired(j, j + 1) * n_seq); // malikap: M1 unparied base score computation
                    }

                    Fast_LogPlusEquals(state.beta, bestM[j + 1][i].beta + (newscore * inv_ktn));
                }
            }
        }

        // beam of M2
        {
            // for every state in M2[j]
            //  1. multi-loop (by extending M2 on the left)
            //  2. M = M2
            
            for (auto &item : bestM2[j]) {
                int i = item.first;
                State &state = item.second;

                // 1. multi-loop (by extending M2 on the left)
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
                        Fast_LogPlusEquals(state.beta, bestMulti[q][p].beta + (newscore * inv_ktn));
                    }
                }

                // 2. M = M2
                Fast_LogPlusEquals(state.beta, bestM[j][i].beta);
            }
        }

        // beam of P
        {
            for(auto &item: bestP[j]) {
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
                
                            Fast_LogPlusEquals(state.beta, bestP[q][p].beta + (newscore * inv_ktn) + get_pscore(p, q));

                            // malikap: generate a next q, and generate a new helix / single_branch
                            q = get_cache_next_position(p, q);
                        }
                    }
                }

                // 2. M = P
                if (i > 0 && j < seq_length - 1) {
                    value_type newscore = 0;
                    for (int s = 0; s < n_seq; s++) {
                        int nuci_1 = (i - 1 > -1) ? s5_fast[i][s] : -1;
                        int nuci = SS_fast[i][s];
                        int nucj = SS_fast[j][s];
                        int nucj1 = (j + 1) < seq_length ? s3_fast[j][s] : -1;
                        newscore += -score_M1(-1, -1, -1, nuci_1, nuci, nucj, nucj1, -1); // no position information needed
                    }

                    Fast_LogPlusEquals(state.beta, bestM[j][i].beta + (newscore * inv_ktn));
                }

                // 3. M2 = M + P
                int k = i - 1;
                if (k > 0 && !bestM[k].empty()) {
                    value_type newscore = 0;
                    for (int s = 0; s < n_seq; s++) {
                        int nuci_1 = (i - 1) > -1 ? s5_fast[i][s] : -1; // TODO, need to check boundary?
                        int nuci = SS_fast[i][s];
                        int nucj = SS_fast[j][s];
                        int nucj1 = (j + 1) < seq_length ? s3_fast[j][s] : -1; // TODO, need to check boundary? it may diff. from RNAalifold
                        newscore += -score_M1(-1, -1, -1, nuci_1, nuci, nucj, nucj1, -1);
                    }

                    value_type m1_alpha = (newscore * inv_ktn);
                    value_type m1_plus_P_alpha = state.alpha + m1_alpha;
                    for (auto &m : bestM[k])
                    {
                        int newi = m.first;
                        State &m_state = m.second;
                        Fast_LogPlusEquals(state.beta, (bestM2[j][newi].beta + m_state.alpha + m1_alpha));
                        Fast_LogPlusEquals(m_state.beta, (bestM2[j][newi].beta + m1_plus_P_alpha));
                    }
                }

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                        value_type newscore = 0;
                        for (int s = 0; s < MSA.size(); s++)
                        {
                            int nuci_1 = (i - 1) > -1 ? s5_fast[i][s] : -1; // external.c line 1165, weird
                            int nuci = SS_fast[i][s];
                            int nucj = SS_fast[j][s];
                            int nucj1 = (j + 1) < seq_length ? s3_fast[j][s] : -1; // external.c line 1165, weird

                            newscore += (-score_external_paired(-1, -1, nuci_1, nuci, nucj, nucj1, -1));
                        }
                        
                        Fast_LogPlusEquals(bestC[k].beta, state.alpha + bestC[j].beta + (newscore * inv_ktn));
                        Fast_LogPlusEquals(state.beta, bestC[k].alpha + bestC[j].beta + (newscore * inv_ktn));
                    } 
                    else {
                        value_type newscore = 0;
                        for (int s = 0; s < MSA.size(); s++)
                        {
                            int nuci = SS_fast[i][s];
                            int nucj = SS_fast[j][s];
                            int nucj1 = (j + 1) < seq_length ? s3_fast[j][s] : -1; // external.c line 1165, weird

                            newscore += -score_external_paired(0, j, -1, nuci, nucj, nucj1, -1);
                        }

                        Fast_LogPlusEquals(state.beta, bestC[j].beta + (newscore * inv_ktn));
                    }
                }
            }
        }


        // beam of Multi
        {
            // for every state in Multi[j]
            //  1. extend (i, j) to (i, jnext)
            //  2. generate P(i, j)
            
            for (auto &item : bestMulti[j]) {
                int i = item.first;
                State &state = item.second;

                // 1. extend (i, j) to (i, jnext)
                int jnext = get_cache_next_position(i, j);
                if (jnext != -1 && (smart_gap[jnext - 1] - smart_gap[j - 1] <= SINGLE_MAX_LEN)) {
                    value_type newscore = 0;
                    if (!multi_approx) {
                        for (int s = 0; s < n_seq; s++)
                            newscore += -score_multi_unpaired(a2s_fast[j - 1][s], a2s_fast[jnext - 1][s]);
                    } else {
                        newscore += -(score_multi_unpaired(j - 1, jnext - 1) * n_seq); // malikap: LAF approximation
                    }

                    Fast_LogPlusEquals(state.beta, bestMulti[jnext][i].beta + (newscore * inv_ktn));
                }

                // 2. generate P(i, j)
                value_type newscore = 0;
                for (int s = 0; s < n_seq; s++) {
                    int nuci = SS_fast[i][s];
                    int nuci1 = s3_fast[i][s];
                    int nucj_1 = s5_fast[j][s];
                    int nucj = SS_fast[j][s];
                    newscore += -score_multi(-1, -1, nuci, nuci1, nucj_1, nucj, -1);
                }

                // newscore += get_pscore(i, j);
                Fast_LogPlusEquals(state.beta, bestP[j][i].beta + (newscore * inv_ktn) + get_pscore(i, j));
            }
        }

        // // beam of H
        // {
        //     for(auto &item: bestH[j]) {
        //         // generate P(i, j)
        //         int i = item.first;
        //         State &state = item.second;
        //         value_type newscore = get_pscore(i, j);
        //         Fast_LogPlusEquals(state.beta, bestP[j][i].beta + (newscore * inv_ktn));
        //     }
        // }
    } // end of for-loop j

    gettimeofday(&end_time, NULL);
    double bpp_elapsed_time = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - start_time.tv_usec) / 1000000.0;
    if (is_verbose)
        fprintf(stdout, "Base Pairing Probabilities Calculation Time: %.2f seconds.\n", bpp_elapsed_time);
    fflush(stdout);
}

void BeamCKYParser::run_lazy_outside() {
    struct timeval parse_starttime, parse_endtime;
    gettimeofday(&parse_starttime, NULL);
    bestC[seq_length - 1].beta = 0.0;
    // bestC[0].alpha = 0.0;
    float global_threshold = bestC[seq_length - 1].alpha - deviation_threshold; // Default: -9.91152;
    int tot = 0, vis = 0, weird = 0;
    int pruned = 0, saved = 0;
    int edge_counts[6] = {0, 0, 0, 0, 0, 0};
    memset(edge_counts, 0, sizeof(edge_counts));
    
    // from right to left
    for (int j = seq_length - 1; j > 0; --j) {
        if (bestC[j].beta > -deviation_threshold) {        
            edge_counts[TYPE_C]++;
            backward_update(0, j, bestC[j], TYPE_C, global_threshold, &pruned, &saved);
            vis++;
        }
        tot++;

        for (int type = TYPE_M; type < TYPE_MAX; type++) { // reverse topol order: C->M->M2->P->Multi
            for (auto &item : best_states[type]->at(j)) {
                int i = item.first;
                State &state = item.second;
                if (state.beta > -deviation_threshold) {
                    edge_counts[type]++;
                    backward_update(i, j, state, (Type)type, global_threshold, &pruned, &saved);
                    vis++;
                }
                tot++;
            }
        }   
    }

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec - parse_starttime.tv_usec) / 1000000.0;

    if (is_verbose) {
        float tot_edges = pruned + saved;
        printf("\nLazy Outside Time: %.2f seconds (%.1f%%) (visited edges: pruned %d, saved %d) (weird %d).\n",
               parse_elapsed_time, 100. * parse_elapsed_time / inside_time,
               pruned, saved, weird);
        printf("nodes: visited %d; total %d (%.1f%%)\n", vis, tot, vis * 100. / tot);
        printf("edge counts: ");
        for (int t = 0; t < 5; t++)
            printf("%s: %d (%.1f%%)  ", type2str[t].c_str(), edge_counts[t], (edge_counts[t] * 100. / tot_edges));
        printf("\n");
    }

    fflush(stdout);
}


void BeamCKYParser::cal_PairProb(State &viterbi)
{
    double kTn = double(kT) * MSA.size();

    Pij.resize(seq_length);
    for (int i = 0; i < seq_length; i++)
        Pij[i].resize(seq_length - (i + 1));

    for (int j = 0; j < seq_length; j++)
    {
        for (auto &item : bestP[j])
        {
            int i = item.first;
            State state = item.second;
            // malikap: only divide by kT since pscore is already divided by n_seq
            value_type temp_prob_inside = state.alpha + state.beta - viterbi.alpha;
            if (temp_prob_inside > value_type(-9.91152))
            {
                value_type prob = Fast_Exp(temp_prob_inside);
                if (prob > value_type(1.0))
                    prob = value_type(1.0);
                if (prob < value_type(bpp_cutoff))
                    continue;
                Pij[i][j - (i + 1)] = prob;
            }
        }
    }

    if (bpp_file_name.size() > 0) {
        output_to_file(bpp_file_name, "a");
    }
}

void BeamCKYParser::threshknot() {
    vector<double> best_prob(seq_length, 0.0);
    vector<pair<int, int>> pairs;
    vector<pair<int, int>> pseudo_pairs1;
    vector<pair<int, int>> pseudo_pairs2;
    vector<pair<int, int>> pseudo_pairs3;
    set<int> visited;

    for(int i = 0; i < seq_length; i++) {
        for(int j = i + turn + 1; j < seq_length; j++) {
            double prob = Pij[i][j - (i + 1)];
            if(prob >= 0.3) {
                best_prob[i] = max(best_prob[i], prob);
                best_prob[j] = max(best_prob[j], prob);
            }
        }
    }

    for(int i = 0; i < seq_length; i++) {
        for(int j = i + turn + 1; j < seq_length; j++) {
            double prob = Pij[i][j - (i + 1)];
            if(prob >= 0.3 && prob == best_prob[i] && prob == best_prob[j]) {
                if (visited.find(i) != visited.end() || visited.find(j) != visited.end())
                    continue;
                
                pairs.push_back(make_pair(i, j));
                visited.insert(i);
                visited.insert(j);
            }
        }
    }

    // check for pseudoknots 1
    for (int i = 0; i < pairs.size(); i++) {
        for (int j = i + 1; j < pairs.size(); j++) {
            if (pairs[i].first < pairs[j].first && pairs[j].first < pairs[i].second && pairs[i].second < pairs[j].second) {
                pseudo_pairs1.push_back(make_pair(pairs[j].first, pairs[j].second));
            }
        }
    }

    // check for pseudoknots 2
    for(int i = 0; i < pseudo_pairs1.size(); i++) {
        for(int j = i + 1; j < pseudo_pairs1.size(); j++) {
            if(pseudo_pairs1[i].first < pseudo_pairs1[j].first && pseudo_pairs1[j].first < pseudo_pairs1[i].second && pseudo_pairs1[i].second < pseudo_pairs1[j].second) {
                pseudo_pairs2.push_back(make_pair(pseudo_pairs1[j].first, pseudo_pairs1[j].second));
            }
        }
    }

    // check for pseudoknots 3
    for(int i = 0; i < pseudo_pairs2.size(); i++) {
        for(int j = i + 1; j < pseudo_pairs2.size(); j++) {
            if(pseudo_pairs2[i].first < pseudo_pairs2[j].first && pseudo_pairs2[j].first < pseudo_pairs2[i].second && pseudo_pairs2[i].second < pseudo_pairs2[j].second) {
                pseudo_pairs3.push_back(make_pair(pseudo_pairs2[j].first, pseudo_pairs2[j].second));
            }
        }
    }

    // struc = ['.' for _ in range(seq_length)]
    string structure(seq_length, '.');

    // normal pairs
    for(int i = 0; i < pairs.size(); i++) {
        structure[pairs[i].first] = '(';
        structure[pairs[i].second] = ')';
    }

    // pseudoknots 1
    for(int i = 0; i < pseudo_pairs1.size(); i++) {
        structure[pseudo_pairs1[i].first] = '[';
        structure[pseudo_pairs1[i].second] = ']';
    }

    // pseudoknots 2
    for(int i = 0; i < pseudo_pairs2.size(); i++) {
        structure[pseudo_pairs2[i].first] = '{';
        structure[pseudo_pairs2[i].second] = '}';
    }

    // pseudoknots 3
    for(int i = 0; i < pseudo_pairs3.size(); i++) {
        structure[pseudo_pairs3[i].first] = '<';
        structure[pseudo_pairs3[i].second] = '>';
    }


    if (threshknot_file_name.size() > 0) {
        FILE *fptr = fopen(threshknot_file_name.c_str(), "w");
        if (fptr == NULL)
        {
            printf("Could not open file!\n");
            return;
        }
        fprintf(fptr, "%s\n\n", structure.c_str());
    } else {
        printf("\nThreshknot Structure:\n");
        printf("%s\n\n", structure.c_str());
    }
}



void BeamCKYParser::get_mea(double gamma)
{
    vector<vector<float>> OPT;
    OPT.resize(seq_length);
    for (int i = 0; i < seq_length; ++i)
        OPT[i].resize(seq_length);

    vector<vector<int>> back_pointer;
    back_pointer.resize(seq_length);
    for (int i = 0; i < seq_length; ++i)
        back_pointer[i].resize(seq_length);

    vector<vector<int>> paired;
    paired.resize(seq_length);

    vector<value_type> Q;
    for (int i = 0; i < seq_length; ++i)
        Q.push_back(value_type(1.0));

    for (int i = 0; i < seq_length; ++i) {
        for (int j = i + turn + 1; j < seq_length; ++j)
            if (Pij[i][j - (i + 1)] > 0) {
                double prob = Pij[i][j - (i + 1)];

                paired[i].push_back(j);
                Q[i] -= prob;
                Q[j] -= prob;
            }
    }

    for (int i = 0; i < seq_length; ++i)
        std::sort(paired[i].begin(), paired[i].end());

    for (int l = 0; l < seq_length; l++)
    {
        for (int i = 0; i < seq_length - l; i++)
        {
            int j = i + l;
            if (i == j)
            {
                OPT[i][j] = Q[i];
                back_pointer[i][j] = -1;
                continue;
            }
            OPT[i][j] = OPT[i][i] + OPT[i + 1][j];
            back_pointer[i][j] = -1;
            for (int k : paired[i])
            {
                if (k > j)
                    break;
                value_type temp_OPT_k1_j;
                if (k < j)
                    temp_OPT_k1_j = OPT[k + 1][j];
                else
                    temp_OPT_k1_j = value_type(0.);
                auto temp_score = 2 * gamma * Pij[i][k - (i + 1)] + OPT[i + 1][k - 1] + temp_OPT_k1_j;
                if (OPT[i][j] < temp_score)
                {
                    OPT[i][j] = temp_score;
                    back_pointer[i][j] = k;
                }
            }
        }
    }

    auto structure = backtrace(0, seq_length - 1, back_pointer);

    if (mea_file_name.size() > 0) {
        FILE *fptr = fopen(mea_file_name.c_str(), "a");
        if (fptr == NULL)
        {
            printf("Could not open MEA file!\n");
            return;
        }
        fprintf(fptr, "\nMEA Structure (gamma = %.2f):\n", gamma);
        fprintf(fptr, "%s\n\n", structure.c_str());
    } else {
        printf("\nMEA Structure (gamma = %.2f):\n", gamma);
        printf("%s\n\n", structure.c_str());
    }
}

string BeamCKYParser::backtrace(const int i, const int j, const vector<vector<int>> &back_pointer, int *bp_num, double *etp)
{
    if (i > j)
        return "";
    if (back_pointer[i][j] == -1)
    {
        if (i == j)
            return ".";
        else
            return "." + backtrace(i + 1, j, back_pointer, bp_num, etp);
    }
    else if (back_pointer[i][j] != 0)
    {
        int k = back_pointer[i][j];
        assert(k + 1 > 0 && k + 1 <= seq_length);
        string temp;
        if (k == j)
            temp = "";
        else
            temp = backtrace(k + 1, j, back_pointer, bp_num, etp);

        if (bp_num != NULL)
            (*bp_num)++;
        if (etp != NULL)
            (*etp) += Pij[i][k - (i + 1)];
            
        return "(" + backtrace(i + 1, k - 1, back_pointer, bp_num, etp) + ")" + temp;
    }
    else {
        return "." + backtrace(i + 1, j, back_pointer, bp_num, etp);
    }
    // assert(false);
    return "";
}



double BeamCKYParser::get_centroid(bool maximize_pmcc) {
    string best_structure = "";
    double max_pmcc = 0.0;
    double best_gamma = 0.0;

    if (maximize_pmcc) {
        double bpp_sum = 0.0;
        for (int i = 0; i < seq_length; ++i)
            for (int j = i + turn + 1; j < seq_length; ++j)
                bpp_sum += Pij[i][j - (i + 1)];

        vector<double> gamma_array;
        for(int i = 1; i<=10; i++){
            gamma_array.push_back(pow(2,i));
            if(i==2) {
                gamma_array.push_back(6);
            }
        }

        for(int i = gamma_array.size() - 1 ; i >= 0 ; i--){
            double gamma = gamma_array[i];
            double bpp_threshold = 1 / (gamma + 1);

            int num_of_bp = 0;
            double etp = 0;
            string structure = predict_centroid(gamma, &num_of_bp, &etp);

            double etn = seq_length * ((seq_length - 1) / 2) - num_of_bp - bpp_sum + etp;
            double efp = num_of_bp - etp;
            double efn = bpp_sum - etp;

            double pmcc = (etp * etn) - (efp * efn);
            pmcc /= sqrt((etp + efp) * (etp + efn) * (etn + efp) * (etn + efn));

            if (pmcc > max_pmcc){
                best_structure = structure;
                max_pmcc = pmcc;
                best_gamma = gamma;
            }
        }
    }
    else {
        best_gamma = mea_gamma;
        best_structure = predict_centroid(mea_gamma);
    }

    if (centroid_file_name.size() > 0) {
        FILE *fptr = fopen(centroid_file_name.c_str(), "w");
        if (fptr == NULL)
        {
            printf("[ERROR] Could not open Centroid file!\n");
            return -1;
        }
        if (maximize_pmcc)
            fprintf(fptr, "Centroid Structure (pmcc = %.4f, gamma = %.2f):\n", max_pmcc, best_gamma);
        else
            fprintf(fptr, "Centroid Structure (gamma = %.2f):\n", mea_gamma);
        fprintf(fptr, "%s\n\n", best_structure.c_str());
    } else {
        if (maximize_pmcc)
            printf("Centroid Structure (pmcc = %.4f, gamma = %.2f):\n", max_pmcc, best_gamma);
        else
            printf("Centroid Structure (gamma = %.2f):\n", mea_gamma);
        printf("%s\n\n", best_structure.c_str());
    }

    return best_gamma;
}

string BeamCKYParser::predict_centroid(double gamma, int* bp_num, double* etp) {
    vector<vector<float>> dp_matrix(seq_length, vector<float>(seq_length, 0.0));
    vector<vector<int>> back_pointer(seq_length, vector<int>(seq_length, 0));
    vector<vector<int>> paired_base(seq_length);

    double bpp_threshold = 1.0 / (gamma + 1);
    for (int i = 0; i < seq_length; ++i) {
        for (int j = i + turn + 1; j < seq_length; ++j) {
            if (Pij[i][j - (i + 1)] > bpp_threshold) {
                paired_base[i].push_back(j);
            } else {
                Pij[i][j - (i + 1)] = 0.0;
            }
        }
    }

    for (int i = 0; i < seq_length; ++i)
        sort(paired_base[i].begin(), paired_base[i].end());

    for (int l = turn + 1; l< seq_length; l++) {
        for (int i = 0; i < seq_length - l; i++) {	  
            int j = i + l;
            dp_matrix[i][j] = dp_matrix[i + 1][j];
            for (int k : paired_base[i]) {
                if (k > j) break;
                float score_kp1_j = 0.0;
                if (k < j){
                    score_kp1_j = dp_matrix[k+1][j];
                }
                float temp_score = (gamma + 1) * Pij[i][k - (i + 1)] - 1 + dp_matrix[i + 1][k - 1] + score_kp1_j;
                if (dp_matrix[i][j] < temp_score){
                    dp_matrix[i][j] = temp_score; 
                    back_pointer[i][j] = k;
                }
            }
        }
    }

    string structure = backtrace(0, seq_length - 1, back_pointer, bp_num, etp);
    return structure;
}


void BeamCKYParser::output_to_file(string file_name, const char *type)
{
    if (!file_name.empty())
    {
        FILE *fptr = fopen(file_name.c_str(), type);
        if (fptr == NULL)
        {
            printf("Could not open file!\n");
            return;
        }
        
        for (int i = 0; i < seq_length; i++)
            for (int j = i + turn + 1; j < seq_length; j++)
                if (Pij[i][j - (i + 1)] > 0)
                    fprintf(fptr, "%d %d %.4e\n", i + 1, j + 1, Pij[i][j - (i + 1)]);

        fprintf(fptr, "\n");
        fclose(fptr);
    }
    return;
}

// cat ../eval/Data/samples/sample07.fasta | ./bin/laf_mfe_vienna 100 1 1 1 -40 1 1 1