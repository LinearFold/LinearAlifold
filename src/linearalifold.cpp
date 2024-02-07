
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <stack>
#include <string>
#include <sys/time.h>
#include <tuple>
#include <unordered_map>

#include "linearalifold.h"
#include "utils/ribo.h"
#include "utils/utility.h"
#include "utils/sparsehash/dense_hash_map"
// #include "log_space.h"

// #define SPECIAL_HP

using namespace std;


bool comparefunc(std::pair<int, State> a, std::pair<int, State> b) {
    return a.first > b.first;
}

void BeamCKYParser::sort_keys(google::dense_hash_map<int, State> &map, std::vector<std::pair<int, State>> &sorted_keys) {
    sorted_keys.clear();

    for (auto &kv : map) {
        sorted_keys.push_back(kv);
    }

    sort(sorted_keys.begin(), sorted_keys.end(), comparefunc);
}


void BeamCKYParser::get_parentheses(char *result, string &seq) {
    memset(result, '.', seq_length);
    result[seq_length] = 0;

    stack<tuple<int, int, State>> stk;
    stk.push(make_tuple(0, seq_length - 1, bestC[seq_length - 1]));

    if (is_verbose) {
        printf(">verbose\n");
    }
    // verbose stuff
    vector<pair<int, int>> multi_todo;
    unordered_map<int, int> mbp; // multi bp
    double total_energy = .0;
    double external_energy = .0;

    while (!stk.empty()) {
        tuple<int, int, State> top = stk.top();
        int i = get<0>(top), j = get<1>(top);
        State &state = get<2>(top);
        stk.pop();

        int k, p, q;

        switch (state.manner) {
        case MANNER_H:
            // this state should not be traced
            break;
        case MANNER_HAIRPIN: {
            result[i] = '(';
            result[j] = ')';
        } break;
        case MANNER_SINGLE: {
            result[i] = '(';
            result[j] = ')';
            p = i + state.trace.paddings.l1;
            q = j - state.trace.paddings.l2;
            stk.push(make_tuple(p, q, bestP[q][p]));
        } break;
        case MANNER_HELIX: {
            result[i] = '(';
            result[j] = ')';
            stk.push(make_tuple(i + 1, j - 1, bestP[j - 1][i + 1]));
        } break;
        case MANNER_MULTI:
            p = i + state.trace.paddings.l1;
            q = j - state.trace.paddings.l2;
            stk.push(make_tuple(p, q, bestM2[q][p]));
            break;
        case MANNER_MULTI_eq_MULTI_plus_U:
            p = i + state.trace.paddings.l1;
            q = j - state.trace.paddings.l2;
            stk.push(make_tuple(p, q, bestM2[q][p]));
            break;
        case MANNER_P_eq_MULTI:
            result[i] = '(';
            result[j] = ')';
            stk.push(make_tuple(i, j, bestMulti[j][i]));
            break;
        case MANNER_M2_eq_M_plus_P:
            k = state.trace.split;
            stk.push(make_tuple(i, k, bestM[k][i]));
            stk.push(make_tuple(k + 1, j, bestP[j][k + 1]));
            break;
        case MANNER_M_eq_M2:
            stk.push(make_tuple(i, j, bestM2[j][i]));
            break;
        case MANNER_M_eq_M_plus_U:
            stk.push(make_tuple(i, j - 1, bestM[j - 1][i]));
            break;
        case MANNER_M_eq_P:
            stk.push(make_tuple(i, j, bestP[j][i]));
            break;
        case MANNER_C_eq_C_plus_U:
            k = j - 1;
            if (k != -1)
                stk.push(make_tuple(0, k, bestC[k]));
            break;
        case MANNER_C_eq_C_plus_P: {
            k = state.trace.split;
            if (k != -1) {
                stk.push(make_tuple(0, k, bestC[k]));
                stk.push(make_tuple(k + 1, j, bestP[j][k + 1]));
            } else {
                stk.push(make_tuple(i, j, bestP[j][i]));
            }
        } break;
        default: // MANNER_NONE or other cases
            printf("wrong manner at %d, %d: manner %d\n", i, j, state.manner);
            fflush(stdout);
            assert(false);
        }
    }

    return;
}

unsigned long quickselect_partition(vector<pair<value_type, int>> &scores, unsigned long lower, unsigned long upper) {
    value_type pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot)
            ++lower;
        while (scores[upper].first > pivot)
            --upper;
        if (scores[lower].first == scores[upper].first)
            ++lower;
        else if (lower < upper)
            swap(scores[lower], scores[upper]);
    }
    return upper;
}

// in-place quick-select
value_type quickselect(vector<pair<value_type, int>> &scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if (lower == upper)
        return scores[lower].first;
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k)
        return scores[split].first;
    else if (k < length)
        return quickselect(scores, lower, split - 1, k);
    else
        return quickselect(scores, split + 1, upper, k - length);
}

value_type BeamCKYParser::beam_prune_partition(google::dense_hash_map<int, State> &beamstep) {
    scores.clear();
    for (auto &item : beamstep) {
        int i = item.first;
        State &cand = item.second;
        int k = i - 1;
        value_type newalpha = (k >= 0 ? bestC[k].alpha : value_type(0.0)) + cand.alpha;
        scores.push_back(make_pair(newalpha, i));
    }
    if (scores.size() <= beam)
        return VALUE_MIN;
    value_type threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold)
            beamstep.erase(p.second);
    }

    return threshold;
}

value_type BeamCKYParser::beam_prune_mfe(google::dense_hash_map<int, State> &beamstep) {
    scores.clear();
    for (auto &item : beamstep) {
        int i = item.first;
        State &cand = item.second;
        int k = i - 1;
        value_type newscore;
        // lisiz: for _V, avoid -inf-int=+inf
        if ((k >= 0) && (bestC[k].alpha == VALUE_MIN))
            newscore = VALUE_MIN;
        else
            newscore = (k >= 0 ? bestC[k].alpha : 0) + cand.alpha;
        scores.push_back(make_pair(newscore, i));
    }
    if (scores.size() <= beam)
        return VALUE_MIN;
    value_type threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold)
            beamstep.erase(p.second);
    }

    return threshold;
}

value_type BeamCKYParser::beam_prune(google::dense_hash_map<int, State> &beamstep) {
    return partition_mode ? beam_prune_partition(beamstep) : beam_prune_mfe(beamstep);
}

void BeamCKYParser::sortM(value_type threshold,
                          google::dense_hash_map<int, State> &beamstep,
                          std::vector<std::pair<value_type, int>> &sorted_stepM) {
    sorted_stepM.clear();
    if (threshold == VALUE_MIN) {
        // no beam pruning before, so scores vector not usable
        for (auto &item : beamstep) {
            int i = item.first;
            State &cand = item.second;
            int k = i - 1;
            value_type newscore;
            // lisiz: constraints may cause all VALUE_MIN, sorting has no use
            if ((k >= 0) && (bestC[k].alpha == VALUE_MIN))
                newscore = cand.alpha;
            else
                newscore = (k >= 0 ? bestC[k].alpha : 0) + cand.alpha;
            sorted_stepM.push_back(make_pair(newscore, i));
        }
    } else {
        for (auto &p : scores) {
            if (p.first >= threshold)
                sorted_stepM.push_back(p);
        }
    }

    sort(sorted_stepM.begin(), sorted_stepM.end(), std::greater<pair<value_type, int>>());
}


double BeamCKYParser::_make_pscores_ij(int i, int j) // malikap: fixed this function
{
    auto n_seq = SS_fast[i].size();

    int pfreq[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    int type;

    for (int s = 0; s < n_seq; s++) {
        if (SS_fast[i][s] == 5 && SS_fast[j][s] == 5) {
            type = 7;
        } else {
            type = get_AUCG_pair_type_CONTRAfold_indices(SS_fast[i][s], SS_fast[j][s]);
        }

        pfreq[type]++;
    }
    if (pfreq[0] * 2 + pfreq[7] > n_seq) {
        return NONE;
    }

    double score = 0.;

    for (int k = 1; k <= 6; k++) { /* ignore pairtype 7 (gap-gap) */
        score += (double)pfreq[k] * (double)(pfreq[k] - 1) / 2 * (double)ribo[k][k];
        for (int l = k + 1; l <= 6; l++) {
            score += (double)pfreq[k] * (double)pfreq[l] * (double)ribo[k][l];
        }
    }

    return (pscore_beta * (100 * score * inv_n - pscore_delta * 100 * (pfreq[0] + pfreq[7] * 0.25))) * (partition_mode ? inv_ktn : inv_n);
}


double BeamCKYParser::get_pscore(int i, int j) {
    if (j < i)
        return get_pscore(j, i);
    if (i < 0 || j >= seq_length || i == j)
        return VALUE_MIN;
    
    int j_adj = j - (i + 1);
    if (pscore[i][j_adj] == VALUE_MIN) {
        pscore[i][j_adj] = _make_pscores_ij(i, j);
    }

    return pscore[i][j_adj];
}

bool BeamCKYParser::check_pairable_ij(int i, int j) {
    return get_pscore(i, j) >= pscore_threshold;
}

vector<pair<int, int>> BeamCKYParser::get_pairs(string &structure) {

    vector<int> stack;
    vector<pair<int, int>> results;
    for (int i = 0; i < structure.size(); i++) {

        if (structure[i] == '(') {
            stack.push_back(i);
        }

        else if (structure[i] == ')') {
            auto left = stack.back();
            stack.pop_back();
            results.push_back(make_pair(left, i));
        }
    }

    return results;
}


int BeamCKYParser::get_cache_prev_position(int l, int r) {
    if (l < 0 || l >= seq_length || r < 0 || r >= seq_length)
        return -1;

    if (prev_pos[l].find(r) == prev_pos[l].end()) {
        prev_pos[l][r] = -1;
        for (int x = r - 1; x > l; x--) {
            if (check_pairable_ij(l, x)) {
                prev_pos[l][r] = x;
                break;
            }
        }
    }

    return prev_pos[l][r];
}

int BeamCKYParser::get_cache_next_position(int l, int r) {
    if (l < 0 || l >= seq_length || r < 0 || r >= seq_length) {
        return -1;
    }

    r -= (l + 1);
    if (next_pos[l][r] == 0) {
        next_pos[l][r] = -1;
        for (int x = l + r + 2; x < seq_length; x++) {  
            if (check_pairable_ij(l, x)) {
                next_pos[l][r] = x;
                break;
            }
        }
    }

    return next_pos[l][r];    
}

void BeamCKYParser::prepare(unsigned num_seq, unsigned seq_length) {

    BeamCKYParser::seq_length = seq_length;

    google::dense_hash_map<int, State> default_state_map;
    default_state_map.set_empty_key(-1);
    default_state_map.set_deleted_key(-2);

    bestH = vector<google::dense_hash_map<int, State>>(seq_length, default_state_map);
    bestM = vector<google::dense_hash_map<int, State>>(seq_length, default_state_map);
    bestM2 = vector<google::dense_hash_map<int, State>>(seq_length, default_state_map);
    bestMulti = vector<google::dense_hash_map<int, State>>(seq_length, default_state_map);
    bestP = vector<google::dense_hash_map<int, State>>(seq_length, default_state_map);
    bestC = vector<State>(seq_length);

#ifdef is_cube_pruning
    sorted_bestM.clear();
    sorted_bestM.resize(seq_length);
#endif

    scores.reserve(seq_length);

    best_states[TYPE_M] = &bestM;
    best_states[TYPE_M2] = &bestM2;
    best_states[TYPE_P] = &bestP;
    best_states[TYPE_MULTI] = &bestMulti;

    sortedP = new vector<int>[seq_length]; // malikap: for M2 = M + P backwards (outside)

    inv_n = 1.0 / num_seq;
    inv_kt = 1.0 / (double)kT;
    inv_ktn = 1.0 / (num_seq * (double)kT);
    
    // malikap: initialize next_position and prev_position 2D vectors
    // next_pos.resize(seq_length);
    // prev_pos.resize(seq_length); 
    // for (int i = 0; i < seq_length; i++) {
    //     next_pos[i].resize(seq_length, 0);
    //     prev_pos[i].resize(seq_length, 0);
    // }

    // initialize next_pos
    next_pos.resize(seq_length);
    for (int i = 0; i < seq_length; ++i) {
        next_pos[i].resize(seq_length - (i + 1));
    }

    // initialize prev_pos
    google::dense_hash_map<int, int> default_prev_pos_map;
    default_prev_pos_map.set_empty_key(-1);      // we'll only use non-negative keys
    default_prev_pos_map.set_deleted_key(-2);    // and -2 is for deleted key
    prev_pos = vector<google::dense_hash_map<int, int>>(seq_length, default_prev_pos_map);
    // next_pos = vector<google::dense_hash_map<int, int>>(seq_length, default_prev_pos_map);

    // initialize pscore matrix
    // google::dense_hash_map<int, double> default_pscore_map;
    // default_pscore_map.set_empty_key(-1);      // we'll only use non-negative keys
    // default_pscore_map.set_deleted_key(-2);    // and -2 is for deleted key
    // pscore = vector<google::dense_hash_map<int, double>>(seq_length, default_pscore_map);
    pscore.resize(seq_length);
    for(int i = 0; i < seq_length; ++i) {
        pscore[i].resize(seq_length - (i + 1), VALUE_MIN);
    }

    if (partition_mode) {
        pscore_threshold *= inv_kt;
    }
}

void BeamCKYParser::postprocess() {
    struct timeval pp_start_time, pp_end_time;
    gettimeofday(&pp_start_time, NULL);
    
    // clear memory
    bestH = vector<google::dense_hash_map<int, State>>();
    bestP = vector<google::dense_hash_map<int, State>>();
    bestM2 = vector<google::dense_hash_map<int, State>>();
    bestM = vector<google::dense_hash_map<int, State>>();
    bestC = vector<State>();
    bestMulti = vector<google::dense_hash_map<int, State>>();
    scores = vector<pair<value_type, int>>();

    next_pos = vector<vector<int>>();
    prev_pos = vector<google::dense_hash_map<int, int>>();
    pscore = vector<vector<double>>();

    a2s_fast = vector<vector<int>>();
    s5_fast = vector<vector<int>>();
    s3_fast = vector<vector<int>>();
    SS_fast = vector<vector<int>>();
    smart_gap = vector<float>();

    // delete sortedP
    delete[] sortedP;
    for (int i = 0; i < 7; ++i) {
        free(ribo[i]);  // Free each float array allocated and pointed to by ribo[i]
    }
    free(ribo);  // Finally, free the array of pointers

    gettimeofday(&pp_end_time, NULL);
    double pp_elapsed_time = pp_end_time.tv_sec - pp_start_time.tv_sec + (pp_end_time.tv_usec - pp_start_time.tv_usec) / 1000000.0;

    if (is_verbose) {
        printf("\nPostprocessing time: %.3f seconds\n", pp_elapsed_time);
    }
}

void BeamCKYParser::parse_alifold(std::vector<std::string> &MSA, float **ribo, vector<vector<int>> &a2s_fast, vector<vector<int>> &s5_fast, vector<vector<int>> &s3_fast, vector<vector<int>> &SS_fast, vector<float> &smart_gap) {
    // number of states
    unsigned long nos_H = 0, nos_P = 0, nos_M2 = 0,
                  nos_M = 0, nos_C = 0, nos_Multi = 0;


    BeamCKYParser::MSA = MSA;
    BeamCKYParser::ribo = ribo;
    BeamCKYParser::SS_fast = SS_fast;
    BeamCKYParser::a2s_fast = a2s_fast;
    BeamCKYParser::s5_fast = s5_fast;
    BeamCKYParser::s3_fast = s3_fast;
    BeamCKYParser::smart_gap = smart_gap;

    prepare(MSA.size(), MSA[0].length());

    if (sampling_size > 0) {
        partition_mode = true;
    }
    
    if (partition_mode) {
        run_inside();

        if (sampling_size == 0) {
            use_lazy_outside ? run_lazy_outside() : outside_alifold();
            cal_PairProb(bestC[seq_length - 1]);
            postprocess();
            PairProb_MEA(MSA[0]);
            threshknot();
        }
        else {
            printf("Running sampling mode\n");
            srand(std::chrono::system_clock::now().time_since_epoch().count());
            bestC[seq_length - 1].set_attributes(0, seq_length - 1, TYPE_C);
            default_random_engine generator(rand());

            for (int i = 0; i < sampling_size; i++) {
                // generate a sample structure
                string structure(seq_length, '.');
                backtrack(&bestC[seq_length - 1], structure, false, &generator);
                printf("%s\n", structure.c_str());
            }
            postprocess();
        }
    }
    else {
        run_inside();
        // string structure(seq_length, '.');
        // bestC[seq_length - 1].set_attributes(0, seq_length - 1, TYPE_C);
        // backtrack(&bestC[seq_length - 1], structure);
        // printf("\nNew Backtrack Structure:\n");
        // printf("%s\n", structure.c_str());
        postprocess();
    }
}


BeamCKYParser::BeamCKYParser(
    int beam_size,      
    bool verbose,
    bool multi_approx,
    bool partition_mode,
    bool use_lazy_outside,
    float pscore_threshold,
    float pscore_beta,
    float pscore_delta,
    float mea_gamma,
    float bpp_cutoff,
    float threshknot_threshold,
    int sampling_size,
    string bpp_file_name,
    string mea_file_name,
    string threshknot_file_name
    ): 
    beam(beam_size),
    is_verbose(verbose),
    multi_approx(multi_approx),
    partition_mode(partition_mode),
    use_lazy_outside(use_lazy_outside),
    pscore_threshold(pscore_threshold),
    pscore_beta(pscore_beta),
    pscore_delta(pscore_delta),
    mea_gamma(mea_gamma),
    bpp_cutoff(bpp_cutoff),
    threshknot_threshold(threshknot_threshold),
    sampling_size(sampling_size),
    bpp_file_name(bpp_file_name),
    mea_file_name(mea_file_name),
    threshknot_file_name(threshknot_file_name)
    {
        initialize();
    }

int main(int argc, char **argv) {
    struct timeval parse_alifold_starttime, parse_alifold_endtime;

    gettimeofday(&parse_alifold_starttime, NULL);

    // general parameters
    int beamsize;
    int energy_model;
    bool is_verbose;
    bool multi_approx;
    bool use_lazy_outside;
    int sampling_size;

    // output parameters
    string bpp_file_name;
    string mea_file_name;
    string threshknot_file_name;
    
    // pairability score parameters
    float pscore_beta;
    float pscore_delta;
    float pscore_threshold;
    float threshknot_threshold;
    
    // partition mode parameter, if true, then partition function is calculated
    bool partition_mode;    

    if (argc >= 1) {
        beamsize = atoi(argv[1]);
        is_verbose = atoi(argv[2]) == 1;
        energy_model = atoi(argv[3]);
        multi_approx = atoi(argv[4]) == 1;
        partition_mode = atoi(argv[5]) == 1;

        pscore_threshold = atoi(argv[6]);
        pscore_beta = atof(argv[7]);
        pscore_delta = atof(argv[8]);
        
        use_lazy_outside = atoi(argv[9]) == 1;
        threshknot_threshold = atof(argv[10]);
        sampling_size = atoi(argv[11]);
        
        bpp_file_name = argv[12];
        mea_file_name = argv[13];
        threshknot_file_name = argv[14];
    }

    std::vector<std::string> MSA;

    for (string seq; getline(cin, seq);) {
        if (seq.length() == 0)
            continue;

        if (seq[0] == ';' || seq[0] == '>') {
            continue;
        }

        // convert to uppercase
        transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

        // convert T to U
        replace(seq.begin(), seq.end(), 'T', 'U');

        MSA.push_back(seq);
    }

    auto n_seq = MSA.size();
    auto MSA_seq_length = MSA[0].size();

    auto ribo = get_ribosum(MSA, n_seq, MSA_seq_length);
    vector<float> smart_gap;
    vector<vector<int>> a2s_fast, s5_fast, s3_fast, SS_fast;
    a2s_prepare_is(MSA, n_seq, MSA_seq_length, a2s_fast, s5_fast, s3_fast, SS_fast, smart_gap);

    BeamCKYParser parser(beamsize, is_verbose, multi_approx, partition_mode, use_lazy_outside, pscore_threshold, pscore_beta, pscore_delta, 3.0, \
                        std::numeric_limits<float>::min(), threshknot_threshold, sampling_size, bpp_file_name, mea_file_name, threshknot_file_name);
    parser.parse_alifold(MSA, ribo, a2s_fast, s5_fast, s3_fast, SS_fast, smart_gap);

    if (is_verbose) {
        printf("\nBeam Size: %d\n", parser.beam);
        printf("Energy Model: %d\n", energy_model);
        printf("Multi Approx: %d\n", parser.multi_approx);
        printf("Partition Mode: %d\n", parser.partition_mode);
        printf("Sampling Mode: %d\n", parser.sampling_size > 0);
        printf("Use Lazy Outside: %d\n", parser.use_lazy_outside);
        
        printf("\nPairability-Score Threshold: %.2f\n", parser.pscore_threshold);
        printf("Pairability-Score Beta: %.2f\n", parser.pscore_beta);
        printf("Pairability-Score Delta: %.2f\n", parser.pscore_delta);
        printf("Threshknot Threshold: %.2f\n", parser.threshknot_threshold);
        printf("Sampling Size: %d\n", parser.sampling_size);

        printf("\nBPP File Name: %s\n", parser.bpp_file_name.size() > 0 ? parser.bpp_file_name.c_str() : "None");
        printf("MEA File Name: %s\n", parser.mea_file_name.size() > 0 ? parser.mea_file_name.c_str() : "None");
        printf("Threshknot File Name: %s\n", parser.threshknot_file_name.size() > 0 ? parser.threshknot_file_name.c_str() : "None");
    }
    return 0;
}
