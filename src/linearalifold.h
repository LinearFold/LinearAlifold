
#ifndef FASTCKY_BEAMCKYPAR_H
#define FASTCKY_BEAMCKYPAR_H

#include <limits>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <random>

#include "utils/sparsehash/dense_hash_map"
#include "utils/energy_model.h"
#include "utils/log_space.h"

using namespace std;

#define MIN_CUBE_PRUNING_SIZE 20
#define kT 61.63207755

#define SINGLE_MAX_LEN 30


enum Manner {
    MANNER_NONE = 0,              // 0: empty
    MANNER_H,                     // 1: hairpin candidate
    MANNER_HAIRPIN,               // 2: hairpin
    MANNER_SINGLE,                // 3: single
    MANNER_HELIX,                 // 4: helix
    MANNER_MULTI,                 // 5: multi = ..M2. [30 restriction on the left and jump on the right]
    MANNER_MULTI_eq_MULTI_plus_U, // 6: multi = multi + U
    MANNER_P_eq_MULTI,            // 7: P = (multi)
    MANNER_M2_eq_M_plus_P,        // 8: M2 = M + P
    MANNER_M_eq_M2,               // 9: M = M2
    MANNER_M_eq_M_plus_U,         // 10: M = M + U
    MANNER_M_eq_P,                // 11: M = P
    /* MANNER_C_eq_U, */
    /* MANNER_C_eq_P, */
    MANNER_C_eq_C_plus_U, // 12: C = C + U
    MANNER_C_eq_C_plus_P, // 13: C = C + P
};


// bool cmp(const tuple<value_type, int, int> &a, const tuple<value_type, int, int> &b) {
//     if (get<0>(a) != get<0>(b))
//         return (get<0>(a) < get<0>(b));

//     else
//         return ((get<1>(a) - get<2>(a)) < (get<1>(b) - get<2>(b)));
// }

// A hash function used to hash a pair of any kind
// struct hash_pair {
//     template <class T1, class T2>
//     size_t operator()(const pair<T1, T2> &p) const {
//         auto hash1 = hash<T1>{}(p.first);
//         auto hash2 = hash<T2>{}(p.second);
//         return hash1 ^ hash2;
//     }
// };

struct pair_hash {
    inline std::size_t operator()(const std::pair<int, int>& v) const {
        std::size_t seed = 0;
        hash_combine(seed, v.first);
        hash_combine(seed, v.second);
        return seed;
    }

    // Generic hash_combine function based on Boost's hash_combine
    template <class T>
    inline void hash_combine(std::size_t& seed, const T& v) const {
        std::hash<T> hasher;
        seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
};


enum Type { // reverse topological order
    TYPE_C = 0,
    TYPE_M = 1,
    TYPE_M2 = 2,
    TYPE_P = 3,
    TYPE_MULTI = 4,
    TYPE_MAX = 5,
};
namespace std {
    template <>
    struct hash<Type> {
        std::size_t operator()(const Type& t) const {
            return std::hash<int>()(static_cast<int>(t));
        }
    };
}

struct State {
    value_type alpha; // inside score
    value_type beta;  // outside score

    Manner manner;
    Type type;
    pair<int, int> idx;

    union TraceInfo {
        int split;
        struct
        {
            int l1;
            int l2;
        } paddings;
    };

    TraceInfo trace;

    State() : manner(MANNER_NONE), alpha(VALUE_MIN), beta(VALUE_MIN){};
    State(value_type s, Manner m) : alpha(s), beta(VALUE_MIN), manner(m){};

    inline void logplus(value_type outside) {
        Fast_LogPlusEquals(this->beta, outside);
    }

    void set_attributes(int i, int j, Type type) {
        this->type = type;
        idx.first = i;
        idx.second = j;
    }

    void set(value_type score_, Manner manner_) {
        this->alpha = score_;
        manner = manner_;
    }

    void set(value_type score_, Manner manner_, int split_) {
        this->alpha = score_;
        manner = manner_;
        trace.split = split_;
    }

    void set(value_type score_, Manner manner_, int l1_, int l2_) {
        this->alpha = score_;
        manner = manner_;
        trace.paddings.l1 = l1_;
        trace.paddings.l2 = l2_;
    }
};

// lazy outside ------------------------------------------------------------------------------------------------------------------------------------------------
// typedef google::dense_hash_map<int, State> *mypointer;


struct StateKey {
    int i, j;
    Type type;
    StateKey(int i, int j, Type type) : i(i), j(j), type(type){};
    bool operator==(const StateKey &other) const {
        return (i == other.i && j == other.j && type == other.type);
    }
};

namespace std {
	template <>
	struct hash<StateKey> {
	    std::size_t operator()(const StateKey& key) const {
		const std::size_t prime = 31; // A prime number for hash combination
		size_t res = 0;

		// Hash individual fields
		auto hash_i = std::hash<int>()(key.i);
		auto hash_j = std::hash<int>()(key.j);
		auto hash_type = std::hash<Type>()(key.type); // Ensure Type is hashable

		// Combine the hashes, taking into account the relationship i < j
		res = hash_i ^ (hash_j << 1); // Shift j's hash to ensure distinctiveness given i < j
		res = res * prime + hash_type;

		return res;
	    }
	};
}

// unified hyperedge
struct HEdge {
    State *left = NULL, *right = NULL; // right=null <=> Edge
    value_type weight;
    HEdge(State *left = NULL, State *right = NULL, value_type weight = VALUE_MIN) : left(left), right(right), weight(weight){}; // default constructor needed

    void set(State *left, State *right, value_type weight) {
        this->left = left;
        this->right = right;
        this->weight = weight;
    }
};

// -------------------------------------------------------------------------------------------------------------------------------------------------------------

class BeamCKYParser {
  public:
    BeamCKYParser(
        int beam_size = 100,
        bool is_verbose = true,
        bool multi_approx = false,
        bool partition_mode = false,        
        bool use_lazy_outside = true,
        float pscore_threshold = -40,
        float pscore_beta = 1.0,
        float pscore_delta = 1.0,
        float mea_gamma = 2.0,
        float bpp_cutoff = 0.0,
        float threshknot_threshold = 0.3,
        int sampling_size = 0,
        string bpp_file_name = "",
        string mea_file_name = "",
        string threshknot_file_name = "",
        string centroid_file_name = ""
    );

    int beam;
    int seq_length;

    string bpp_file_name;
    string mea_file_name;
    string threshknot_file_name;
    string centroid_file_name;

    float bpp_cutoff = 0.0;
    float mea_gamma = 2.0;
    float pscore_beta = 1.0;
    float pscore_delta = 1.0;
    float pscore_threshold = -40;
    float threshknot_threshold = 0.3;
    int sampling_size = 0;
    
    bool partition_mode = false;
    bool use_lazy_outside = true;
    bool multi_approx = false;
    bool is_verbose = true;

    void parse_alifold(vector<string> &MSA, float **ribo, vector<vector<int>> &a2s_fast, vector<vector<int>> &s5_fast, vector<vector<int>> &s3_fast, vector<vector<int>> &SS_fast, vector<float> &smart_gap);

    void run_inside();
    void outside_alifold();
    void run_lazy_outside();
    void threshknot();

  private:
    void get_parentheses(char *result, string &seq);
    void backtrack(State *state, string &structure, bool best_only=true, default_random_engine *generator=NULL);
    vector<HEdge> get_incoming_hedges(State *state, double multiplier, bool best_only);
    vector<pair<int, int>> get_pairs(string &structure);

    vector<vector<double>> Pij;
    void output_to_file(string file_name, const char *type);
    void cal_PairProb(State &viterbi);

    std::unordered_map<StateKey, vector<HEdge>, hash<StateKey>> state_hedges_cache;

    string backtrace(const int i, const int j, const vector<vector<int>> &back_pointer, int *bp_num = NULL, double *etp = NULL);
    void get_mea(double gamma);

    // centroid prediction methods
    double get_centroid(bool maximize_pmcc=true);
    string predict_centroid(double gamma, int* bp_num = NULL, double* etp = NULL);

    void postprocess();


    double _make_pscores_ij(int i, int j);
    bool check_pairable_ij(int i, int j);
    double get_pscore(int i, int j);


    vector<google::dense_hash_map<int, State>> bestH, bestP, bestM2, bestMulti, bestM;
    vector<google::dense_hash_map<int, State>> *best_states[5]; // pointing to (bestC), bestM, bestM2, bestP, bestMulti

    // special hairpin stuff
    vector<string> seq_MSA_no_gap;
    vector<vector<int>> if_tetraloops_MSA;
    vector<vector<int>> if_hexaloops_MSA;
    vector<vector<int>> if_triloops_MSA;

    // malikap: add ribo, pscore, and MSA
    float **ribo;
    vector<vector<double>> pscore;
    vector<string> MSA;

    // multiplier
    double inv_n;
    double inv_kt;
    double inv_ktn;

    // malikap: add a2s_fast, s5_fast, s3_fast, SS_fast as class members
    vector<vector<int>> a2s_fast, s5_fast, s3_fast, SS_fast;

    // malikap: add smartgap as class member
    float smart_gap_threshold = 0.5;
    vector<float> smart_gap;

    // malikap: utility function to get next and prev positions in cache
    vector<vector<int>> next_pos;
    vector<google::dense_hash_map<int, int>> prev_pos;
    int get_cache_next_position(int i, int j);
    int get_cache_prev_position(int i, int j);

    // same as bestM, but ordered
    vector<vector<pair<value_type, int>>> sorted_bestM;

    // hzhang: sort keys in each beam to avoid randomness
    vector<pair<int, State>> keys;

    // hzhang: sort keys in each beam to avoid randomness
    void sort_keys(google::dense_hash_map<int, State> &map, vector<pair<int, State>> &sorted_keys);

    void sortM(value_type threshold,
               google::dense_hash_map<int, State> &beamstep,
               vector<pair<value_type, int>> &sorted_stepM);

    vector<State> bestC;

    void prepare(unsigned num_seq, unsigned seq_length);

    void update_if_better(State &state, value_type newscore, Manner manner) {
        if (state.alpha < newscore)
            state.set(newscore, manner);
    };

    void update_if_better(State &state, value_type newscore, Manner manner, int split) {
        if (state.alpha < newscore || state.manner == MANNER_NONE)
            state.set(newscore, manner, split);
    };

    void update_if_better(State &state, value_type newscore, Manner manner, int l1, int l2) {
        if (state.alpha < newscore || state.manner == MANNER_NONE)
            state.set(newscore, manner, l1, l2);
    };

    value_type beam_prune_mfe(google::dense_hash_map<int, State> &beamstep);
    value_type beam_prune_partition(google::dense_hash_map<int, State> &beamstep);
    value_type beam_prune(google::dense_hash_map<int, State> &beamstep);

    // vector to store the scores at each beam temporarily for beam pruning
    vector<pair<value_type, int>> scores;

    // lazy outside
    void backward_update(int i, int j, State &state, Type type, float global_threshold, int* pruned, int* saved);


    vector<int> *sortedP; // sorted by -i for P(i,j); used for M2=M+P backwards

    float inside_time;
    float deviation_threshold = 9.91152;     // inside + outside > global - ..
    // float deviation_threshold = numeric_limits<float>::lowest();
    // float deviation_threshold = numeric_limits<float>::max();
    // float deviation_threshold = 0;
};

#endif // FASTCKY_BEAMCKYPAR_H
