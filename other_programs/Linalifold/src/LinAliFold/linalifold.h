#ifndef LINALIFOLD_H
#define LINALIFOLD_H

#include <getopt.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <stack>
#include <limits>
#include <algorithm>
#include <unordered_map>
#include <cstring>
#include <cmath>
#include <sys/time.h>

#include "energy_par.h"
#include "intloops.h"

using namespace std;

struct hash_pair{ 
  template <class T1, class T2> 
  size_t operator()(const pair<T1, T2>& p) const{ 
    auto hash1 = hash<T1>{}(p.first); 
    auto hash2 = hash<T2>{}(p.second); 
    return hash1 ^ hash2; 
  } 
};

enum Manner {
  NONE = 0,
  HAIRPIN,
  INTERNAL,
  STEM,
  MULTI,
  MULTI_eq_MULTI_plus_UNPAIRED,
  S_eq_MULTI,
  MULTI2_eq_MULTI1_plus_STEM,
  MULTI1_eq_MULTI2,
  MULTI1_eq_MULTI1_plus_UNPAIRED, 
  MULTI1_eq_STEM, 
  OUTER_eq_OUTER_plus_UNPAIRED,
  OUTER_eq_OUTER_plus_STEM, 
};

struct State {
  double score;
  Manner type;
  int split;
  int l1;
  int l2;

  State(){
    score = MIN_VALUE;
    type = NONE;
  };

  void set(double i, Manner m){
    score = i; type = m;
  }

  void set(double i, Manner m, int s){
    score = i; type = m; split = s;;
  }

  void set(double i, Manner m, int a, int b){
    score = i; type = m; l1 = a; l2 = b;
  }
};

class LinAliFold{
  public:
  LinAliFold(){
    _beam_size = 100;
    _beta = 1.0;
    _delta  = 1.0;
    _conservation_score_threshold = NEG_INF;
    _input_file_name = "";
    _ribosum_flag = 0;
  }
  void Run();
  void SetParameters(int argc, char* argv[]);

private:
  void ReadData();
  void Parse();
  void Initialize();
  void TraceBack();
  double AllowBasePairs(int h, int i);
  int PruneBeam(unordered_map<int, State> &candidate_list);
  int QuickSelect(vector<pair<int, int> >& scores, int lower,int upper, int k);
  int QuickSelectPartition(vector<pair<int, int> >& scores, int lower, int upper);
  void UpdateIfBetter(State& s, double new_score, Manner m);
  void UpdateIfBetter(State &s, double new_score, Manner m, int split);
  void UpdateIfBetter(State &s, double new_score, Manner m, int l1, int l2);
  int HairpinEnergy(int i, int j, int nuc_i, int nuc_i_p1, int nuc_j_m1, int nuc_j);
  int LoopEnergy(int type, int type2, int h, int i, int p, int q, int seq_id);
  int DangleEnergy(int type, int a, int b, int seq_id);
  void SortMulti1(int threshold, unordered_map<int, State> &candidate, vector<pair<int, int> > &sorted_multi1);
  void GetRibosumMatrix();
  void MakeRibosumData();

  int _ribosum_flag;
  int _beam_size;
  int _num_of_seq;
  int _alignment_length;
  double _beta;
  double _delta;
  double _conservation_score_threshold;
  string _input_file_name;
  string _structure;
  vector<string> _gapped_seq_list;
  vector<string> _name_list;
  vector<vector<int> > _a2s_index; //alignment to sequences
  vector<vector<int> > _int_seq_list; //ungapped
  vector<vector<double> > _allowed_base_pair;
  unordered_map<pair<int, int>, vector<vector <float> > ,hash_pair> _all_ribosum_data;
  vector<vector <float> > _ribosum_matrix;

  vector<unordered_map<int, State> > _stem_best;
  vector<unordered_map<int, State> > _multi2_best;
  vector<unordered_map<int, State> > _multi1_best;
  vector<unordered_map<int, State> > _multi_best;
  vector<vector<pair<int, int> > > _sorted_multi1_best;
  vector<State> _outer_best;
};

#endif
