#ifndef LINEAR_PARTITION_H
#define LINEAR_PARTITION_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <limits>
#include "energy_par.h"
#include "utils.h"

using namespace std;

class LinearPartition{
 public:
  LinearPartition(const vector<int>& v, const int b, const float t){
    _int_seq = v;
    _beam_size = b;
    _bpp_threshold = t;
    _utils.Set(_beam_size);
  }
  void Run(unordered_map<pair<int,int>, float, hash_pair>& bpp_matrix);
  

private:
  void Initialize();
  void Inside();
  void Outside();
  void CalcBPP(unordered_map<pair<int,int>, float, hash_pair>& bpp_matrix);

  void Stem(int i, int flag);
  void Multi(int i, int flag);
  void Multi1(int i, int flag);
  void Multi2(int i, int flag);

  Utils _utils;
  int _beam_size;
  int _seq_length;
  float _beta_pf;
  float _bpp_threshold;
  vector<int> _int_seq;
  vector<vector<int> > _next_pair;
  
  vector<unordered_map<int, State> > _hairpin_sum;
  vector<unordered_map<int, State> > _stem_sum;
  vector<unordered_map<int, State> > _multi2_sum;
  vector<unordered_map<int, State> > _multi1_sum;
  vector<unordered_map<int, State> > _multi_sum;
  vector<State> _outer_sum;
};

#endif
