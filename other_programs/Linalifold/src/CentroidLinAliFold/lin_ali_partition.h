#ifndef LIN_ALI_PARTITION_H
#define LIN_ALI_PARTITION_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <limits>

#include "energy_par.h"
#include "utils.h"

using namespace std;

class LinAliPartition{
 public:
  LinAliPartition(vector<vector<int> >& v, int b, int r, double e, double d, float t, vector<vector<int> >& a2s, vector<string> &g){
    _int_seq_list = v;
    _beam_size = b;
    _ribosum_flag = r;
    _beta = e;
    _delta = d;
    _bpp_threshold = t;
    _a2s_index = a2s;
    _gapped_seq_list = g;
    _cons_score_threshold = NEG_INF;
    
    _alignment_length = _gapped_seq_list[0].size();
    _num_of_seq = _int_seq_list.size();
    _utils.Set(_beam_size);
  }
  void Run(unordered_map<pair<int,int>, float, hash_pair>& bpp_matrix);
  

private:
  void Initialize();
  void Inside();
  void Outside();
  void CalcBPP(unordered_map<pair<int,int>, float, hash_pair>& bpp_matrix);
  
  double AllowBasePairs(int h, int i);
  void GetRibosumMatrix();
  void MakeRibosumData();
  
  void Stem(int i, int flag);
  void Multi(int i, int flag);
  void Multi1(int i, int flag);
  void Multi2(int i, int flag);
  void Clear();
  
  Utils _utils;
  int _ribosum_flag;
  int _beam_size;
  int _alignment_length;
  int _num_of_seq;
  float _beta_pf;
  float _bpp_threshold;
  double _beta;
  double _delta;
  double _cons_score_threshold;
  vector<string> _gapped_seq_list;
  vector<vector<int> > _int_seq_list;
  vector<vector<double> > _allowed_base_pair;
  vector<vector<int> > _a2s_index;
  unordered_map<pair<int, int>, vector<vector <float> > ,hash_pair> _all_ribosum_data;
  vector<vector <float> > _ribosum_matrix;
  
  vector<unordered_map<int, State> > _stem_sum;
  vector<unordered_map<int, State> > _multi2_sum;
  vector<unordered_map<int, State> > _multi1_sum;
  vector<unordered_map<int, State> > _multi_sum;
  vector<State> _outer_sum;
};

#endif
