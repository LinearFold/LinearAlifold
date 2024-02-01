#ifndef UTILS_H
#define UTILS_H

#include <unordered_map>
#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>
#include <cstring>
#include "intloops.h"
#include "energy_par.h"

using namespace std;

struct hash_pair{ 
  template <class T1, class T2> 
  size_t operator()(const pair<T1, T2>& p) const{ 
    auto hash1 = hash<T1>{}(p.first); 
    auto hash2 = hash<T2>{}(p.second); 
    return hash1 ^ hash2; 
  } 
};

struct State {
  float alpha;
  float beta;
  
  State(){
    alpha = MIN_VALUE;
    beta = MIN_VALUE;
  };
};

class Utils{
 public:
  Utils(){}
  int HairpinEnergy(int i, int j, int nuc_i, int nuc_i_p1, int nuc_j_m1, int nuc_j);
  int DangleEnergy(int type, int a, int b, int len, vector<int>& int_seq);
  int LoopEnergy(int type, int type2, int h, int i, int p,int q, vector<int>& int_seq);
  void Set(int b);
  void FastLogPlusEquals(float &x, float y);
  float FastExp(float x);
  void PruneBeam(unordered_map<int, State> &candidate_list, vector<State>& outer_sum);

 private:
  int _beam_size;
  
  float FastLogExpPlusOne(float x);
  int QuickSelectPartition(vector<pair<float, int> >& scores, int lower, int upper);
  float QuickSelect(vector<pair<float, int> >& scores, int lower,int upper, int k);
};

#endif
