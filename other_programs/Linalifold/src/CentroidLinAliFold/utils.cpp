#include "utils.h"
#include <iostream>

void Utils::Set(int b){
  _beam_size = b;
}

int Utils::HairpinEnergy(int i, int j, int nuc_i, int nuc_i_p1, int nuc_j_m1, int nuc_j){
  int size = j-i-1;
  
  int type = BP_pair[nuc_i][nuc_j];
  int energy = size <= 30 ? hairpin37[size] : hairpin37[30] + (int)lxc37*log( size/30.);
  
  if(size != 3){
    energy += mismatchH37[type][nuc_i_p1][nuc_j_m1];
  }else{
    if(type > 2 || type == 0 ){energy += TerminalAU;}
  }
  
  return energy;
}

int Utils::DangleEnergy(int type, int a, int b, int len, vector<int>& int_seq){
  int x = 0;
 
  if (a >= 0 && int_seq[a] > 0) x += dangle5_37[type][int_seq[a]];
  if (b < len && int_seq[b] > 0) x += dangle3_37[type][int_seq[b]];
  if(type>2 || type == 0){
    x += TerminalAU;
  }
  
  return(x);
}

int Utils::LoopEnergy(int type, int type2, int h, int i, int p,int q, vector<int>& int_seq){
  int z=0;
  int u1 = h-p-1;
  int u2 = q-i-1;
  
  if ((u1==0) && (u2==0)){
    z = stack37[type][type2];
  }else{
    if ((u1==0)||(u2==0)) {
      int u;
      u = u1 == 0 ? u2 : u1;
      z = u <=30 ? bulge37[u] : bulge37[30]+lxc37*log( u/30.);
      
      if (u == 1){
	z += stack37[type][type2];
      }else {
	if (type>2 || type == 0 ){ z += TerminalAU;}
	if (type2>2  || type2 == 0 ){ z += TerminalAU;}
      }
    }else{
      if(type == 0){ type = 7;}
      if(type2 == 0){ type2 = 7;}
      
      if (u1+u2==2) {
	z = int11_37[type][type2][int_seq[p+1]][int_seq[q-1]];
      }else if ((u1==1) && (u2==2)){
	z = int21_37[type][type2][int_seq[p+1]][int_seq[i+1]][int_seq[q-1]];
      }else if ((u1==2) && (u2==1)){
	z = int21_37[type2][type][int_seq[i+1]][int_seq[p+1]][int_seq[h-1]];
      }else if ((u1==2) && (u2==2)){
	z = int22_37[type][type2][int_seq[p+1]][int_seq[h-1]][int_seq[i+1]][int_seq[q-1]];
      }else{
	if(type == 7){type = 0;}
	if(type2 == 7){type2 = 0;}
	z = internal_loop37[u1+u2]+mismatchI37[type][int_seq[p+1]][int_seq[q-1]]+mismatchI37[type2][int_seq[i+1]][int_seq[h-1]];
	z += min(MAX_NINIO, abs(u1-u2)*F_ninio37);
      }
    }
  }
  return z;
}

// log space: borrowed from CONTRAfold

float Utils::FastLogExpPlusOne(float x){
  if (x < float(3.3792499610)){
    if (x < float(1.6320158198)){
      if (x < float(0.6615367791)){
	return ((float(-0.0065591595)*x+float(0.1276442762))*x+float(0.4996554598))*x+float(0.6931542306);
      }
      return ((float(-0.0155157557)*x+float(0.1446775699))*x+float(0.4882939746))*x+float(0.6958092989);
    }
    if (x < float(2.4912588184)){
      return ((float(-0.0128909247)*x+float(0.1301028251))*x+float(0.5150398748))*x+float(0.6795585882);
    }
    return ((float(-0.0072142647)*x+float(0.0877540853))*x+float(0.6208708362))*x+float(0.5909675829);
  }
  if (x < float(5.7890710412)){
    if (x < float(4.4261691294)){
      return ((float(-0.0031455354)*x+float(0.0467229449))*x+float(0.7592532310))*x+float(0.4348794399);
    }
    return ((float(-0.0010110698)*x+float(0.0185943421))*x+float(0.8831730747))*x+float(0.2523695427);
  }
  if (x < float(7.8162726752)){
    return ((float(-0.0001962780)*x+float(0.0046084408))*x+float(0.9634431978))*x+float(0.0983148903);
  }
  return ((float(-0.0000113994)*x+float(0.0003734731))*x+float(0.9959107193))*x+float(0.0149855051);
}

void Utils::FastLogPlusEquals(float &x, float y){
  if (x < y){
    std::swap(x, y);
  }
  if (y > float(NEG_INF/2) && x-y < float(11.8624794162)){
    x = FastLogExpPlusOne(x-y) + y;
  }

  //exact logsumexp
  //float temp = x > y ? x + log(exp(y-x) + 1.0) : y + log(exp(x-y) + 1.0) ;
  //x = temp;
}

float Utils::FastExp(float x){    
  if (x < float(-2.4915033807)){
    if (x < float(-5.8622823336)){
      if (x < float(-9.91152)){
	return float(0);
      }
      return ((float(0.0000803850)*x+float(0.0021627428))*x+float(0.0194708555))*x+float(0.0588080014);
    }
    if (x < float(-3.8396630909)){
      return ((float(0.0013889414)*x+float(0.0244676474))*x+float(0.1471290604))*x+float(0.3042757740);
    }
    return ((float(0.0072335607)*x+float(0.0906002677))*x+float(0.3983111356))*x+float(0.6245959221);
  }
  
  if (x < float(-0.6725053211)){
    if (x < float(-1.4805375919)){
      return ((float(0.0232410351)*x+float(0.2085645908))*x+float(0.6906367911))*x+float(0.8682322329);
    }
    return ((float(0.0573782771)*x+float(0.3580258429))*x+float(0.9121133217))*x+float(0.9793091728);
  }
  if (x < float(0)){
    return ((float(0.1199175927)*x+float(0.4815668234))*x+float(0.9975991939))*x+float(0.9999505077);
  }
  return (x > float(46.052) ? float(1e20) : expf(x));
}

int Utils::QuickSelectPartition(vector<pair<float, int> >& scores, int lower, int upper){
  float pivot = scores[upper].first;
  while (lower < upper) {
    while (scores[lower].first < pivot) ++lower;
    while (scores[upper].first > pivot) --upper;
    if (scores[lower].first == scores[upper].first) ++lower;
    else if (lower < upper) swap(scores[lower], scores[upper]);
    
  }
  return upper;
}

float Utils::QuickSelect(vector<pair<float, int> >& scores, int lower,int upper, int k){
  if ( lower == upper ) return scores[lower].first;
  int split = QuickSelectPartition(scores, lower, upper);
  int length = split - lower + 1;
  if (length == k) return scores[split].first;
  else if (k  < length) return QuickSelect(scores, lower, split-1, k);
  else return QuickSelect(scores, split+1, upper, k - length);
}

void Utils::PruneBeam(unordered_map<int, State> &candidate_list, vector<State>& outer_sum){
  vector<pair<float, int> > scores;
  
  for (auto &it : candidate_list) {
    int i = it.first;
    State &cand = it.second;
    float new_alpha = (i >= 1 ? outer_sum[i-1].alpha : 0) + cand.alpha;
   
    scores.push_back(make_pair(new_alpha, i));
  }
  if (scores.size() <= _beam_size){
    return;
  }
  float threshold = QuickSelect(scores, 0, scores.size() - 1, scores.size() - _beam_size);
  
  for (auto &p : scores) {
    if (p.first < threshold){
      candidate_list.erase(p.second);
    }
  }
  
  return;
}
