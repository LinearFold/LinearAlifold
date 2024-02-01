#include "linear_partition.h"

void LinearPartition::Run(unordered_map<pair<int,int>, float, hash_pair>& bpp_matrix){
  Initialize();
  Inside();
  Outside();
  CalcBPP(bpp_matrix);

  return;
}

void LinearPartition::Inside(){  
  if(_seq_length > 0){
    _outer_sum[0].alpha = 0.0;
  }
  if(_seq_length > 1){
    _outer_sum[1].alpha = 0.0;
  }

  for(int i = 0; i < _seq_length; i++){
    int nuc_i = _int_seq[i];
    int nuc_i_p1 = (i+1) < _seq_length ? _int_seq[i+1] : -1;

    //Hairpin Search
    if(_hairpin_sum[i].size() > _beam_size){
      _utils.PruneBeam(_hairpin_sum[i], _outer_sum);
    }
    
    int i_pair = _next_pair[nuc_i][i];
    while(i_pair - i <= TURN && i_pair != -1){
      i_pair = _next_pair[nuc_i][i_pair];
    }
    
    if (i_pair != -1) {
      int nuc_i_pair = _int_seq[i_pair];
      int nuc_i_pair_m1 = i_pair > 0 ? _int_seq[i_pair-1] : -1;
      float new_score = _utils.HairpinEnergy(i, i_pair, nuc_i, nuc_i_p1, nuc_i_pair_m1, nuc_i_pair);
      _utils.FastLogPlusEquals(_hairpin_sum[i_pair][i].alpha, -new_score/kT);
    }

    for (auto &it : _hairpin_sum[i]) {
      int h = it.first;
      State &state = it.second;
      _utils.FastLogPlusEquals(_stem_sum[i][h].alpha, state.alpha);
      
      int nuc_h = _int_seq[h];
      int h_pair = _next_pair[nuc_h][i];
      if (h_pair != -1) {
	int nuc_h_p1 = (h+1) < _seq_length ? _int_seq[h+1] : -1;
	int nuc_h_pair = _int_seq[h_pair];
	int nuc_h_pair_m1 = h_pair > 0 ? _int_seq[h_pair-1] : -1;
	
	float new_score = _utils.HairpinEnergy(h, h_pair, nuc_h, nuc_h_p1, nuc_h_pair_m1, nuc_h_pair);	
	_utils.FastLogPlusEquals(_hairpin_sum[h_pair][h].alpha, -new_score/kT);
      }
     
    }    
    if (i == 0){
      continue;
    }
    if(_multi_sum[i].size() > _beam_size){
      _utils.PruneBeam(_multi_sum[i], _outer_sum);
    }
    Multi(i, 0);    
    
    if(_stem_sum[i].size() > _beam_size){
      _utils.PruneBeam(_stem_sum[i], _outer_sum);
    }
    Stem(i, 0);
    
    if(_multi2_sum[i].size() > _beam_size){
      _utils.PruneBeam(_multi2_sum[i], _outer_sum);
    }    
    Multi2(i, 0);
    
    if (_multi1_sum[i].size() > _beam_size){
      _utils.PruneBeam(_multi1_sum[i], _outer_sum);
    }
    Multi1(i, 0);

    if (i < _seq_length-1) {
      _utils.FastLogPlusEquals(_outer_sum[i+1].alpha, _outer_sum[i].alpha);
    }
  }
  
}

void LinearPartition::Outside(){
  _outer_sum[_seq_length-1].beta = 0.0;

  _beta_pf = MIN_VALUE;
  for(int i = _seq_length-1; i >= 0; i--) {
    //Outer
    if(i < _seq_length-1){
      _utils.FastLogPlusEquals(_outer_sum[i].beta,  _outer_sum[i+1].beta);
    }
    
    Multi1(i, 1);
    Multi2(i, 1);
    Stem(i, 1);
    Multi(i, 1);
  }
 
  _utils.FastLogPlusEquals(_beta_pf, _outer_sum[0].beta);
  //cout << -kT *_outer_sum[_seq_length-1].alpha/100.0 << endl;
  //cout << -kT *_beta_pf/100.0 << endl;
  return;
}

void LinearPartition::CalcBPP(unordered_map<pair<int,int>, float, hash_pair>& bpp_matrix){
  for(int i=0; i < _seq_length; i++){
    for(auto &it : _stem_sum[i]){
      int h = it.first;
      State state = it.second;
            
      float log_prob = state.alpha + state.beta - _beta_pf;
      if (log_prob > float(-9.91152)) {
	float prob = _utils.FastExp(log_prob);
	
	if(prob > float(1.0)){
	  prob = float(1.0);
	}
	if(prob < float(_bpp_threshold)){
	  continue;
	}
	bpp_matrix[make_pair(h, i)] = prob;
      }
    }
  }
  return;
}

void LinearPartition::Stem(int i, int flag){
  int nuc_i = _int_seq[i];
  for(auto & it : _stem_sum[i]){
    int h = it.first;
    State& state = it.second;
    int nuc_h = _int_seq[h];
    
    if(h > 0 && i <_seq_length-1){
      for (int p = h-1;p >= max(h-MAX_LOOP-1, 0); p--){
	int nuc_p = _int_seq[p];
	int nuc_p_p1 = _int_seq[p+1];
	int q = _next_pair[nuc_p][i];
	  
	while(q != -1 && ((h - p) + (q - i) - 2 <= MAX_LOOP)){
	  int nuc_q = _int_seq[q];
	  int nuc_q_m1 = _int_seq[q-1];

	  int type = BP_pair[nuc_p][nuc_q];
	  int type2 = BP_pair[nuc_i][nuc_h];

	  float new_score = _utils.LoopEnergy(type, type2, h, i, p, q, _int_seq);
	  if(flag == 0){
	    _utils.FastLogPlusEquals(_stem_sum[q][p].alpha, state.alpha - new_score/kT);
	  }else{
	    _utils.FastLogPlusEquals(state.beta, _stem_sum[q][p].beta - new_score/kT);
	  }
	  q = _next_pair[nuc_p][q];
	}
      }
    }
      
    if(h > 0 && i < _seq_length-1){
      float new_score = ML_intern37 + _utils.DangleEnergy(BP_pair[nuc_h][nuc_i], h-1, i+1, _seq_length, _int_seq);
      if(flag == 0){
	_utils.FastLogPlusEquals(_multi1_sum[i][h].alpha, state.alpha - new_score/kT);
      }else{
	_utils.FastLogPlusEquals(state.beta, _multi1_sum[i][h].beta - new_score/kT);
      }
    }
      
    if (h > 1 && !_multi1_sum[h-1].empty()) {

      float multi_energy = ML_intern37 + _utils.DangleEnergy(BP_pair[nuc_h][nuc_i], h-1, i+1, _seq_length, _int_seq);
      for (auto &m : _multi1_sum[h-1]) {
	int new_h = m.first;
	State& m_state = m.second;
	if(flag == 0){
	  _utils.FastLogPlusEquals(_multi2_sum[i][new_h].alpha, m_state.alpha + state.alpha - multi_energy/kT);
	}else{
	  _utils.FastLogPlusEquals(state.beta, _multi2_sum[i][new_h].beta + m_state.alpha - multi_energy/kT);
	  _utils.FastLogPlusEquals(m_state.beta, _multi2_sum[i][new_h].beta + state.alpha - multi_energy/kT);
	}
      }
    }
      
    if(h >= 1){	
      float new_score = _utils.DangleEnergy(BP_pair[nuc_h][nuc_i], h-1, i+1, _seq_length, _int_seq);
      if(flag == 0){
	_utils.FastLogPlusEquals(_outer_sum[i].alpha, _outer_sum[h-1].alpha + state.alpha - new_score/kT);
      }else{
	_utils.FastLogPlusEquals(_outer_sum[h-1].beta, state.alpha + _outer_sum[i].beta - new_score/kT);
	_utils.FastLogPlusEquals(state.beta, _outer_sum[h-1].alpha + _outer_sum[i].beta - new_score/kT);
      }
    } else {	
      float new_score = _utils.DangleEnergy(BP_pair[_int_seq[0]][nuc_i], -1, i+1, _seq_length, _int_seq);
      if(flag == 0){
	_utils.FastLogPlusEquals(_outer_sum[i].alpha, state.alpha - new_score/kT);
      }else{
	_utils.FastLogPlusEquals(_beta_pf, state.alpha+ _outer_sum[i].beta - new_score/kT);
	_utils.FastLogPlusEquals(state.beta, _outer_sum[i].beta - new_score/kT);
      }
    }
      
  }
}

void LinearPartition::Multi1(int i, int flag){
  for(auto& it : _multi1_sum[i]) {
    int h = it.first;
    State& state = it.second;
    if (i < _seq_length-1) {
      float new_score = ML_BASE37;
      if(flag == 0){
	_utils.FastLogPlusEquals(_multi1_sum[i+1][h].alpha, state.alpha - new_score/kT);
      }else{
	_utils.FastLogPlusEquals(state.beta, _multi1_sum[i+1][h].beta - new_score/kT);
      }
    }
  }

}

void LinearPartition::Multi2(int i, int flag){
  for(auto& it : _multi2_sum[i]) {
    int h = it.first;
    State& state = it.second;
    if(flag == 0){
      _utils.FastLogPlusEquals(_multi1_sum[i][h].alpha, state.alpha);
    }else{
      _utils.FastLogPlusEquals(state.beta, _multi1_sum[i][h].beta);
    }
    
    for(int g = h-1; g >= max(0, h - MAX_LOOP -1); g--){
      int nuc_g = _int_seq[g];
      int g_pair = _next_pair[nuc_g][i];
      if(g_pair != -1){
	float new_score = ML_BASE37*(h-g-1) + ML_BASE37*(g_pair-i-1);
	if(flag == 0){
	  _utils.FastLogPlusEquals(_multi_sum[g_pair][g].alpha, state.alpha - new_score/kT);
	}else{
	  _utils.FastLogPlusEquals(state.beta, _multi_sum[g_pair][g].beta - new_score/kT);
	}
      }
    }
  }
  
}

void LinearPartition::Multi(int i, int flag){
  int nuc_i = _int_seq[i];
  
  for(auto& it : _multi_sum[i]){
    int h = it.first;
    State& state = it.second;
    
    int nuc_h = _int_seq[h];
    float new_score = ML_intern37 + ML_closing37 + _utils.DangleEnergy(BP_pair[nuc_i][nuc_h], i-1, h+1, _seq_length, _int_seq);
    if(flag == 0){
      _utils.FastLogPlusEquals(_stem_sum[i][h].alpha, state.alpha - new_score/kT);
    }else{
      _utils.FastLogPlusEquals(state.beta, _stem_sum[i][h].beta - new_score/kT);
    }
    
    int h_pair = _next_pair[nuc_h][i];
    
    if (h_pair != -1 && h_pair - i <= MAX_LOOP+1) {
      new_score = ML_BASE37*(h_pair-i);
      if(flag == 0){
	_utils.FastLogPlusEquals(_multi_sum[h_pair][h].alpha, state.alpha - new_score/kT);
      }else{
	_utils.FastLogPlusEquals(state.beta, _multi_sum[h_pair][h].beta - new_score/kT);
      }
    }
  }
}

void LinearPartition::Initialize(){
  _seq_length = _int_seq.size();

  _hairpin_sum.clear();
  _hairpin_sum.resize(_seq_length);
  _stem_sum.clear();
  _stem_sum.resize(_seq_length);
  _multi2_sum.clear();
  _multi2_sum.resize(_seq_length);
  _multi1_sum.clear();
  _multi1_sum.resize(_seq_length);
  _multi_sum.clear();
  _multi_sum.resize(_seq_length);
  _outer_sum.clear();
  _outer_sum.resize(_seq_length);

  _next_pair.resize(NUM_OF_NUC_TYPES, vector<int>(_seq_length, -1));
  for(int nuc_i = 0; nuc_i < NUM_OF_NUC_TYPES; nuc_i++){
    int next = -1;
    for(int j = _seq_length-1; j >=0; j--){
      _next_pair[nuc_i][j] = next;
      if(BP_pair[nuc_i][_int_seq[j]]){
	next = j;
      }
    }
  }
}
