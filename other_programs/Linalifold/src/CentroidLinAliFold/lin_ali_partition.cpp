#include "lin_ali_partition.h"

void LinAliPartition::Run(unordered_map<pair<int,int>, float, hash_pair>& bpp_matrix){
  Initialize();
  Inside();
  Outside();
  CalcBPP(bpp_matrix);
  Clear();
  return;
}

void LinAliPartition::Clear(){
  for(int i = 0; i < _allowed_base_pair.size(); i++){
    _allowed_base_pair[i].clear();
    _allowed_base_pair[i].shrink_to_fit();
  }
  _allowed_base_pair.clear();
  _allowed_base_pair.shrink_to_fit();

  _all_ribosum_data.clear();
  
  for(int i = 0; i < _alignment_length; i++){
    _stem_sum[i].clear();
    _multi2_sum[i].clear();
    _multi1_sum[i].clear();
    _multi_sum[i].clear();
  }
  _stem_sum.clear();
  _stem_sum.shrink_to_fit();
  _multi2_sum.clear();
  _multi2_sum.shrink_to_fit();
  _multi1_sum.clear();
  _multi1_sum.shrink_to_fit();
  _multi_sum.clear();
  _multi_sum.shrink_to_fit();
  _outer_sum.clear();
  _outer_sum.shrink_to_fit();
}

void LinAliPartition::Inside(){
  if(_alignment_length > 0){
    _outer_sum[0].alpha = 0.0;
  }
  if(_alignment_length > 1){
    _outer_sum[1].alpha = 0.0;
  }

  for(int i = 0; i < _alignment_length; i++){
    // Hairpin Search
    for(int h = max(0, i-_beam_size); h < i-TURN; h++){
      double cons_score = AllowBasePairs(h, i);
      
      if(cons_score > _cons_score_threshold){
	float new_score = 0;
	
	for(int s = 0; s < _num_of_seq; s++){
	  int index_h = _a2s_index[s][h];
	  int index_i = _a2s_index[s][i];
	  
	  if(index_h < 0 || index_i < 0 || index_i - index_h <= TURN){
	    new_score += ILLEGAL_HAIRPIN_SCORE;	   
	  }else{
	    int nuc_h = _int_seq_list[s][index_h];
	    int nuc_i = _int_seq_list[s][index_i];	    	    
	    int nuc_h_p1 = _int_seq_list[s][index_h+1];
	    int nuc_i_m1 = _int_seq_list[s][index_i-1];
	    
	    new_score += _utils.HairpinEnergy(h, i, nuc_h, nuc_h_p1, nuc_i_m1, nuc_i);
	  }	  
	}
	new_score = new_score/_num_of_seq - cons_score;	  
	_utils.FastLogPlusEquals(_stem_sum[i][h].alpha, -new_score/kT);
      }
    }
    if(i == 0){
      continue;
    }    
 
    if(_multi_sum[i].size() > _beam_size){
      _utils.PruneBeam(_multi_sum[i], _outer_sum);
    }
    Multi(i,0);
    
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
    
    //Outer
    if(i < _alignment_length-1){
      _utils.FastLogPlusEquals(_outer_sum[i+1].alpha, _outer_sum[i].alpha);
    }
  }
}

void LinAliPartition::Outside(){
  _outer_sum[_alignment_length-1].beta = 0.0;

  _beta_pf = MIN_VALUE;
  for(int i = _alignment_length-1; i >= 0; i--) {
    //Outer
    if(i < _alignment_length-1){
      float new_score = _outer_sum[i+1].beta; // add outer unpaired base (energy:0);
      _utils.FastLogPlusEquals(_outer_sum[i].beta, new_score);
    }
    
    Multi1(i, 1);
    Multi2(i, 1);
    Stem(i, 1);
    Multi(i, 1);
  }
 
  _utils.FastLogPlusEquals(_beta_pf, _outer_sum[0].beta);
  //cout << -kT *_outer_sum[_alignment_length-1].alpha/100.0 << endl;
  //cout << -kT *_beta_pf/100.0 << endl;
  return;
}

void LinAliPartition::CalcBPP(unordered_map<pair<int,int>, float, hash_pair>& bpp_matrix){
  for(int i=0; i< _alignment_length; i++){
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

void LinAliPartition::Stem(int i, int flag){  
  for(auto &it : _stem_sum[i]){      
    int h = it.first;
    State& state = it.second;
    
    if(h > 0 && i < _alignment_length-1){
      float new_score = ML_intern37*_num_of_seq;
      for(int s = 0; s < _num_of_seq; s++){
	int index_h = _a2s_index[s][h];
	int index_i = _a2s_index[s][i];
	
	if(index_h >= 0 && index_i >= 0){
	  int nuc_h = _int_seq_list[s][index_h];
	  int nuc_i = _int_seq_list[s][index_i];
	  new_score += _utils.DangleEnergy(BP_pair[nuc_h][nuc_i], index_h-1, index_i+1, _int_seq_list[s].size(), _int_seq_list[s]);
	}
      }
      new_score = new_score/_num_of_seq;
      if(flag == 0){
	_utils.FastLogPlusEquals(_multi1_sum[i][h].alpha, state.alpha - new_score/kT);
      }else{
	_utils.FastLogPlusEquals(state.beta, _multi1_sum[i][h].beta - new_score/kT);
      }
    }
    
    
    if (h >= 1 && !_multi1_sum[h-1].empty()) {
      float multi1_score =  ML_intern37*_num_of_seq;
      for(int s = 0; s < _num_of_seq; s++){
	int index_h = _a2s_index[s][h];
	int index_i = _a2s_index[s][i];
	
	if(index_h >= 0 && index_i >= 0){
	  int nuc_h = _int_seq_list[s][index_h];
	  int nuc_i = _int_seq_list[s][index_i];
	  multi1_score += _utils.DangleEnergy(BP_pair[nuc_h][nuc_i], index_h-1, index_i+1, _int_seq_list[s].size(), _int_seq_list[s]);
	}
      }
      multi1_score = multi1_score/_num_of_seq;
      for (auto &m : _multi1_sum[h-1]) {
	int new_h = m.first;
	State& m_state = m.second;
	
	if(flag == 0){
	  _utils.FastLogPlusEquals(_multi2_sum[i][new_h].alpha, m_state.alpha + state.alpha - multi1_score/kT);
	}else{
	  _utils.FastLogPlusEquals(state.beta, _multi2_sum[i][new_h].beta + m_state.alpha - multi1_score/kT);
	  _utils.FastLogPlusEquals(m_state.beta, _multi2_sum[i][new_h].beta + state.alpha - multi1_score/kT);
	}
      }
    }

    {
      float new_score = 0;      
      for(int s = 0; s < _num_of_seq; s++){
	int index_h = _a2s_index[s][h];
	int index_i = _a2s_index[s][i];
	if(index_h >= 0 && index_i >= 0){
	  int nuc_h = _int_seq_list[s][index_h];
	  int nuc_i = _int_seq_list[s][index_i];
	  
	  new_score += _utils.DangleEnergy(BP_pair[nuc_h][nuc_i], index_h-1, index_i+1, _int_seq_list[s].size(), _int_seq_list[s]);
	}
      }
      new_score = new_score/_num_of_seq;
      
      if(h >= 1){
	if(flag == 0){
	  _utils.FastLogPlusEquals(_outer_sum[i].alpha, _outer_sum[h-1].alpha + state.alpha - new_score/kT);
	}else{
	  _utils.FastLogPlusEquals(_outer_sum[h-1].beta, state.alpha + _outer_sum[i].beta - new_score/kT);
	  _utils.FastLogPlusEquals(state.beta, _outer_sum[h-1].alpha + _outer_sum[i].beta - new_score/kT);
	}
      } else {
	if(flag == 0){
	  _utils.FastLogPlusEquals(_outer_sum[i].alpha, state.alpha - new_score/kT);
	}else{	  	  
	  _utils.FastLogPlusEquals(_beta_pf, state.alpha + _outer_sum[i].beta - new_score/kT);
	  _utils.FastLogPlusEquals(state.beta, _outer_sum[i].beta - new_score/kT);
	}
      }
    }
   
    if(h > 0 && i < _alignment_length-1){
      for (int p = h-1; p >= max(h-MAX_LOOP-1, 0); p--){
	for(int q = i+1; (h - p) + (q - i) - 2 <= MAX_LOOP && q < _alignment_length ;q++){
	  double cons_score = AllowBasePairs(p, q);
	  if(cons_score > _cons_score_threshold){
	    float new_score = 0;
	    for(int s = 0; s < _num_of_seq; s++){
	      int index_h = _a2s_index[s][h];
	      int index_i = _a2s_index[s][i];
	      int index_p = _a2s_index[s][p];
	      int index_q = _a2s_index[s][q];
	      
	      if(index_h >= 0 && index_i >= 0 && index_p >= 0 && index_q >= 0){
		int nuc_h = _int_seq_list[s][index_h];
		int nuc_i = _int_seq_list[s][index_i];
		int nuc_p = _int_seq_list[s][index_p];
		int nuc_q = _int_seq_list[s][index_q];
		
		int type = BP_pair[nuc_p][nuc_q];
		int type2 = BP_pair[nuc_i][nuc_h];
		
		new_score += _utils.LoopEnergy(type, type2, index_h, index_i, index_p, index_q, _int_seq_list[s]);
	      }
	    }
	    new_score = new_score/_num_of_seq - cons_score;
	    if(flag == 0){
	      _utils.FastLogPlusEquals(_stem_sum[q][p].alpha, state.alpha - new_score/kT);
	    }else{
	      _utils.FastLogPlusEquals(state.beta, _stem_sum[q][p].beta - new_score/kT);
	    }
	  }
	}	 
      }
    }
  }
}

void LinAliPartition::Multi1(int i, int flag){
  for (auto &it : _multi1_sum[i]){
    int h = it.first;
    State& state = it.second;
    if(i+1 != _alignment_length){
      float new_score = 0;
      for(int s = 0; s < _num_of_seq; s++){
	int index_h = _a2s_index[s][h];
	if(index_h >= 0){
	  new_score += ML_BASE37;
	}  
      }
      new_score = new_score/_num_of_seq;
      if(flag == 0){
	_utils.FastLogPlusEquals(_multi1_sum[i+1][h].alpha, state.alpha - new_score/kT); 
      }else{
	_utils.FastLogPlusEquals(state.beta, _multi1_sum[i+1][h].beta - new_score/kT);
      }
    }
  }
}

void LinAliPartition::Multi2(int i, int flag){ 
  for(auto &it : _multi2_sum[i]){      
    int h = it.first;
    State& state = it.second;
    if(flag == 0){
      _utils.FastLogPlusEquals(_multi1_sum[i][h].alpha, state.alpha);
    }else{
      _utils.FastLogPlusEquals(state.beta, _multi1_sum[i][h].beta);
    }
    
    for(int g = h-1; g >= max(0, h-MAX_LOOP-1); g--){
      for(int g_pair = i+1; g_pair < _alignment_length; g_pair++){
	double cons_score = AllowBasePairs(g, g_pair);
	
	if(cons_score > _cons_score_threshold){
	  //approximation
	  float new_score = ML_BASE37*((h-g-1)+(g_pair-i-1));
	  if(flag == 0){
	    _utils.FastLogPlusEquals(_multi_sum[g_pair][g].alpha, state.alpha - new_score/kT);
	  }else{
	    _utils.FastLogPlusEquals(state.beta, _multi_sum[g_pair][g].beta - new_score/kT);
	  }
	  break;
	}	  
      }
    }
  }
}

void LinAliPartition::Multi(int i, int flag){
  unordered_map<int, State>& stem_candidate = _stem_sum[i];
  
  for(auto &it : _multi_sum[i]){
    int h = it.first;
    State& state = it.second;
    
    double cons_score = AllowBasePairs(h, i);
    float new_score = (ML_intern37 + ML_closing37)*_num_of_seq;
    
    for(int s = 0; s < _num_of_seq; s++){
      int index_h = _a2s_index[s][h];
      int index_i = _a2s_index[s][i];
      
      if(index_h >= 0 && index_i >= 0){
	int nuc_h = _int_seq_list[s][index_h];
	int nuc_i = _int_seq_list[s][index_i];
	new_score += _utils.DangleEnergy(BP_pair[nuc_i][nuc_h], index_i-1, index_h+1, _int_seq_list[s].size(), _int_seq_list[s]);
      }
    }
    new_score = new_score/_num_of_seq - cons_score;
    if(flag == 0){
      _utils.FastLogPlusEquals(_stem_sum[i][h].alpha, state.alpha - new_score/kT);
    }else{
      _utils.FastLogPlusEquals(state.beta, _stem_sum[i][h].beta - new_score/kT);
    }

    
    for(int h_pair = i+1; h_pair <= min(_alignment_length-1, i+MAX_LOOP+1); h_pair++){
      double cons_score = AllowBasePairs(h, h_pair);
      if(cons_score > _cons_score_threshold){
	//approximation
	float new_score = ML_BASE37*(h_pair-i);
	if(flag == 0){
	  _utils.FastLogPlusEquals(_multi_sum[h_pair][h].alpha, state.alpha - new_score/kT);
	}else{
	  _utils.FastLogPlusEquals(state.beta, _multi_sum[h_pair][h].beta - new_score/kT);
	}
	break;
      }
    }
  }    
}

double LinAliPartition::AllowBasePairs(int h, int i){
  if(_allowed_base_pair[h][i-h] != 0.0){
    return(_allowed_base_pair[h][i-h]);
  }
  int count = 0;
  vector<int> type_freq(8,0);
  
  for(int s = 0; s < _num_of_seq; s++){
    int index_h = _a2s_index[s][h];
    int index_i = _a2s_index[s][i];    
    
    if(index_h >= 0 && index_i >= 0){
      int nuc_h = _int_seq_list[s][index_h];
      int nuc_i = _int_seq_list[s][index_i];
      int type = BP_pair[nuc_h][nuc_i];
      type_freq[type]++;
      if(!type || index_i - index_h <= TURN){
	count++;
      }
    }else{
      if(index_h < 0 && index_i < 0){
	type_freq[7]++;
      }else{
	type_freq[0]++;
      }
      count++;
    }
  }
  

  double score = 0.0;
  for (int j = 1; j<=6; j++){
    for(int k = j; k<=6; k++){
      score += type_freq[j]*type_freq[k]*_ribosum_matrix[j][k];
    }
  }    

  if (2*type_freq[0]+type_freq[7] < _num_of_seq) { 
    score = _beta * ((SCALE*score)/_num_of_seq - _delta*SCALE*(type_freq[0] + type_freq[7]*0.25));
    score = score/_num_of_seq;
    // //average      
  }else{
    score = NEG_INF;
  }
  _allowed_base_pair[h][i-h] = score;
  return(score);
}

void LinAliPartition::Initialize(){
  _stem_sum.clear();
  _stem_sum.resize(_alignment_length);
  _multi2_sum.clear();
  _multi2_sum.resize(_alignment_length);
  _multi1_sum.clear();  
  _multi1_sum.resize(_alignment_length);
  _multi_sum.clear();
  _multi_sum.resize(_alignment_length);
  _outer_sum.clear();
  _outer_sum.resize(_alignment_length);

  _allowed_base_pair.resize(_alignment_length);
  for(int i = 0; i < _alignment_length; i++){
    _allowed_base_pair[i].resize(_alignment_length-i, 0.0);
  }

  if(_ribosum_flag == 1){
    GetRibosumMatrix();
  }else{    
    _ribosum_matrix =
      {{0,0,0,0,0,0,0}, 
       {0,0,2,2,1,2,2},
       {0,2,0,1,2,2,2},
       {0,2,1,0,2,1,2},
       {0,1,2,2,0,2,1},
       {0,2,2,1,2,0,2},
       {0,2,2,2,1,2,0}};
  }
}
