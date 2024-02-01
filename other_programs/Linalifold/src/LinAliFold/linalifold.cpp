#include "linalifold.h"

void LinAliFold::Run(){
  // for runtime statistics
  struct timeval parse_starttime, parse_endtime;
  gettimeofday(&parse_starttime, NULL);

  ReadData();
  
  Initialize(); 
  Parse();
  TraceBack();
 
  cout << fixed;
  cout << "score:" << setprecision(2) << -_outer_best[_alignment_length-1].score/100.0 << endl;
  cout << _structure << endl;

  gettimeofday(&parse_endtime, NULL);
  double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec - parse_starttime.tv_usec) / 1000000.0;
  fprintf(stdout, "LinAliFold Calculation Time: %.3f seconds.\n", parse_elapsed_time);
}

void LinAliFold::Parse(){
  if(_alignment_length > 0){
    _outer_best[0].set(0, OUTER_eq_OUTER_plus_UNPAIRED);
  }
  if(_alignment_length > 1){
    _outer_best[1].set(0, OUTER_eq_OUTER_plus_UNPAIRED);
  }

  for(int i = 0; i < _alignment_length; i++){
    unordered_map<int, State>& stem_candidate = _stem_best[i];
    unordered_map<int, State>& multi_candidate = _multi_best[i];
    unordered_map<int, State>& multi2_candidate = _multi2_best[i];
    unordered_map<int, State>& multi1_candidate = _multi1_best[i];
    State& outer_candidate = _outer_best[i];
    
    // Hairpin Search
    for(int h = max(0, i-_beam_size); h < i-TURN; h++){
      double cons_score = AllowBasePairs(h, i);

      if(cons_score > _conservation_score_threshold){
	double new_score = 0;
	
	for(int s = 0; s < _num_of_seq; s++){
	  int index_h = _a2s_index[s][h];
	  int index_i = _a2s_index[s][i];
	  
	  if(index_h < 0 || index_i < 0 || index_i - index_h <= TURN){
	    new_score -= ILLEGAL_HAIRPIN_SCORE;	   
	  }else{
	    int nuc_h = _int_seq_list[s][index_h];
	    int nuc_i = _int_seq_list[s][index_i];	    	    
	    int nuc_h_p1 = _int_seq_list[s][index_h+1];
	    int nuc_i_m1 = _int_seq_list[s][index_i-1];
	    new_score -= HairpinEnergy(index_h, index_i, nuc_h, nuc_h_p1, nuc_i_m1, nuc_i);
	  }
	}
	UpdateIfBetter(stem_candidate[h], cons_score+new_score/_num_of_seq, HAIRPIN);	
      }
    }
    // Multi Beam Search

    if(multi_candidate.size() > _beam_size){
      PruneBeam(multi_candidate);
    }
    for(auto &it : multi_candidate){
      int h = it.first;
      State& state = it.second;

      double cons_score = AllowBasePairs(h, i);
      double new_score = -(ML_intern37 + ML_closing37)*_num_of_seq;
      for(int s = 0; s < _num_of_seq; s++){
	int index_h = _a2s_index[s][h];
	int index_i = _a2s_index[s][i];
	
	if(index_h >= 0 && index_i >= 0){
	  int nuc_h = _int_seq_list[s][index_h];
	  int nuc_i = _int_seq_list[s][index_i];
	  new_score -= DangleEnergy(BP_pair[nuc_i][nuc_h], index_i-1, index_h+1, s);
	}
      }
      
      UpdateIfBetter(stem_candidate[h], state.score + cons_score + new_score/_num_of_seq, S_eq_MULTI);

      for(int h_pair = i+1; h_pair <= min(_alignment_length-1, i+MAXLOOP+1); h_pair++){
	double cons_score = AllowBasePairs(h, h_pair);
	if(cons_score > _conservation_score_threshold){
	  int new_l1 = state.l1;
	  int new_l2 = state.l2+h_pair - i;
	  //approximation
	  int new_score = -ML_BASE37*_num_of_seq*(h_pair-i);	 
	  UpdateIfBetter(_multi_best[h_pair][h], state.score + new_score/_num_of_seq, MULTI_eq_MULTI_plus_UNPAIRED, new_l1, new_l2);
	  break;
	}
      }
    }
    // Stem Beam Search
    
    if(stem_candidate.size() > _beam_size){
      PruneBeam(stem_candidate);
    }

    bool use_cube_pruning = (_beam_size > CUBE_PRUNING_SIZE) && (stem_candidate.size() > CUBE_PRUNING_SIZE);

    for(auto &it : stem_candidate){      
      int h = it.first;      
      State& state = it.second;
      if(h > 0 && i < _alignment_length-1){
	double new_score = -ML_intern37*_num_of_seq;
	for(int s = 0; s < _num_of_seq; s++){
	  int index_h = _a2s_index[s][h];
	  int index_i = _a2s_index[s][i];
	  
	  if(index_h >= 0 && index_i >= 0){
	    int nuc_h = _int_seq_list[s][index_h];
	    int nuc_i = _int_seq_list[s][index_i];
	    new_score -= DangleEnergy(BP_pair[nuc_h][nuc_i], index_h-1, index_i+1, s);
	  }
	}
	UpdateIfBetter(multi1_candidate[h], state.score + new_score/_num_of_seq, MULTI1_eq_STEM);
	
      }
      
      if(!use_cube_pruning){
	if (h >= 1 && !_multi1_best[h-1].empty()) {
	  double multi1_score =  -ML_intern37*_num_of_seq;
	  for(int s = 0; s < _num_of_seq; s++){
	    int index_h = _a2s_index[s][h];
	    int index_i = _a2s_index[s][i];

	    if(index_h >= 0 && index_i >= 0){
	      int nuc_h = _int_seq_list[s][index_h];
	      int nuc_i = _int_seq_list[s][index_i];
	      multi1_score -= DangleEnergy(BP_pair[nuc_h][nuc_i], index_h-1, index_i+1, s);
	    }
	  }
	  multi1_score = state.score + multi1_score/_num_of_seq;
	  auto multi2_iter = multi2_candidate.find(h);
	  if (multi2_iter==multi2_candidate.end() || multi1_score > multi2_iter->second.score) {
	    for (auto &m : _multi1_best[h-1]) {
	      int new_h = m.first;
	      double new_score = multi1_score + m.second.score;
	      
	      UpdateIfBetter(multi2_candidate[new_h], new_score, MULTI2_eq_MULTI1_plus_STEM, h-1);
	    
	    }
	  }
	}
      }
	
      if(h >= 1){
	State& outer_prefix = _outer_best[h-1];	
	if (outer_prefix.type != NONE) {
	  double new_score = 0;
	  for(int s = 0; s < _num_of_seq; s++){
	    int index_h = _a2s_index[s][h];
	    int index_i = _a2s_index[s][i];
	    if(index_h >= 0 && index_i >= 0){
	      int nuc_h = _int_seq_list[s][index_h];
	      int nuc_i = _int_seq_list[s][index_i];
	      
	      new_score -= DangleEnergy(BP_pair[nuc_h][nuc_i], index_h-1, index_i+1, s);
	    }
	  }
	  
	  UpdateIfBetter(outer_candidate,  outer_prefix.score + state.score + new_score/_num_of_seq, OUTER_eq_OUTER_plus_STEM, h-1);
	}
      } else {
	double new_score = 0;
	for(int s = 0; s < _num_of_seq; s++){
	  int index_h = _a2s_index[s][0];
	  int index_i = _a2s_index[s][i];
	  if(index_i >= 0 && index_h >= 0){
	    int nuc_h = _int_seq_list[s][index_h];
	    int nuc_i = _int_seq_list[s][index_i];
	    new_score -= DangleEnergy(BP_pair[nuc_h][nuc_i], -1, index_i+1, s);	    
	  }	  
	}
	
	UpdateIfBetter(outer_candidate,  state.score +new_score/_num_of_seq, OUTER_eq_OUTER_plus_STEM, -1);	
      }

      if(h > 0 && i < _alignment_length-1){
	for (int p = h-1; p >= max(h-MAXLOOP-1, 0); p--){
	  for(int q = i+1; (h - p) + (q - i) - 2 <= MAXLOOP && q < _alignment_length ;q++){
	    double cons_score = AllowBasePairs(p, q);
	    if(cons_score > _conservation_score_threshold){
	      double new_score = 0;
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

		  new_score -= LoopEnergy(type, type2, index_h, index_i, index_p, index_q, s);
		}
	      }
	      
	   
	      if(p == h-1 && q == i+1){
		UpdateIfBetter(_stem_best[q][p], cons_score + state.score + new_score/_num_of_seq, STEM);		  
	      }else{
		UpdateIfBetter(_stem_best[q][p], cons_score + state.score + new_score/_num_of_seq, INTERNAL, h-p, q-i);
		
	      }

	      
	      
	    }
	  }	 
	}
      }      
    }
 
    if(use_cube_pruning){
      vector<int> valid_Stems;
      vector<int> multi1_scores;
      
      for (auto &it: stem_candidate) {
	int h = it.first;
	State &state = it.second;

	if (h >= 1 && !_multi1_best[h-1].empty()) {
	  double multi1_score =  -ML_intern37*_num_of_seq;
	  for(int s = 0; s < _num_of_seq; s++){
	    int index_h = _a2s_index[s][h];
	    int index_i = _a2s_index[s][i];

	    if(index_h >= 0 && index_i >= 0){
	      int nuc_h = _int_seq_list[s][index_h];
	      int nuc_i = _int_seq_list[s][index_i];
	      multi1_score -= DangleEnergy(BP_pair[nuc_h][nuc_i], index_h-1, index_i+1, s);
	    }	    
	  }
	  multi1_score = state.score + multi1_score/_num_of_seq;

	  auto multi2_iter = multi2_candidate.find(h);
	  if (multi2_iter==multi2_candidate.end() || multi1_score > multi2_iter->second.score) {
	    valid_Stems.push_back(h);
	    multi1_scores.push_back(multi1_score);
	  }
	}
      }

      
      vector<pair<int, pair<int, int> > > max_heap;
      //(heuristic score, (index of valid_stems(stem:h and i), rank of sorted_multi1[h-1])
      for (int s = 0; s < valid_Stems.size(); s++) {
	int h = valid_Stems[s];
	max_heap.push_back(make_pair(multi1_scores[s] + _sorted_multi1_best[h-1][0].first, make_pair(s, 0)));
	push_heap(max_heap.begin(), max_heap.end());
      }
	
      int filled = 0;
      int prev_score = MIN_VALUE;
      int current_score = MIN_VALUE;
      while ((filled < _beam_size || current_score == prev_score) && !max_heap.empty()) {
	auto &top = max_heap.front();
	prev_score = current_score;
	current_score = top.first;
	int index_S = top.second.first;
	int index_M = top.second.second;
	
	int h = valid_Stems[index_S];
	int new_h = _sorted_multi1_best[h-1][index_M].second;
	double new_score = multi1_scores[index_S] + _multi1_best[h-1][new_h].score;
	pop_heap(max_heap.begin(), max_heap.end());
	max_heap.pop_back();
	
	if (multi2_candidate[new_h].type == NONE) {
	  filled++;
	  UpdateIfBetter(multi2_candidate[new_h], new_score, MULTI2_eq_MULTI1_plus_STEM, h-1);
	}
	
	index_M++;
	while (index_M < _sorted_multi1_best[h-1].size()) {
	  int candidate_score = multi1_scores[index_S] + _sorted_multi1_best[h-1][index_M].first;
	  int candidate_new_h = _sorted_multi1_best[h-1][index_M].second;
	  if (multi2_candidate.find(candidate_new_h) == multi2_candidate.end()){
	    max_heap.push_back(make_pair(candidate_score, make_pair(index_S, index_M)));
	    push_heap(max_heap.begin(), max_heap.end());
	    break;
	  } else {
	    index_M++;
	  }
	}
      }
    }
 
    // Multi2 Beam Search    
    if(multi2_candidate.size() > _beam_size){
      PruneBeam(multi2_candidate);
    }

    for(auto &it : multi2_candidate){      
      int h = it.first;
      State& state = it.second;	  
      UpdateIfBetter(multi1_candidate[h], state.score, MULTI1_eq_MULTI2);

    

      for(int g = h-1; g >= max(0, h-MAXLOOP-1); g--){
	for(int g_pair = i+1; g_pair < _alignment_length; g_pair++){
	  double cons_score = AllowBasePairs(g, g_pair);
	  
	  if(cons_score > _conservation_score_threshold){
	    //approximation
	    double new_score = -ML_BASE37*_num_of_seq*((h-g-1)+(g_pair-i-1));
	    UpdateIfBetter(_multi_best[g_pair][g], state.score + new_score/_num_of_seq, MULTI, h-g, g_pair-i);
	    break;
	  }
	}
      }
    }
    // Multi1BeamSearch    
    int beam_prune_threshold = MIN_VALUE;
    if (multi1_candidate.size() > _beam_size){
      beam_prune_threshold = PruneBeam(multi1_candidate);
    }
    SortMulti1(beam_prune_threshold, multi1_candidate, _sorted_multi1_best[i]);

    for (auto &it : multi1_candidate){
      int h = it.first;
      State& state = it.second;
      if(i+1 != _alignment_length){
	double new_score = 0;
	for(int s = 0; s < _num_of_seq; s++){
	  int index_h = _a2s_index[s][h];
	  if(index_h >= 0){
	    new_score -= ML_BASE37;
	  }  
	}
	UpdateIfBetter(_multi1_best[i+1][h], state.score + new_score/_num_of_seq, MULTI1_eq_MULTI1_plus_UNPAIRED);
      }
    }
  
    //Outer
    if(i < _alignment_length-1){
      double new_score = outer_candidate.score; // add outer unpaired base (energy:0);
      UpdateIfBetter(_outer_best[i+1], new_score, OUTER_eq_OUTER_plus_UNPAIRED);
    }
  }
}

void LinAliFold::TraceBack(){
  stack<tuple<int, int, State> > stk;
  stk.push(make_tuple(0, _alignment_length-1, _outer_best[_alignment_length-1]));
  
  vector<pair<int,int> > multi_todo;
  unordered_map<int,int> mbp; // multi bp
  int k,p,q;

  double sum = 0.0;
  while(!stk.empty()){
    tuple<int, int, State> top = stk.top();
    int i = get<0>(top), j = get<1>(top);
    
    State& state = get<2>(top);
    stk.pop();
    switch(state.type){
    case HAIRPIN:
      _structure[i] = '(';
      _structure[j] = ')';
      sum += _allowed_base_pair[i][j-i];      
      break;
    case INTERNAL:
      _structure[i] = '(';
      _structure[j] = ')';
       sum += _allowed_base_pair[i][j-i];   
      p = i + state.l1;
      q = j - state.l2;
      stk.push(make_tuple(p, q, _stem_best[q][p]));
      break;
    case STEM:
      _structure[i] = '(';
      _structure[j] = ')';
       sum += _allowed_base_pair[i][j-i];   
      stk.push(make_tuple(i+1, j-1, _stem_best[j-1][i+1]));
      break;
    case MULTI:      
      p = i + state.l1;
      q = j - state.l2;
      stk.push(make_tuple(p, q, _multi2_best[q][p]));
      break;
    case MULTI_eq_MULTI_plus_UNPAIRED:
      p = i + state.l1;
      q = j - state.l2;
      stk.push(make_tuple(p, q, _multi2_best[q][p]));
      break;
    case S_eq_MULTI:
      _structure[i] = '(';
      _structure[j] = ')';
      stk.push(make_tuple(i, j, _multi_best[j][i]));
      break;
    case MULTI2_eq_MULTI1_plus_STEM:
      k = state.split;
      stk.push(make_tuple(i, k, _multi1_best[k][i]));
      stk.push(make_tuple(k+1, j, _stem_best[j][k+1]));
      break;
    case MULTI1_eq_MULTI2:
      stk.push(make_tuple(i, j, _multi2_best[j][i]));
      break;
    case MULTI1_eq_MULTI1_plus_UNPAIRED:
      stk.push(make_tuple(i, j-1, _multi1_best[j-1][i]));
      break;
    case MULTI1_eq_STEM:
      stk.push(make_tuple(i, j, _stem_best[j][i]));
      break;
    case OUTER_eq_OUTER_plus_UNPAIRED:
      k = j - 1;
      if (k != -1){
	stk.push(make_tuple(0, k, _outer_best[k]));
      }
      break;
    case OUTER_eq_OUTER_plus_STEM:
      {
	k = state.split;
	if (k != -1) {
	  stk.push(make_tuple(0, k, _outer_best[k]));
	  stk.push(make_tuple(k+1, j, _stem_best[j][k+1]));
	}else {
	  stk.push(make_tuple(i, j, _stem_best[j][i]));
	}
      }
      break;
    case NONE:
      cout << "TraceBack error." << endl;
      exit(1);
    }
  }
  
  return;
}

double LinAliFold::AllowBasePairs(int h, int i){  
  if(_allowed_base_pair[h][i-h] != 0.0){
    return(_allowed_base_pair[h][i-h]);
  }
  vector<int> type_freq(8,0);
  
  for(int s = 0; s < _num_of_seq; s++){
    int index_h = _a2s_index[s][h];
    int index_i = _a2s_index[s][i];    
    
    if(index_h >= 0 && index_i >= 0){
      int nuc_h = _int_seq_list[s][index_h];
      int nuc_i = _int_seq_list[s][index_i];
      int type = BP_pair[nuc_h][nuc_i];
      type_freq[type]++;
    }else{
      if(index_h < 0 && index_i < 0){
	type_freq[7]++;
      }else{
	type_freq[0]++;
      }
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
    score = score/_num_of_seq; //average
    
  }else{
    score = NEG_INF;
  }
  _allowed_base_pair[h][i-h] = score;    
  return(score);
}

int LinAliFold::DangleEnergy(int type, int a, int b, int seq_id){
  int x = 0;
  int len = _int_seq_list[seq_id].size();
  
  if (a >= 0 && _int_seq_list[seq_id][a] > 0) x += dangle5_37[type][_int_seq_list[seq_id][a]];
  if (b < len && _int_seq_list[seq_id][b] > 0) x += dangle3_37[type][_int_seq_list[seq_id][b]];
  if(type>2 || type == 0){
    x += TerminalAU;
  }
  return(x);
}

int LinAliFold::HairpinEnergy(int i, int j, int nuc_i, int nuc_i_p1, int nuc_j_m1, int nuc_j){
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

int LinAliFold::LoopEnergy(int type, int type2, int h, int i, int p, int q, int seq_id){
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
	if (type>2 || type == 0){ z += TerminalAU;}
	if (type2>2 || type2 == 0){ z += TerminalAU;}
      }
    }else{
      vector<int>& s = _int_seq_list[seq_id];
      if(type == 0){ type = 7;}
      if(type2 == 0){ type2 = 7;}
      
      if (u1+u2==2) {
	z = int11_37[type][type2][s[p+1]][s[q-1]];
      }else if ((u1==1) && (u2==2)){
	z = int21_37[type][type2][s[p+1]][s[i+1]][s[q-1]];
      }else if ((u1==2) && (u2==1)){
	z = int21_37[type2][type][s[i+1]][s[p+1]][s[h-1]];
      }else if ((u1==2) && (u2==2)){
	z = int22_37[type][type2][s[p+1]][s[h-1]][s[i+1]][s[q-1]];
      }else{
	if(type == 7){type = 0;}
	if(type2 == 7){type2 = 0;}
	z = internal_loop37[u1+u2]+mismatchI37[type][s[p+1]][s[q-1]]+mismatchI37[type2][s[i+1]][s[h-1]];
	z += min(MAX_NINIO, abs(u1-u2)*F_ninio37);
      }
    }
  }
  return z;
}

void LinAliFold::UpdateIfBetter(State& s, double new_score, Manner m){
  
  if (s.score < new_score){
    s.set(new_score, m);
  }
}

void LinAliFold::UpdateIfBetter(State &s, double new_score, Manner m, int split){
  if (s.score < new_score || s.type == NONE){
    s.set(new_score, m, split);
  }
};

void LinAliFold::UpdateIfBetter(State &s, double new_score, Manner m, int l1, int l2){
  if (s.score < new_score || s.type == NONE){
    s.set(new_score, m, l1, l2);
  }
};

void LinAliFold::Initialize(){
  _stem_best.resize(_alignment_length);
  _multi2_best.resize(_alignment_length);
  _multi1_best.resize(_alignment_length);
  _multi_best.resize(_alignment_length);
  _outer_best.resize(_alignment_length);
  _sorted_multi1_best.resize(_alignment_length);

  _allowed_base_pair.resize(_alignment_length);
  for(int i = 0; i < _alignment_length; i++){
    _allowed_base_pair[i].resize(_alignment_length-i);
  }
}

int LinAliFold::QuickSelectPartition(vector<pair<int, int> >& scores, int lower, int upper){
  int pivot = scores[upper].first;
  while (lower < upper) {
    while (scores[lower].first < pivot) ++lower;
    while (scores[upper].first > pivot) --upper;
    if (scores[lower].first == scores[upper].first) ++lower;
    else if (lower < upper) swap(scores[lower], scores[upper]);
    
  }
  return upper;
}

int LinAliFold::QuickSelect(vector<pair<int, int> >& scores, int lower,int upper, int k){
  if ( lower == upper ) return scores[lower].first;
  int split = QuickSelectPartition(scores, lower, upper);
  int length = split - lower + 1;
  if (length == k) return scores[split].first;
  else if (k  < length) return QuickSelect(scores, lower, split-1, k);
  else return QuickSelect(scores, split+1, upper, k - length);
}

void LinAliFold::SortMulti1(int threshold, unordered_map<int, State> &candidate, vector<pair<int, int> > &sorted_multi1){
  sorted_multi1.clear();
  
  for (auto &it : candidate) {
    int i = it.first;
    State &cand = it.second;
    int new_score = (i >= 1 ? _outer_best[i-1].score : 0) + cand.score;
    sorted_multi1.push_back(make_pair(new_score, i));
  }
  sort(sorted_multi1.begin(), sorted_multi1.end(), greater<pair<int, int> >());
}

int LinAliFold::PruneBeam(unordered_map<int, State> &candidate_list){
  vector<pair<int, int> > scores;
  
  for (auto &it : candidate_list) {
    int i = it.first;
    State &cand = it.second;
    int new_score;

    if ((i >= 1) && (_outer_best[i-1].score == MIN_VALUE)){
      new_score = MIN_VALUE;
    }else{
      new_score = (i >= 1 ? _outer_best[i-1].score : 0) + cand.score;
    }
    scores.push_back(make_pair(new_score, i));
  }
  if (scores.size() <= _beam_size){
    return MIN_VALUE;
  }
  int threshold = QuickSelect(scores, 0, scores.size() - 1, scores.size() - _beam_size);
  for (auto &p : scores) {
    if (p.first < threshold){
      candidate_list.erase(p.second);
    }
  }
  
  return threshold;
}

void LinAliFold::ReadData(){
  ifstream fp;
  string buffer;
  fp.open(_input_file_name.c_str(),ios::in);
  if (!fp){
    cout << "Error: can't open input_file:"+_input_file_name+"." <<endl;
    exit(1);
  }

  getline(fp,buffer);
  _name_list.push_back(buffer.substr(1,buffer.size()-1));
  
  string temp_seq = "";
  int count = 0;
  while(getline(fp,buffer)){
    if(buffer[0] == '>'){
      _name_list.push_back(buffer.substr(1,buffer.size()-1));      
      transform(temp_seq.begin(), temp_seq.end(), temp_seq.begin(), ::toupper);
      replace(temp_seq.begin(), temp_seq.end(), 'T', 'U');      
      _gapped_seq_list.push_back(temp_seq);
      temp_seq = "";
    }else{
      if(buffer.size()>=2){
	if(buffer.substr(buffer.size()-2,2) == "\r\n"){
	  buffer.erase(buffer.size()-2,2);
	}
      }
      if(buffer[buffer.size()-1] == '\r' || buffer[buffer.size()-1] == '\n'){
	buffer.erase(buffer.size()-1,1);
      }
      count += 1;
      temp_seq += buffer;
    }
  }
  transform(temp_seq.begin(), temp_seq.end(), temp_seq.begin(), ::toupper);
  replace(temp_seq.begin(), temp_seq.end(), 'T', 'U');
  _gapped_seq_list.push_back(temp_seq);
  fp.close();
  
  _num_of_seq = _gapped_seq_list.size();
  _alignment_length = _gapped_seq_list[0].size();
  _structure.resize(_alignment_length, '.');;
  _a2s_index.resize(_num_of_seq, vector<int>(_alignment_length,0));

  _int_seq_list.resize(_num_of_seq);
   
  for(int i = 0; i < _num_of_seq; i++){
    if(_alignment_length != _gapped_seq_list[i].size()){
      cout << "The alignment length is not identical for all sequences." << endl;
      exit(1);
    }
    
    int index = 0;
    vector<int> int_seq;
    
    for(int j = 0; j < _alignment_length; j++){
      char c = _gapped_seq_list[i][j];
      if((c != '-') && (c != '_') && (c != '~') && (c != '.')){
	int_seq.push_back(GET_ACGU_NUM(c));
	_a2s_index[i][j] = index;
	index++;	
      }else{
	_a2s_index[i][j] = -1;
      }
    }
    _int_seq_list[i] = int_seq;
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

void LinAliFold::SetParameters(int argc,char* argv[]){
  int c;
  extern char *optarg;
  while ((c = getopt(argc, argv, "i:b:d:e:r:")) != -1) {
    
    switch (c) {
    case 'i':
      _input_file_name = optarg;
      break;

    case 'b':
      _beam_size = atoi(optarg);
      if(_beam_size <= 0){
	cout << "The beam size have to be greater than zero." << endl;
	exit(1);
      }
      break;

    case 'd':
      _delta = atof(optarg);
      break;

    case 'e':
      _beta = atof(optarg);
      break;

    case 'r':
      _ribosum_flag = atoi(optarg);
      break;

    default:
      cerr << "The argument is an invalid command." << endl;
      exit(1); 
    }      
  }
}
