#include "centroid_lin_ali_fold.h"
#include "linear_partition.h"
#include "lin_ali_partition.h"
#include <time.h>
#include <sstream>

void CentroidLinAliFold::Run(){
  // for runtime statistics
  struct timeval parse_starttime, parse_endtime;
  gettimeofday(&parse_starttime, NULL);

  ReadData();
  LinAliPartition lin_ali_partition(_int_seq_list, _beam_size, _ribosum_flag, _beta, _delta, _bpp_threshold, _a2s_index, _gapped_seq_list);
  unordered_map<pair<int,int>, float, hash_pair> lin_ali_bpp_matrix;
  lin_ali_partition.Run(lin_ali_bpp_matrix);
  
  for(auto &it : lin_ali_bpp_matrix){
    pair<int, int> temp_pair = make_pair(it.first.first, it.first.second);
    float prob = it.second*(1-_mix_weight);   
    if(_bpp_matrix.find(temp_pair) == _bpp_matrix.end()){
      _bpp_matrix[temp_pair] = prob;
    }else{
      _bpp_matrix[temp_pair] += prob;
    }
  }
  lin_ali_bpp_matrix.clear();
  
  for(int i = 0; i < _num_of_seq; i++){
    LinearPartition linear_partition(_int_seq_list[i], _beam_size, _bpp_threshold);
    unordered_map<pair<int,int>, float, hash_pair> temp_bpp_matrix;
    linear_partition.Run(temp_bpp_matrix);
    for(auto &it : temp_bpp_matrix){
      int p = it.first.first;
      int q = it.first.second;
      int new_p = _s2a_index[i][p];
      int new_q = _s2a_index[i][q];

      float prob = it.second*_mix_weight/_num_of_seq;

      pair<int, int> temp_pair = make_pair(new_p, new_q);
      if(_bpp_matrix.find(temp_pair) == _bpp_matrix.end()){
	_bpp_matrix[temp_pair] = prob;
      }else{
	_bpp_matrix[temp_pair] += prob;
      }
    }
    temp_bpp_matrix.clear();
  }

  if(_bpp_output_flag == 1){
    for(auto it = _bpp_matrix.begin(); it != _bpp_matrix.end();++it){
      auto i = it->first.first;
      auto j = it->first.second;
      auto score = it->second;
    
      if(score > _bpp_threshold){
	cout << i << " " << j << " " << score << endl;	
      }
    }
  }
  
  if(_gamma != -1 && _pea_flag == 0){
    _bpp_threshold = 1/(_gamma+1);
    int num_of_bp = 0;
    double etp = 0;
      
    vector<unordered_map<int, int> > back_pointer(_alignment_length);
    Predict(back_pointer);
    string structure = TraceBack(0, _alignment_length-1, back_pointer, num_of_bp, etp);
    cout << structure << endl;
  }else{
    double bpp_sum = 0.0;
    if(_pea_flag == 1){
      for(auto& it : _bpp_matrix){
	auto score = it.second;
	bpp_sum += score;
      }      
    }

    vector<string> string_array;
    double max_pmcc = 0.0;
    double temp_gamma = 0.0;
    string temp_structure = "";
    for(int i = _gamma_array.size()-1 ; i >= 0 ; i--){
      _gamma = _gamma_array[i];
      _bpp_threshold = 1/(_gamma+1);
      int num_of_bp = 0;
      double etp = 0;
      
      vector<unordered_map<int, int> > back_pointer(_alignment_length);
      Predict(back_pointer);
      stringstream ss;
      string structure = TraceBack(0, _alignment_length-1, back_pointer, num_of_bp, etp);

      if(_pea_flag == 0){
	ss << structure << " (g=" << _gamma << ")";
	string_array.push_back(ss.str());
      }else{
	double etn = _alignment_length*(_alignment_length-1)/2 - num_of_bp - bpp_sum + etp;
	double efp = num_of_bp - etp;
	double efn = bpp_sum - etp;

	double pseudo_mcc = etp*etn - efp*efn;
	pseudo_mcc /= sqrt((etp+efp)*(etp+efn)*(etn+efp)*(etn+efn));
	if(pseudo_mcc > max_pmcc){
	  max_pmcc = pseudo_mcc;
	  temp_gamma = _gamma;
	  temp_structure = structure;
	}
      }
    }
    if(_pea_flag == 0){
      for(int i = string_array.size()-1; i>=0 ; i--){
	cout << string_array[i] << endl;
      }
    }else{
      cout << temp_structure << " (g=" << temp_gamma << ", pmcc=" << max_pmcc << ")" << endl;
    }
  }

  gettimeofday(&parse_endtime, NULL);
  double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec - parse_starttime.tv_usec) / 1000000.0;
  fprintf(stdout, "CentroidLinAliFold Calculation Time: %.3f seconds.\n", parse_elapsed_time);
}

string CentroidLinAliFold::TraceBack(const int i, const int j, const vector<unordered_map<int, int> >& back_pointer, int& num_of_bp, double& expected_tp){
  if(i>j){
    return "";
  }
  auto it = back_pointer[i].find(j);
  if(it == back_pointer[i].end()){
    if (i == j){
      return ".";
    }else{
      return "." + TraceBack(i+1,j, back_pointer, num_of_bp, expected_tp);
    }
  }else{
    int k = it->second;
    string temp;
    if (k == j){
      temp = "";
    }else{
      temp = TraceBack(k+1,j, back_pointer, num_of_bp, expected_tp);
    }
    num_of_bp++;
    expected_tp += _bpp_matrix[make_pair(i,k)];
    return "(" + TraceBack(i+1,k-1, back_pointer, num_of_bp, expected_tp) + ")" + temp;
  }
  return "";
}

void CentroidLinAliFold::Predict(vector<unordered_map<int, int> >& back_pointer){
  vector<vector<float> > dp_matrix(_alignment_length, vector<float>(_alignment_length, 0.0));
  vector<unordered_map<int, float> > bpp(_alignment_length);  
  vector<vector<int> > paired_base(_alignment_length);

  for(auto it = _bpp_matrix.begin(); it != _bpp_matrix.end();){
    auto i = it->first.first;
    auto j = it->first.second;
    auto score = it->second;
    
    if(score > _bpp_threshold){
      bpp[i][j] = score;    
      paired_base[i].push_back(j);
      ++it;
    }else{      
      it = _bpp_matrix.erase(it);
    }
  }
  for (int i = 0; i < _alignment_length; ++i){
    sort(paired_base[i].begin(), paired_base[i].end());
  }

  for (int l = TURN+1; l< _alignment_length; l++){
    for (int i = 0; i< _alignment_length - l; i++){	  
      int j = i + l;
      
      dp_matrix[i][j] = dp_matrix[i+1][j];
      for (int k : paired_base[i]){
	if (k>j){
	  break;
	}
	float score_kp1_j = 0.0;
	if (k<j){
	  score_kp1_j = dp_matrix[k+1][j];
	}
	float temp_score = (_gamma+1) * bpp[i][k] - 1 + dp_matrix[i+1][k-1] + score_kp1_j;
	if (dp_matrix[i][j] < temp_score){
	  dp_matrix[i][j] = temp_score;
	  back_pointer[i][j] = k;
	}
      }
    }
  }    
}

void CentroidLinAliFold::ReadData(){
  ifstream fp;
  string buffer;
  fp.open(_input_file_name.c_str(),ios::in);
  if (!fp){
    cout << "Error: can't open input_file:"+_input_file_name+"." <<endl;
    exit(1);
  }

  getline(fp,buffer);
  
  string temp_seq = "";
  int count = 0;
  while(getline(fp,buffer)){
    if(buffer[0] == '>'){
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
  _a2s_index.resize(_num_of_seq, vector<int>(_alignment_length,0));
  _int_seq_list.resize(_num_of_seq);
   
  for(int i = 0; i < _num_of_seq; i++){
    if(_alignment_length != _gapped_seq_list[i].size()){
      cout << "The alignment length is not identical for all sequences." << endl;
      exit(1);
    }
    
    int index = 0;
    vector<int> int_seq;
    vector<int> s2a;
    
    for(int j = 0; j < _alignment_length; j++){
      char c = _gapped_seq_list[i][j];
      if((c != '-') && (c != '_') && (c != '~') && (c != '.')){
	int_seq.push_back(GET_ACGU_NUM(c));
	_a2s_index[i][j] = index;
	s2a.push_back(j);
	index++;
      }else{
	_a2s_index[i][j] = -1;
      }
    }
 
    _s2a_index.push_back(s2a);
    _int_seq_list[i] = int_seq;
  }
}

void CentroidLinAliFold::SetParameters(int argc,char* argv[]){
  int c;
  extern char *optarg;
  while ((c = getopt(argc, argv, "i:o:b:d:e:g:t:w:r:p:")) != -1) {
    
    switch (c) {
    case 'i':
      _input_file_name = optarg;
      break;

    case 'o':
      _bpp_output_flag = atoi(optarg);
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

    case 'g':
      _gamma = atof(optarg);
      if(_gamma < 0 && _gamma != -1){
	cout << "The gamma is an invalid value." << endl;
	exit(1);
      }
      break;

    case 't':
      _bpp_threshold = atof(optarg);
      break;

    case 'w':
      _mix_weight = atof(optarg);
      break;

    case 'r':
      _ribosum_flag = atoi(optarg);
      break;

    case 'p':
      _pea_flag = atoi(optarg);
      break;

    default:
      cerr << "The argument is an invalid command." << endl;
      exit(1); 
    }      
  }
}
