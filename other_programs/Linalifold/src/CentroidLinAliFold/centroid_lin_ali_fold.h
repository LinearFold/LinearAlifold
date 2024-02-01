#ifndef CENTROID_LIN_ALI_FOLD_H
#define CENTROID_LIN_ALI_FOLD_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <getopt.h>
#include <limits>
#include <sys/time.h>

#include "energy_par.h"
#include "intloops.h"
#include "utils.h"


using namespace std;

class CentroidLinAliFold{
public:
  CentroidLinAliFold(){
    _input_file_name = "";
    _beam_size = 100;
    _gamma = 1.0;
    _beta = 1.0;
    _delta  = 1.0;
    _mix_weight = 0.5;
    _ribosum_flag = 0;
    _pea_flag = 0;
    _bpp_output_flag = 0;
    _bpp_threshold = 0.0001;
    for(int i = -5; i<=10; i++){
      _gamma_array.push_back(pow(2,i));
      if(i==2){
	_gamma_array.push_back(6);
      }
    }
  }
  void Run();
  void SetParameters(int argc, char* argv[]);
  
private:
  void ReadData();
  void Predict(vector<unordered_map<int, int> >& back_pointer);
  string TraceBack(const int i, const int j, const vector<unordered_map<int, int> >& back_pointer, int& num_of_bp, double& expected_tp);
  
private:
  string _input_file_name;
  int _bpp_output_flag;
  int _ribosum_flag;
  int _pea_flag;
  int _beam_size;
  int _num_of_seq;
  int _alignment_length;
  float _bpp_threshold;
  double _delta;
  double _beta;
  double _gamma;
  double _mix_weight;
  vector<double> _gamma_array;
  vector<vector<int> > _a2s_index;
  vector<vector<int> > _s2a_index;
  unordered_map<pair<int,int>, float, hash_pair> _bpp_matrix;

  vector<string> _gapped_seq_list;
  vector<vector<int> > _int_seq_list;
};

#endif
