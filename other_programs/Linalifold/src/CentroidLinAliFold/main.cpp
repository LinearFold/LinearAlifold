#include "centroid_lin_ali_fold.h"

void PrintUsage() {
  cout << "CentroidLinAliFold version 0.1.0 - Linear-time consensu structure predictin using generalized centroid estimators for RNA alignments" << endl;
}

int main(int argc, char** argv){
  if(argc == 1 || strcmp(argv[1],"-h") == 0){
    PrintUsage();
    exit(1);
  }
  CentroidLinAliFold centroid_lin_ali_fold;
  centroid_lin_ali_fold.SetParameters(argc, argv);
  centroid_lin_ali_fold.Run();
}
