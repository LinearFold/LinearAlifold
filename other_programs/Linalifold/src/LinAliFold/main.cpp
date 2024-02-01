#include "linalifold.h"

void PrintUsage() {
  cout << "LinAliFold version 0.1.0 - Linear-time consensus structure prediction for RNA alignments." << endl;
}

int main(int argc, char** argv){
  if(argc == 1 || strcmp(argv[1],"-h") == 0){
    PrintUsage();
    exit(1);
  }
  
  LinAliFold lin_ali_fold;
  lin_ali_fold.SetParameters(argc, argv);
  lin_ali_fold.Run();
}
