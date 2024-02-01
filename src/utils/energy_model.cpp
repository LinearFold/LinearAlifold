#include "energy_model.h"

char Triloops[241] =
    "CAACG "
    "GUUAC ";

char Tetraloops[281] =
    "CAACGG "
    "CCAAGG "
    "CCACGG "
    "CCCAGG "
    "CCGAGG "
    "CCGCGG "
    "CCUAGG "
    "CCUCGG "
    "CUAAGG "
    "CUACGG "
    "CUCAGG "
    "CUCCGG "
    "CUGCGG "
    "CUUAGG "
    "CUUCGG "
    "CUUUGG ";

char Hexaloops[361] =
    "ACAGUACU "
    "ACAGUGAU "
    "ACAGUGCU "
    "ACAGUGUU ";

#ifdef EM_Vienna
bool mismatch1nI37_exists = true;
bool mismatch23I37_exists = true;
bool mismatchM37_exists = true;
bool mismatchExt37_exists = true;
bool specialHP_exists = true;
#endif

// For bl* energy model, these parameters are not defined
#ifdef EM_BL_Star
bool mismatch1nI37_exists = false;
bool mismatch23I37_exists = false;
bool mismatchM37_exists = false;
bool mismatchExt37_exists = false;
bool specialHP_exists = false;

int Triloop37[];
int Tetraloop37[];
int Hexaloop37[];

int mismatch1nI37[8][5][5];
int mismatch23I37[8][5][5];
int mismatchM37[8][5][5];
int mismatchExt37[8][5][5];
#endif
