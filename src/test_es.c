#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include "math_lib.h"

int main() {
  MAT2D m;

  printf("---- TEST 2D EIGENSYSTEM MATH_LIB ----\n");

  m = set_mat2d( 0.0,  1.0,
                -2.0, -3.0);

  EigenSys es = eigs_2x2_sym(m);

  print_eigs(es, "non symmetric matrix");

  return 0;
}
