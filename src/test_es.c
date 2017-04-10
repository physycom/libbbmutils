#include <stdio.h>
#include <math.h>

#include "libbbmutils/bbmutils.h"

int main() {
  MAT2D m;

  printf("---- TEST 2D EIGENSYSTEM MATH_LIB ----\n");

  set_mat2d(&m, 0.0,  1.0,
               -2.0, -3.0);

  EigenSys es = eigs_2x2_sym(m);

  print_eigs(es, "non symmetric matrix");

  return 0;
}
