#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include "math_lib.h"

int main() {
  VEC2D v;
  MAT2D rot, m;
  double angle;

  printf("---- TEST 2D MATH_LIB ----\n");

  v = set_vec2d(1.0, 1.0);

  m = set_mat2d(0.0, 5.0,
                5.0, -2.0);

  angle = M_PI / 4;
  rot = make_rotation_2d(angle);

  printf("angle = %5.2f\n", angle*180.0 / M_PI);
  print_mat2d(rot, "rotation = R");
  print_mat2d(transpose_mat2d(rot), "rotation (test transpose) = Rt");
  print_vec2d(v, "\n2d vector = v");
  print_vec2d(rotate_vec2d(rot, v), "2d vector (rotated) = Rv");
  print_mat2d(m, "\n2d matrix = A");
  print_mat2d(rotate_mat2d(rot, m), "2d matrix (rotated, test function) = f(A)");
  print_mat2d(product_mat2d(rot, product_mat2d(m, transpose_mat2d(rot))), "2d matrix (rotated, test product) = RARt");
  print_vec2d(rotate_vec2d(m, v), "\n2d vector = w = Av");
  print_vec2d(rotate_vec2d(rot, rotate_vec2d(m, v)), "\n2d vector (rotated) = Rw");
  print_vec2d(rotate_vec2d( rotate_mat2d(rot, m), rotate_vec2d(rot, v) ), "\n w (rotated, test associativity) = (RvRt) (Rv)");


  return 0;
}
