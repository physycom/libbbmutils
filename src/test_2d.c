#include <stdio.h>
#include <math.h>

#include "libbbmutils/bbmutils.h"

int main() {
  printf("---- TEST 2D MATH_LIB ----\n");
  double angle;
  VEC2D v, vr, u, ur;
  MAT2D m, mt, m2, rot, mr;

  set_vec2d(&v, 1.0,  1.0);
  print_vec2d(&v, "\n2d vector = v");
  
  set_mat2d(&m, 0.0,  4.0,
                5.0, -2.0);
  print_mat2d(&m, "\n2d matrix = A");
  
  rotate_vec2d(&u, &m, &v);
  print_vec2d(&u, "\n2d vector (test product) = u = Av");

  transpose_mat2d(&mt, &m);
  print_mat2d(&mt, "\n2d matrix (test transpose) = At");

  product_mat2d(&m2, &m, &m);
  print_mat2d(&m2, "\n2d matrix (test product) = A A");
  
  angle = M_PI / 4;
  make_rotation_2d(&rot, angle);
  printf("\nangle = %5.2f\n", angle*180.0 / M_PI);
  print_mat2d(&rot, "\n2d matrix (rotation) = R");

  rotate_vec2d(&vr, &rot, &v);
  print_vec2d(&vr, "\n2d vector (test rotation) = Rv");

  rotate_mat2d(&mr, &rot, &m);
  print_mat2d(&mr, "\n2d matrix (test rotation) = B = R A Rt");
  
  rotate_vec2d(&ur, &rot, &u);
  print_vec2d(&ur, "\n2d vector = Ru");
  rotate_vec2d(&ur, &mr, &vr);
  print_vec2d(&ur, "\n2d vector = Bu");

  return 0;
}
