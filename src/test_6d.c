#include <stdio.h>
#include <math.h>

#include "libbbmutils/bbmutils.h"

int main(){
  VEC6D v, vr;
  VEC3D axis;
  MAT3D m00, m01, m10, m11, rot;
  MAT6D m, mr;

	printf("---- TEST 6D MATH_LIB ----\n");

  set_vec6d(&v, 1.0,1.0,0.0,  1.0,0.0,1.0);
  print_vec6d(&v, "\n6d vector = v");

  set_mat3d(&m00,  1.0, 0.0, 0.0, 0.0,  2.0, 0.0, 0.0, 0.0,  3.0);
  set_mat3d(&m01,  4.0, 0.0, 0.0, 0.0,  5.0, 0.0, 0.0, 0.0,  6.0);
  set_mat3d(&m10,  7.0, 0.0, 0.0, 0.0,  8.0, 0.0, 0.0, 0.0,  9.0);
  set_mat3d(&m11, 10.0, 0.0, 0.0, 0.0, 11.0, 0.0, 0.0, 0.0, 12.0);
  set_mat6d(&m, &m00, &m01, &m10, &m11);
  print_mat6d(&m, "\n6d matrix = A");

  set_vec3d(&axis, 0.0, 1.0, 0.0);
  make_rotation(&rot, &axis, M_PI/2);
  print_mat3d(&rot, "\n3d matrix (rotation) = R");

  rotate_vec6d(&vr, &rot, &v);
  print_vec6d(&vr, "\n6d vector (test rotation) = Rv");

  rotate_mat6d(&mr, &rot, &m);
  print_mat6d(&mr, "\n6d matrix (test rotation) = R A Rt");
  
  return 0;
}
