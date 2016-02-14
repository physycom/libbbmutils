#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include "math_lib.h"

int main() {
  VEC3D v, axis;
  MAT3D rot, m;
  double angle;

  printf("---- TEST 3D MATH_LIB ----\n");

  v = set_vec3d(1.0, 1.0, 0.0);

  m = set_mat3d(1.0, 0.0, 0.0, 
                0.0, 2.0, 0.0, 
                0.0, 0.0, 3.0);

  angle = M_PI / 2;
  axis = set_vec3d(1.0, 0.0, 0.0);
  rot = make_rotation(axis, angle);

  printf("angle = %5.2f\n", angle*180.0 / M_PI);
  print_vec3d(axis, "n");
  print_mat3d(rot, "rot");
  print_mat3d(transpose_mat3d(rot), "rot (test transpose)");
  print_vec3d(v, "\nv");
  print_vec3d(rotate_vec3d(rot, v), "v_rotated");
  print_mat3d(m, "\nm");
  print_mat3d(rotate_mat3d(rot, m), "m_rotated");
  print_mat3d(product_mat3d(product_mat3d(rot, m), transpose_mat3d(rot)), "m_rotated (test product)");


  return 0;
}
