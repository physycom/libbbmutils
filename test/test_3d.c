#include <stdio.h>
#include <math.h>

#include <bbmutils.h>

int main() {
  VEC3D v, u, vr, ur, axis;
  MAT3D m, mt, m2, mr, rot;
  double angle;

  printf("---- TEST 3D MATH_LIB ----\n");

  set_vec3d(&v, 1.0, 1.0, 0.0);
  print_vec3d(&v, "\n3d vector = v");

  set_mat3d(&m, 1.0, 0.0, -1.0,
                0.0, 2.0,  0.0,
                0.0, 0.0,  3.0);
  print_mat3d(&m, "\n3d matrix = A");

  rotate_vec3d(&u, &m, &v);
  print_vec3d(&u, "\n3d vector (test product) = u = Av");

  transpose_mat3d(&mt, &m);
  print_mat3d(&mt, "\n3d matrix (test transpose) = At");

  product_mat3d(&m2, &m, &m);
  print_mat3d(&m2, "\n3d matrix (test product) = A A");

  angle = M_PI / 2;
  set_vec3d(&axis, 1.0, 0.0, 0.0);
  make_rotation(&rot, &axis, angle);
  print_mat3d(&rot, "\n3d matrix (rotation) = R");

  rotate_vec3d(&vr, &rot, &v);
  print_vec3d(&vr, "\n3d vector (test rotation) = Rv");

  rotate_mat3d(&mr, &rot, &m);
  print_mat3d(&mr, "\n3d matrix (test rotation) = B = R A Rt");

  rotate_vec3d(&ur, &rot, &u);
  print_vec3d(&ur, "\n3d vector = Ru");
  rotate_vec3d(&ur, &mr, &vr);
  print_vec3d(&ur, "\n3d vector = Bu");

  return 0;
}
