#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include "math_lib.h"

int main(){
  VEC6D v;
  MAT6D m;
  MAT3D rot1, rot2;

	printf("---- TEST 6D MATH_LIB ----\n");

  v = set_vec6d(1.0,1.0,0.0,  1.0,0.0,1.0);

  rot1 = make_rotation(set_vec3d(0.0, 1.0, 0.0), M_PI);

  m = set_mat6d( set_mat3d(  1.0, 0.0, 0.0, 0.0,  2.0, 0.0, 0.0, 0.0,  3.0 ) ,
                 set_mat3d(  4.0, 0.0, 0.0, 0.0,  5.0, 0.0, 0.0, 0.0,  6.0 ),
                 set_mat3d(  7.0, 0.0, 0.0, 0.0,  8.0, 0.0, 0.0, 0.0,  9.0 ),
                 set_mat3d( 10.0, 0.0, 0.0, 0.0, 11.0, 0.0, 0.0, 0.0, 12.0 ));

  rot2 = make_rotation(set_vec3d(0.0, 1.0, 0.0), M_PI/2.0);

  print_vec6d(v, "v");
  print_vec6d(rotate_vec6d(rot1, v), "v_rotated");
  print_mat6d(m, "\nm");
  print_mat6d(rotate_mat6d(rot2, m), "m_rotated");

  
  return 0;
}
