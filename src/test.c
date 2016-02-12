#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif

#include "matrix_lib.h"

int main(){
	VEC3D v1, v2, axis;
	MAT3D rot, id;
	double angle;

	printf("---- 3D ROTATIONS LIB ----\n");

	v1 = set_vec3d(1.0,1.0,0.0);
	print_vec3d(v1, "v1");
	v2 = set_vec3d(0.0,1.0,1.0);
	print_vec3d(v2, "v2");
	id = set_mat3d(1.0,0.0,0.0,0.0,2.0,0.0,0.0,0.0,3.0);
	print_mat3d(id, "id");

	angle = M_PI / 2;
	printf("angle = %5.2f\n", angle*180.0/M_PI);
	axis = set_vec3d(1.0,0.0,0.0);
	print_vec3d(axis, "n");
	rot = make_rotation(axis, angle);
  print_mat3d(rot, "rot");
  print_mat3d(transpose_mat3d(rot), "test transpose");


	print_vec3d(rotate_vec3d(v1, rot), "v1_rotated");
	print_vec3d(rotate_vec3d(v2, rot), "v2_rotated");
	print_mat3d(rotate_mat3d(id, rot), "id_rotated");
	print_mat3d(product_mat3d(rot, product_mat3d(id, transpose_mat3d(rot))), "test product");


	return 0;
}
