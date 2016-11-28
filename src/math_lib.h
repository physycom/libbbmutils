#define ENABLE_DISPLAY_FUNCTIONS

#ifdef __cplusplus
extern "C" {
#endif

  // 2d vector algebra
  typedef struct VEC2D {
    double x, y, mod;
  } VEC2D;

  void set_vec2d(VEC2D * v, const double x, const double y);

  void setmod_vec2d(VEC2D * v);

  void normalize_vec2d(VEC2D * v);

  double prod_dot_2d(const VEC2D * v, const VEC2D * u);

  double prod_cross_2d(const VEC2D * v, const VEC2D * u);            // returns z-component on a right-handed frame

  void multiply_vec2d(VEC2D * v, const double alpha);


  // 2d matrix algebra
  typedef struct MAT2D {
    double xx, xy, yx, yy;
  } MAT2D;

  void set_mat2d(MAT2D * m, double xx, double xy, double yx, double yy);
  
  void transpose_mat2d(MAT2D * mt, const MAT2D * m);

  void product_mat2d(MAT2D * result, const MAT2D * a, const MAT2D * b);

  void multiply_mat2d(double, MAT2D *);

  void make_rotation_2d_cs(MAT2D * rot, const double c, const double s);

  void make_rotation_2d(MAT2D * rot, const double angle);

  void rotate_vec2d(VEC2D * result, const MAT2D * rot, const VEC2D * v);

  void rotate_mat2d(MAT2D * result, const MAT2D * rot, const MAT2D * m);


  // 3d vector algebra
  typedef struct VEC3D {
    double x, y, z, mod;
  } VEC3D;

  void set_vec3d(VEC3D * v, const double x, const double y, const double z);

  void setmod_vec3d(VEC3D * v);

  void normalize_vec3d(VEC3D * v);

  double prod_dot_3d(const VEC3D * a, const VEC3D * b);

  void prod_cross_3d(VEC3D * result, const VEC3D * a, const VEC3D * b);

  void multiply_vec3d(VEC3D * v, const double alpha);


  // 3d matrix algebra
  typedef struct MAT3D {
    double xx, xy, xz, yx, yy, yz, zx, zy, zz;
  } MAT3D;

  void set_mat3d(MAT3D * m, const double xx, const double xy, const double xz, const double yx, const double yy, const double yz, const double zx, const double zy, const double zz);

  void transpose_mat3d(MAT3D * result, const MAT3D * m);

  void product_mat3d(MAT3D * result, const MAT3D * a, const MAT3D * b);

  void multiply_mat3d(MAT3D * m, const double alpha);

  void make_rotation_cs(MAT3D * rot, const VEC3D * v, const double c, const double s);

  void make_rotation(MAT3D * rot, const VEC3D * v, const double angle);

  void rotate_vec3d(VEC3D * result, const MAT3D * rot, const VEC3D * v);

  void rotate_mat3d(MAT3D * result, const MAT3D * rot, const MAT3D * m);


  // 6d vector
  typedef struct VEC6D {
    VEC3D a, g;
  } VEC6D;

  void set_vec6d(VEC6D * v, const double ax, const double ay, const double az, const double gx, const double gy, const double gz);

  void rotate_vec6d(VEC6D * result, const MAT3D * rot, const VEC6D * v);


  // 6d matrix
  typedef struct MAT6D {
    MAT3D A[2][2];
  } MAT6D;

  void set_mat6d(MAT6D * m, const MAT3D * A00, const MAT3D * A01, const MAT3D * A10, const MAT3D * A11);

  void rotate_mat6d(MAT6D * result, const MAT3D * rot, const MAT6D * m);


  // 2d eigenvalue problem
  typedef struct EigenSys {
    MAT2D A;
    double l1, l2;
    VEC2D v1, v2, u1, u2;
  } EigenSys;

  EigenSys eigs_2x2_sym_normalized(double, double);

  EigenSys eigs_2x2_sym(MAT2D);

  int check_eigs(EigenSys);

#ifdef ENABLE_DISPLAY_FUNCTIONS
  // display function
  void print_vec2d(const VEC2D * v, const char * name);
  
  void print_mat2d(const MAT2D * m, const char * name);

  void print_vec3d(const VEC3D * v, const char * name);

  void print_mat3d(const MAT3D * m, const char * name);

  void print_vec6d(const VEC6D * v, const char * name);

  void print_mat6d(const MAT6D * m, const char * name);

  void print_eigs(EigenSys, const char *);
#endif

#ifdef __cplusplus
}
#endif
