/**
 * @file   math_lib.c
 * @Author A. Fabbri (alessandro.fabbri27@unibo.it), S. Sinigardi (stefano.sinigardi@unibo.it)
 * @date   February, 2016
 * @brief  This source file contains the definitions of functions for the mathematical framework, as well as some preprocessor macros.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>     // memset
#include <math.h>

#include "math_lib/math_lib.h"

/*! \brief \f$ 180/\pi \f$, to convert from radians to degree. */
#define RAD_TO_DEG    57.2957795131             /* 180/pi */

#ifndef M_PI
/*! \brief \f$ \pi \f$. */
#define M_PI          3.14159265358979323846    /* pi */
#endif

/*! \brief Zero threshold for numerically null vectors. */
#define EPSILON       1e-8

// 2d vector algebra
void set_vec2d(VEC2D * v, const double x, const double y) {
  v->x = x;
  v->y = y;
  v->mod = sqrt(x*x + y*y);
  if (v->mod < EPSILON) {
    v->x = 0.0;
    v->y = 0.0;
    v->mod = 0.0;
  }
}

void setmod_vec2d(VEC2D * v) {
  v->mod = sqrt(v->x*v->x + v->y*v->y);
  if (v->mod < EPSILON) v->mod = 0.0;
}

void normalize_vec2d(VEC2D * v) {
  v->x /= v->mod;
  v->y /= v->mod;
  v->mod = 1.0;
}

double prod_dot_2d(const VEC2D * v, const VEC2D * u) {
  return v->x * u->x + v->y * u->y;
}

double prod_cross_2d(const VEC2D * v, const VEC2D * u) {
  return v->x * u->y - v->y * u->x;
}

void multiply_vec2d(VEC2D * v, const double alpha) {
  v->x *= alpha;
  v->y *= alpha;
  v->mod *= fabs(alpha);
}


// 2d matrix algebra
void set_mat2d(MAT2D * m, double xx, double xy, double yx, double yy) {
  m->xx = xx;
  m->xy = xy;
  m->yx = yx;
  m->yy = yy;
}

void transpose_mat2d(MAT2D * mt, const MAT2D * m){
  mt->xx = m->xx;
  mt->xy = m->yx;
  mt->yx = m->xy;
  mt->yy = m->yy;
}

void product_mat2d(MAT2D * result, const MAT2D * a, const MAT2D * b) {
  result->xx = a->xx * b->xx + a->xy * b->yx;
  result->xy = a->xx * b->xy + a->xy * b->yy;
  result->yx = a->yx * b->xx + a->yy * b->yx;
  result->yy = a->yx * b->xy + a->yy * b->yy;
}

void multiply_mat2d(double alpha, MAT2D * A) {
  A->xx *= alpha;
  A->xy *= alpha;
  A->yx *= alpha;
  A->yy *= alpha;
}

inline void make_rotation_2d_cs(MAT2D * rot, const double c, const double s) {
  rot->xx = c;
  rot->xy = -s;
  rot->yx = s;
  rot->yy = c;
}

void make_rotation_2d(MAT2D * rot, const double angle) {
  double c = cos(angle), s = sin(angle);
  make_rotation_2d_cs(rot, c, s);
}

void rotate_vec2d(VEC2D * result, const MAT2D * rot, const VEC2D * v) {
  result->x = rot->xx * v->x + rot->xy * v->y;
  result->y = rot->yx * v->x + rot->yy * v->y;
}

void rotate_mat2d(MAT2D * result, const MAT2D * rot, const MAT2D * m) {
  MAT2D temp;
  transpose_mat2d(result, rot);        // result = rot^T
  product_mat2d(&temp, m, result);     // temp   = m rot^T
  product_mat2d(result, rot, &temp);   // result = rot m rot^T
}


// 3d vector algebra
void set_vec3d(VEC3D * v, const double x, const double y, const double z) {
  v->x = x;
  v->y = y;
  v->z = z;
  v->mod = sqrt(x*x + y*y + z*z);
  if (v->mod < EPSILON) {
    v->mod = 0.0;
    v->x = x;
    v->y = y;
    v->z = z;
  }
}

void setmod_vec3d(VEC3D * v) {
  v->mod = sqrt(v->x*v->x + v->y*v->y + v->z*v->z);
}

void normalize_vec3d(VEC3D * v) {
  v->x /= v->mod;
  v->y /= v->mod;
  v->z /= v->mod;
  v->mod = 1.0;
}

double prod_dot_3d(const VEC3D * v, const VEC3D * u) {
  return v->x * u->x + v->y * u->y + v->z * u->z;
}

void prod_cross_3d(VEC3D * result, const VEC3D * v, const VEC3D * u) {
    result->x = v->y * u->z - v->z * u->y;
    result->y = v->z * u->x - v->x * u->z;
    result->z = v->x * u->y - v->y * u->x;
}

void multiply_vec3d(VEC3D * v, const double alpha) {
  v->x *= alpha;
  v->y *= alpha;
  v->z *= alpha;
  v->mod *= fabs(alpha);
}


// 3d matrix algebra
void set_mat3d(MAT3D * m, const double xx, const double xy, const double xz, const double yx, const double yy, const double yz, const double zx, const double zy, const double zz) {
  m->xx = xx;
  m->xy = xy;
  m->xz = xz;
  m->yx = yx;
  m->yy = yy;
  m->yz = yz;
  m->zx = zx;
  m->zy = zy;
  m->zz = zz;
}

void transpose_mat3d(MAT3D * result, const MAT3D * m) {
  result->xx = m->xx;
  result->xy = m->yx;
  result->xz = m->zx;
  result->yx = m->xy;
  result->yy = m->yy;
  result->yz = m->zy;
  result->zx = m->xz;
  result->zy = m->yz;
  result->zz = m->zz;
}

void product_mat3d(MAT3D * result, const MAT3D * a, const MAT3D * b) {
  result->xx = a->xx * b->xx + a->xy * b->yx + a->xz * b->zx;
  result->xy = a->xx * b->xy + a->xy * b->yy + a->xz * b->zy;
  result->xz = a->xx * b->xz + a->xy * b->yz + a->xz * b->zz;
  result->yx = a->yx * b->xx + a->yy * b->yx + a->yz * b->zx;
  result->yy = a->yx * b->xy + a->yy * b->yy + a->yz * b->zy;
  result->yz = a->yx * b->xz + a->yy * b->yz + a->yz * b->zz;
  result->zx = a->zx * b->xx + a->zy * b->yx + a->zz * b->zx;
  result->zy = a->zx * b->xy + a->zy * b->yy + a->zz * b->zy;
  result->zz = a->zx * b->xz + a->zy * b->yz + a->zz * b->zz;
}

void multiply_mat3d(MAT3D * m, const double alpha) {
  m->xx *= alpha;
  m->xy *= alpha;
  m->xz *= alpha;
  m->yx *= alpha;
  m->yy *= alpha;
  m->yz *= alpha;
  m->zx *= alpha;
  m->zy *= alpha;
  m->zz *= alpha;
}

inline void make_rotation_cs(MAT3D * rot, const VEC3D * v, const double c, const double s) {
  VEC3D axis;
  set_vec3d(&axis, v->x, v->y, v->z);

  if (axis.mod < EPSILON) {
    memset(rot, 0, sizeof(*rot));
    return;
  }
  else normalize_vec3d(&axis);

  rot->xx = axis.x*axis.x*(1 - c) + c;   	    rot->xy = axis.x*axis.y*(1 - c) - axis.z*s; rot->xz = axis.x*axis.z*(1 - c) + axis.y*s;
  rot->yx = axis.y*axis.x*(1 - c) + axis.z*s; rot->yy = axis.y*axis.y*(1 - c) + c;	      rot->yz = axis.y*axis.z*(1 - c) - axis.x*s;
  rot->zx = axis.z*axis.x*(1 - c) - axis.y*s; rot->zy = axis.z*axis.y*(1 - c) + axis.x*s; rot->zz = axis.z*axis.z*(1 - c) + c;
}

void make_rotation(MAT3D * rot, const VEC3D * v, const double angle) {
  double c = cos(angle), s = sin(angle);
  make_rotation_cs(rot, v, c, s);
}

void rotate_vec3d(VEC3D * result, const MAT3D * rot, const VEC3D * v) {
  result->x = rot->xx * v->x + rot->xy * v->y + rot->xz * v->z;
  result->y = rot->yx * v->x + rot->yy * v->y + rot->yz * v->z;
  result->z = rot->zx * v->x + rot->zy * v->y + rot->zz * v->z;
}

void rotate_mat3d(MAT3D * result, const MAT3D * rot, const MAT3D * m) {
  MAT3D temp;
  transpose_mat3d(result, rot);        // result = rot^T
  product_mat3d(&temp, m, result);     // temp   = m rot^T
  product_mat3d(result, rot, &temp);   // result = rot m rot^T
}


// 6d vector
void set_vec6d(VEC6D * v, const double ax, const double ay, const double az, const double gx, const double gy, const double gz) {
  set_vec3d(&(v->a), ax, ay, az);
  set_vec3d(&(v->g), gx, gy, gz);
}

void rotate_vec6d(VEC6D * result, const MAT3D * rot, const VEC6D * v) {
  rotate_vec3d(&(result->a), rot, &(v->a));
  rotate_vec3d(&(result->g), rot, &(v->g));
}


// 6d matrix
void set_mat6d(MAT6D * m, const MAT3D * A00, const MAT3D * A01, const MAT3D * A10, const MAT3D * A11) {
  m->A[0][0] = *A00;
  m->A[0][1] = *A01;
  m->A[1][0] = *A10;
  m->A[1][1] = *A11;
}

void rotate_mat6d(MAT6D * result, const MAT3D * rot, const MAT6D * m) {
  rotate_mat3d(&(result->A[0][0]), rot, &(m->A[0][0]));
  rotate_mat3d(&(result->A[0][1]), rot, &(m->A[0][1]));
  rotate_mat3d(&(result->A[1][0]), rot, &(m->A[1][0]));
  rotate_mat3d(&(result->A[1][1]), rot, &(m->A[1][1]));
}


// 2d eigenvalue problem

// Returns the eigensystem of a 2x2 matrix of the form
// | 1 a |
// | a b |
// for real a,b
EigenSys eigs_2x2_sym_normalized(double a, double b) {
  EigenSys es;
  es.A.xx = 1.0;
  es.A.xy = a;
  es.A.yx = a;
  es.A.yy = b;
  es.l1 = (1 + b)*0.5 + 0.5*sqrt((b - 1)*(b - 1) + 4 * a*a);
  es.l2 = (1 + b)*0.5 - 0.5*sqrt((b - 1)*(b - 1) + 4 * a*a);
  set_vec2d(&(es.v1), 1.0, (es.l1 - 1) / a);
  set_vec2d(&(es.v2), 1.0, (es.l2 - 1) / a );
  if ( prod_cross_2d(&(es.v1), &(es.v2)) < 0) { // "right-handed", in the sense that z-axis is pointing "outward" the screen
    multiply_vec2d(&(es.v2), -1.0);
  }
  es.u1 = es.v1; normalize_vec2d(&es.u1);
  es.u2 = es.v2; normalize_vec2d(&es.u2);
  return es;
}

// Returns the eigensystem of a 2x2 matrix of the form
// | a b |
// | b c |
// for real a,b,c
EigenSys eigs_2x2_sym(MAT2D A) {
  double trA, detA;
  EigenSys es;
  set_mat2d(&(es.A), A.xx, A.xy, A.yx, A.yy);
  trA = es.A.xx + es.A.yy;
  detA = es.A.xx * es.A.yy - es.A.xy * es.A.yx;
  es.l1 = (trA + sqrt(trA*trA - 4 * detA))*0.5;
  es.l2 = (trA - sqrt(trA*trA - 4 * detA))*0.5;
  set_vec2d(&(es.v1), 1.0, (es.l1 - es.A.xx) / es.A.xy);         // A.xx * 1 + A.xy * y = l1 * 1
  set_vec2d(&(es.v2), 1.0, (es.l2 - es.A.xx) / es.A.xy);         // A.xx * 1 + A.xy * y = l2 * 1
  if (prod_cross_2d(&(es.v1), &(es.v2)) < 0) { // "right-handed", in the sense that z-axis is pointing "outward" the screen
    multiply_vec2d(&(es.v2), -1.0);
  }
  es.u1 = es.v1; normalize_vec2d(&es.u1);
  es.u2 = es.v2; normalize_vec2d(&es.u2);
  return es;
}

int check_eigs(EigenSys es) {
  int check = 0;
  // check eigenvalues
  VEC2D lhs, rhs;
  rotate_vec2d(&lhs, &(es.A), &(es.u1));
  rhs = es.u1; 
  multiply_vec2d(&rhs, es.l1);
  if (fabs(lhs.x - rhs.x) < EPSILON) check += 1;
  if (fabs(lhs.y - rhs.y) < EPSILON) check += 10;
  rotate_vec2d(&lhs, &(es.A), &(es.u2));
  rhs = es.u2; multiply_vec2d(&rhs, es.l2);
  if (fabs(lhs.x - rhs.x) < EPSILON) check += 100;
  if (fabs(lhs.y - rhs.y) < EPSILON) check += 1000;
  return check;
}

#ifdef ENABLE_DISPLAY_FUNCTIONS
// Display functions
void print_vec2d(const VEC2D * v, const char * name) {
  printf("%s = \n| %6.3f  %6.3f |\n", name, v->x, v->y);
}

void print_mat2d(const MAT2D * m, const char * name) {
  printf("%s = \n| %6.3f %6.3f |\n| %6.3f %6.3f |\n", name, m->xx, m->xy, m->yx, m->yy);
}

void print_vec3d(const VEC3D * v, const char * name) {
  printf("%s = \n| %6.3f  %6.3f  %6.3f |\n", name, v->x, v->y, v->z);
}

void print_mat3d(const MAT3D * m, const char * name) {
  printf("%s = \n| %6.3f %6.3f %6.3f |\n| %6.3f %6.3f %6.3f |\n| %6.3f %6.3f %6.3f |\n", name, m->xx, m->xy, m->xz, m->yx, m->yy, m->yz, m->zx, m->zy, m->zz);
}

void print_vec6d(const VEC6D * v, const char * name) {
  printf("%s = \n| %6.3f %6.3f %6.3f  ,  %6.3f %6.3f %6.3f |\n", name, v->a.x, v->a.y, v->a.z, v->g.x, v->g.y, v->g.z);
}

void print_mat6d(const MAT6D * m, const char * name) {
  printf("%s = \n", name);
  printf("| %6.3f %6.3f %6.3f \t %6.3f %6.3f %6.3f |\n", m->A[0][0].xx, m->A[0][0].xy, m->A[0][0].xz, m->A[0][1].xx, m->A[0][1].xy, m->A[0][1].xz);
  printf("| %6.3f %6.3f %6.3f \t %6.3f %6.3f %6.3f |\n", m->A[0][0].yx, m->A[0][0].yy, m->A[0][0].yz, m->A[0][1].yx, m->A[0][1].yy, m->A[0][1].yz);
  printf("| %6.3f %6.3f %6.3f \t %6.3f %6.3f %6.3f |\n", m->A[0][0].zx, m->A[0][0].zy, m->A[0][0].zz, m->A[0][1].zx, m->A[0][1].zy, m->A[0][1].zz);
  printf("| %8s                                    |\n", " ");
  printf("| %6.3f %6.3f %6.3f \t %6.3f %6.3f %6.3f |\n", m->A[1][0].xx, m->A[1][0].xy, m->A[1][0].xz, m->A[1][1].xx, m->A[1][1].xy, m->A[1][1].xz);
  printf("| %6.3f %6.3f %6.3f \t %6.3f %6.3f %6.3f |\n", m->A[1][0].yx, m->A[1][0].yy, m->A[1][0].yz, m->A[1][1].yx, m->A[1][1].yy, m->A[1][1].yz);
  printf("| %6.3f %6.3f %6.3f \t %6.3f %6.3f %6.3f |\n", m->A[1][0].zx, m->A[1][0].zy, m->A[1][0].zz, m->A[1][1].zx, m->A[1][1].zy, m->A[1][1].zz);
}

void print_eigs(EigenSys es, const char * description) {
  printf("%s\n", description);
  print_mat2d(&(es.A), "A ");
  printf("\nl1 = %6.3f\n", es.l1);
  print_vec2d(&(es.v1), "v1 ");
  print_vec2d(&(es.u1), "u1 ");
  printf("\nl2 = %6.3f\n", es.l2);
  print_vec2d(&(es.u2), "u2 ");
  print_vec2d(&(es.v2), "v2 ");
  printf("\ncheck (1111 is ok) = %04d\n", check_eigs(es));
}

#endif
