/**
 * @file   math_lib.h
 * @Author A. Fabbri (alessandro.fabbri27@unibo.it), S. Sinigardi (stefano.sinigardi@unibo.it)
 * @date   February, 2016
 * @brief  This header file contains the declaration of structs and functions for the mathematical framework.
 * 
 * This header file contains the declaration of the basic types of objects which implements the 2D/3D/6D vector algebra, the 2D/3D/6D matrix algebra and some utilities for the solution of a 2x2 eigenvalue problem for symmetric matrices (based on analytical solutions).
 *
 */

#ifndef PHYSYCOM_MATH_LIB_H
#define PHYSYCOM_MATH_LIB_H

/*! Symbol to include the built-in display functions. */
#define ENABLE_DISPLAY_FUNCTIONS

#ifdef __cplusplus
extern "C" {
#endif

/*! \struct VEC2D
 *  \brief Two-dimensional vector object
 *
 *  This object holds the data for representing a 2D vector
 */
  typedef struct VEC2D {
    double x;   /*!< The 2D vector \f$x\f$ coordinate */
    double y;   /*!< The 2D vector \f$y\f$ coordinate */
    double mod; /*!< The 2D vector modulus \f$\sqrt{x^2 + y^2}\f$ */
  } VEC2D;

/*! \name 2D Vector algebra functions
 *  Utilities to handle with elementary algebraic 2D operations on vectors
 */
///@{
/*! \function
 *  \brief Function to populate 2D vector.
 *
 *  This routine sets the values of the coordinates and also evaluates the modulus, setting it to zero below a threshold of #EPSILON.
 *
 *  \param v a pointer to the vector to be populated.
 *  \param x a double representing the \f$x\f$ coordinate. 
 *  \param y a double representing the \f$y\f$ coordinate. 
 */
  void set_vec2d(VEC2D * v, const double x, const double y);

/*! \function
 *  \brief Function to automatically set vector modulus from its coordinates.
 *
 *  This routine sets the value of the 2D vector modulus calculating it from the coordinates values, setting it to zero below a threshold of #EPSILON.
 *
 *  \param v a pointer to the vector.
 */
  void setmod_vec2d(VEC2D * v);

/*! \function
 *  \brief Function to normalize the vector.
 *
 *  This routine normalizes the 2D vector, by dividing VEC2D::x and VEC2D::y by VEC2D::mod. Then VEC2D::mod is set to \f$1.0\f$, so a small numerical mismatch between last value and the square root of coordintes' squares sum can occur.
 *
 *  \param v a pointer to the vector.
 */
  void normalize_vec2d(VEC2D * v);

/*! \function
 *  \brief Function to evaluate the scalar product.
 *
 *  This routine evaluates the scalar product between two 2D vectors, by means of \f$\vec{a} \cdot \vec{b} = a_x b_x + a_y b_y\f$.
 *
 *  \param v a pointer to the first vector.
 *  \param u a pointer to the second vector.
 *  \return Scalar product value.
 */
  double prod_dot_2d(const VEC2D * v, const VEC2D * u);

/*! \function
 *  \brief Function to evaluate the cross product.
 *
 *  This routine evaluates the \f$z\f$ coordinate of cross product between two 2D vectors, by means of \f$\vec{a} \times \vec{b} \cdot \hat{k} = a_x b_y - a_y b_x\f$. The positive direction of the \f$z\f$ axis is assumed to be pointing outward the screen.
 *
 *  \param v a pointer to the first vector.
 *  \param u a pointer to the second vector.
 *  \return The \f$z\f$ coordinate of the cross product.
 */
  double prod_cross_2d(const VEC2D * v, const VEC2D * u);

/*! \function
 *  \brief Function to multiply a vector by a scalar.
 *
 *  This routine implements the multiplication of a 2D vector by a scalar, the original vector coordinates are overwritten in the process.
 *
 *  \param v a pointer to the vector to be multiplied.
 *  \param alpha the scalar multiplier.
 */
  void multiply_vec2d(VEC2D * v, const double alpha);
///@}

/*! \struct MAT2D
 *  \brief Two-dimensional matrix object
 *
 *  This object holds the data for representing a 2D square matrix
 */
  typedef struct MAT2D {
    double xx;  /*! The element \f$a_{00}\f$ of the matrix. */
    double xy;  /*! The element \f$a_{01}\f$ of the matrix. */
    double yx;  /*! The element \f$a_{10}\f$ of the matrix. */
    double yy;  /*! The element \f$a_{11}\f$ of the matrix. */
  } MAT2D;

/*! \name 2D Matrix algebra functions
 *  Utilities to handle with elementary algebraic 2D operations on square matrices
 */
///@{
/*! \function set_mat2d
 *  \brief Function to populate 2D matrix.
 *
 *  This routine sets the value of the matrix elements.
 *
 *  \param m a pointer to the matrix to be populated.
 *  \param xx a double representing the element \f$a_{00}\f$ of the matrix. 
 *  \param xy a double representing the element \f$a_{01}\f$ of the matrix. 
 *  \param yx a double representing the element \f$a_{10}\f$ of the matrix. 
 *  \param yy a double representing the element \f$a_{11}\f$ of the matrix. 
 */
  void set_mat2d(MAT2D * m, double xx, double xy, double yx, double yy);

/*! \function
 *  \brief Function to transpose a 2D matrix.
 *
 *  This routine evaluates the transpose of a given matrix, i.e. a matrix \f$B\f$ whose elements are given by \f$ b_{ij} = a_{ji}\f$ in terms of the elements of the original matrix \f$A\f$.
 *
 *  \param mt a pointer to a matrix to store the transposition result.
 *  \param m a pointer to the original matrix.
 */
  void transpose_mat2d(MAT2D * mt, const MAT2D * m);

/*! \function
 *  \brief Function to multiply two 2D matrix.
 *
 *  This routine evaluates the product of two given matrices \f$A\f$ and \f$B\f$ according to the standard rule \f$c_{ij} = a_{ik} b_{kj}\f$ (Einstein summation convention assumed).
 *
 *  \param result a pointer to a matrix to store the product result.
 *  \param a a pointer to the left matrix.
 *  \param a a pointer to the right matrix.
 */
  void product_mat2d(MAT2D * result, const MAT2D * a, const MAT2D * b);

/*! \function
 *  \brief Function to multiply a 2D matrix by a scalar.
 *
 *  This routine implements the multiplication of a given matrix by a given scalar, the original matrix elements are overwritten in the process.
 *
 *  \param alpha a double representing the scalar factor.
 *  \param m a pointer to matrix to be multiplied.
 */
  void multiply_mat2d(double alpha, MAT2D * m);

/*! \function
 *  \brief Function to generate a 2D rotation.
 *
 *  This routine populates a matrix with the values corresponding to a counter-clockwise rotation around the origin of an angle whose sine and cosine values are provided, namely \f$\left( \begin{array}{cc} \cos \theta & -\sin \theta \\ \sin \theta & \cos \theta \end{array} \right)\f$. It is implemented as an inline function.
 *
 *  \param rot a pointer to the rotation matrix to be populated.
 *  \param c the cosine of the rotation angle \f$\cos \theta\f$.
 *  \param s the sine of the rotation angle \f$\sin \theta\f$.
 */
  void make_rotation_2d_cs(MAT2D * rot, const double c, const double s);

/*! \function
 *  \brief Function to generate a 2D rotation.
 *
 *  This routine populates a matrix with the values corresponding to a counter-clockwise rotation around the origin of an angle in radians, namely \f$\left( \begin{array}{cc} \cos \theta & -\sin \theta \\ \sin \theta & \cos \theta \end{array} \right)\f$.
 *
 *  \param rot a pointer to the rotation matrix to be populated.
 *  \param angle the angle of rotation in radians.
 */
  void make_rotation_2d(MAT2D * rot, const double angle);

/*! \function
 *  \brief Function to apply a rotation to a 2D vector.
 *
 *  This routine rotates a given vector with a given rotation matrix, it actually implements the product of a 2D matrix with a 2D vector, namely \f$u_i = a_{ij}v_j\f$ (Einstein summation convention assumed).
 *
 *  \param result a pointer to 2D vector to store the result.
 *  \param rot a pointer to the (rotation) matrix to be applied.
 *  \param v a pointer to the original vector.
 */
  void rotate_vec2d(VEC2D * result, const MAT2D * rot, const VEC2D * v);

/*! \function
 *  \brief Function to apply a rotation to a 2D matrix.
 *
 *  This routine rotates a given matrix with a given rotation matrix, it actually implements the matrix transformation \f$A \rightarrow R^T A R\f$, orthogonal matrix version of a matrix similarity transformation.
 *
 *  \param result a pointer to 2D matrix to store the result.
 *  \param rot a pointer to the (rotation) matrix to be applied.
 *  \param m a pointer to the original matrix.
 */
  void rotate_mat2d(MAT2D * result, const MAT2D * rot, const MAT2D * m);
///@}

/*! \struct VEC3D
 *  \brief Three-dimensional vector object.
 *
 *  This object holds the data for representing a 3D vector.
 */
  typedef struct VEC3D {
    double x, y, z, mod;
  } VEC3D;

/*! \name 3D Matrix algebra functions
 *  Utilities to handle with elementary algebraic 3D operations on vectors.
 */
///@{
/*! \function 
 *  \brief Function to populate 3D vector.
 *  \param v a pointer to the vector to be populated.
 *  \param x a double representing the \f$x\f$ coordinate. 
 *  \param y a double representing the \f$y\f$ coordinate. 
 *  \param z a double representing the \f$z\f$ coordinate. 
 */
  void set_vec3d(VEC3D * v, const double x, const double y, const double z);

  void setmod_vec3d(VEC3D * v);

  void normalize_vec3d(VEC3D * v);

  double prod_dot_3d(const VEC3D * a, const VEC3D * b);

  void prod_cross_3d(VEC3D * result, const VEC3D * a, const VEC3D * b);

  void multiply_vec3d(VEC3D * v, const double alpha);
///@}

/*! \struct MAT3D
 *  \brief Three-dimensional square matrix object.
 *
 *  This object holds the data for representing a 3D square matrix.
 */
  typedef struct MAT3D {
    double xx, xy, xz, yx, yy, yz, zx, zy, zz;
  } MAT3D;

/*! \name 3D Matrix algebra functions
 *  Utilities to handle with elementary algebraic 3D operations on square matrices.
 */
///@{
/*! \function 
 *  \brief Function to populate 3D vector.
 *  \param v a pointer to the vector to be populated.
 *  \param xx a double representing the element \f$a_{00}\f$ of the matrix. 
 *  \param xy a double representing the element \f$a_{01}\f$ of the matrix. 
 *  \param yz a double representing the element \f$a_{02}\f$ of the matrix. 
 *  \param yx a double representing the element \f$a_{10}\f$ of the matrix. 
 *  \param yy a double representing the element \f$a_{11}\f$ of the matrix. 
 *  \param yz a double representing the element \f$a_{12}\f$ of the matrix. 
 *  \param zx a double representing the element \f$a_{20}\f$ of the matrix. 
 *  \param zy a double representing the element \f$a_{21}\f$ of the matrix. 
 *  \param zz a double representing the element \f$a_{22}\f$ of the matrix. 
 */
  void set_mat3d(MAT3D * m, const double xx, const double xy, const double xz, const double yx, const double yy, const double yz, const double zx, const double zy, const double zz);

  void transpose_mat3d(MAT3D * result, const MAT3D * m);

  void product_mat3d(MAT3D * result, const MAT3D * a, const MAT3D * b);

  void multiply_mat3d(MAT3D * m, const double alpha);

  void make_rotation_cs(MAT3D * rot, const VEC3D * v, const double c, const double s);

  void make_rotation(MAT3D * rot, const VEC3D * v, const double angle);

  void rotate_vec3d(VEC3D * result, const MAT3D * rot, const VEC3D * v);

  void rotate_mat3d(MAT3D * result, const MAT3D * rot, const MAT3D * m);
///@}

/*! \struct VEC6D
 *  \brief Six-dimensional vector object.
 *
 *  This object holds the data for representing a 6D vector. More precisely, this is not an element of the vector space \f$\mathbb{R}^6\f$ rather it represent a collection of two 3D vectors and provides utilities function to manipulate them at once.
 */
  typedef struct VEC6D {
    VEC3D a;  /*!< The first 6D vector sub-object */ 
    VEC3D g;  /*!< The second 6D vector sub-object */ 
  } VEC6D;

/*! \name 6D Vector algebra functions
 *  Utilities to handle with elementary algebraic 6D operations on vectors.
 */
///@{
/*! \function 
 *  \brief Function to populate 6D vector.
 *  \param v a pointer to the vector to be populated.
 *  \param ax a double representing the \f$x\f$ coordinate of first vector. 
 *  \param ay a double representing the \f$y\f$ coordinate of first vector. 
 *  \param az a double representing the \f$z\f$ coordinate of first vector. 
 *  \param gx a double representing the \f$x\f$ coordinate of second vector. 
 *  \param gy a double representing the \f$y\f$ coordinate of second vector. 
 *  \param gz a double representing the \f$z\f$ coordinate of second vector. 
 */
  void set_vec6d(VEC6D * v, const double ax, const double ay, const double az, const double gx, const double gy, const double gz);

  void rotate_vec6d(VEC6D * result, const MAT3D * rot, const VEC6D * v);
///@}

/*! \struct MAT6D
 *  \brief Six-dimensional square matrix object.
 *
 *  This object holds the data for representing a 6D square matrix. 
 */
  typedef struct MAT6D {
    MAT3D A[2][2];
  } MAT6D;

/*! \name 6D Square matrix algebra functions
 *  Utilities to handle with elementary algebraic 6D operations on square matrices.
 */
///@{
/*! \function 
 *  \brief Function to populate 6D matrix.
 *  \param m a pointer to the matrix to be populated.
 *  \param A00 a pointer to 3D matrix representing the \f$A_{00}\f$ sub-matrix. 
 *  \param A01 a pointer to 3D matrix representing the \f$A_{01}\f$ sub-matrix.
 *  \param A10 a pointer to 3D matrix representing the \f$A_{10}\f$ sub-matrix.
 *  \param A11 a pointer to 3D matrix representing the \f$A_{11}\f$ sub-matrix.
 */
  void set_mat6d(MAT6D * m, const MAT3D * A00, const MAT3D * A01, const MAT3D * A10, const MAT3D * A11);

  void rotate_mat6d(MAT6D * result, const MAT3D * rot, const MAT6D * m);
///@}

/*! \struct EigenSys
 *  \brief Object containing the data for a 2D eigenvalue problem.
 *
 *  This object holds the data related to a two-dimensional real-valued eigenvalue problem for symmetric matrices. 
 */
  typedef struct EigenSys {
    MAT2D A;    /*!< The 2D symmetric matrix to be diagonalized. */  
    double l1;  /*!< Largest eigenvalue. */ 
    double l2;  /*!< Smallest eigenvalue. */ 
    VEC2D v1;   /*!< Eigenvector corresponding to the largest eigenvalue. */
    VEC2D v2;   /*!< Eigenvector corresponding to the smallest eigenvalue. */
    VEC2D u1;   /*!< Normalized eigenvector corresponding to the largest eigenvalue. */
    VEC2D u2;   /*!< Normalized eigenvector corresponding to the smallest eigenvalue. */
  } EigenSys;

/*! \name 2D Eigenproblems utilities
 *  Tools to set up and solve the eigenvalue problem for input 2D symmetric real-valued matrices.
 */
///@{
/*! \function 
 *  \brief Return the eigensystem of a 2D matrix.
 *
 *  Routine that solves, by means of exact analytical formulas, the eigenvalue for a matrix of the form \f$\left( \begin{array}{cc} 1 & a \\ a & b\end{array} \right)\f$, for real \f$a,b\f$.
 *
 *  \param a element \f$a_{01} = a_{10}\f$ of the matrix.
 *  \param b element \f$a_{11}\f$ of the matrix.
 */
  EigenSys eigs_2x2_sym_normalized(double a, double b);

  EigenSys eigs_2x2_sym(MAT2D m);

  int check_eigs(EigenSys es);
///@}

#ifdef ENABLE_DISPLAY_FUNCTIONS
/*! \name Display functions
 *  Support function to display the value of the various mathematical objects.
 */
///@{
/*! \function 
 *  \brief Display 2D vector.
 *
 *  This routine displays to stdout the coordinates of a given vector along with a given name for it.
 *
 *  \param v pointer to the vector to be displayed.
 *  \param name a char pointer containing the name of the vector.
 */
  void print_vec2d(const VEC2D * v, const char * name);
  
  void print_mat2d(const MAT2D * m, const char * name);

  void print_vec3d(const VEC3D * v, const char * name);

  void print_mat3d(const MAT3D * m, const char * name);

  void print_vec6d(const VEC6D * v, const char * name);

  void print_mat6d(const MAT6D * m, const char * name);

  void print_eigs(EigenSys, const char *);
///@}
#endif

#ifdef __cplusplus
}
#endif

#endif
