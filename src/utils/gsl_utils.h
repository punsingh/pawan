/*! PArticle Wake ANalysis
 * @file gsl_utils.h
 * @brief GNU Scientific Library utilities for PAWAN
 *
 * @author Puneet Singh
 * @date 03/28/2021
 *
 */

#ifndef GSL_UTILS_H_
#define GSL_UTILS_H_

#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>


/*
 *
 * PRINT OPERATIONS
 *
 */

/*! \fn DISP
 * \brief Print string and complex number
 * \param	s	String
 * \param	c	Complex number
 */
inline void DISP(std::string s, const gsl_complex &c, std::ostream &os = std::cout){
	os << "\t" << s << " = "<< GSL_REAL(c) << " +i " << GSL_IMAG(c) << std::endl;
};

/*! \fn DISP
 * \brief Print string and long array of values
 * \param	s	String
 * \param	v	gsl vector
 */
inline void DISP(std::string s, const gsl_vector *v, std::ostream &os = std::cout){
	if(v->size==0){
		os << "\t" << s << " is empty."<< std::endl;
	}
	else{
		os << "\t" << s << " = "<< std::endl;
		for(int i = 0; i<v->size; ++i){
		       	os << "\t" << gsl_vector_get(v,i) << std::endl;
       		}
	}
};

/*! \fn DISP
 * \brief Print string and a limited array of values
 * \param	s	String
 * \param	v	gsl vector
 * \param	n	int
 */
inline void DISP(std::string s, const gsl_vector *v, const int &n, std::ostream &os = std::cout){
	int m = n<=v->size?n:v->size;
	if(m==0){
		os << "\t" << s << " is empty."<< std::endl;
	}
	else{
		os << "\t" << s << " = "<< std::endl;
		for(int i = 0; i<m; ++i){
		       	os << "\t" << gsl_vector_get(v,i) << std::endl;
       		}
	}
};

/*! \fn DISPT
 * \brief Print string and long array of values Transposed
 * \param	s	String
 * \param	v	gsl vector
 */
inline void DISPT(std::string s, const gsl_vector *v, std::ostream &os = std::cout){
	if(v->size==0){
		os << "\t" << s << " is empty."<< std::endl;
	}
	else{
		os << "\t" << s << " =";
		for(int i = 0; i<v->size; i++){
		       	os << "\t" << gsl_vector_get(v,i);
       		}
		os << std::endl;
	}
};

/*! \fn DISPT
 * \brief Print string and limited array of values transposed
 * \param	s	String
 * \param	v	gsl vector
 * \param	n	int
 */
inline void DISPT(std::string s, const gsl_vector *v, const int &n, std::ostream &os = std::cout){
	int m = n<=v->size?n:v->size;
	if(m==0){
		os << "\t" << s << " is empty."<< std::endl;
	}
	else{
		os << "\t" << s << " =";
		for(int i = 0; i<m; ++i){
		       	os << "\t" << gsl_vector_get(v,i);
       		}
		os << std::endl;
	}
};

/*! \fn DISP
 * \brief Print string and matrix of values
 * \param	s	String
 * \param	m	gsl matrix
 */
inline void DISP(std::string s, const gsl_matrix *m, std::ostream &os = std::cout){
	if(m->size1==0 && m->size2==0){
		os << "\t" << s << " is empty."<< std::endl;
	}
	else{
		os << "\t" << s << " = "<< std::endl;
		for(int i = 0; i<m->size1; ++i){
			os << "\t";
			for(int j = 0; j<m->size2; ++j){
		       		os << "\t" << gsl_matrix_get(m,i,j);
			}
			os << std::endl;
       		}
	}
};

/*
 *
 * VECTOR OPERATIONS
 *
 */

/*! \fn get_gsl_vector
 * \brief	Write gsl vector to double array
 * \param	x	gsl vector
 * \param	A	double array
 * \param	n	Length of vector
 */
inline void get_gsl_vector(const gsl_vector *x, double *A, const int &n){
	for(int i = 0; i<n; ++i){
		A[i] = gsl_vector_get(x,i);
	}
};

/*! \fn get_gsl_vector
 * \brief	Write gsl vector to double array
 * \param	x	gsl vector
 * \param	A	double array
 * \param	n	Length of vector
 * \param	offset	Vector offset
 */
inline void get_gsl_vector(const gsl_vector *x, double *A, const int &n, const int &offset){
	for(int i = 0; i<n; ++i){
		A[i] = gsl_vector_get(x,i+offset);
	}
};

/*! \fn set_gsl_vector
 * \brief	Write double array to gsl vector
 * \param	A	double array
 * \param	x	gsl vector
 * \param	n	Length of vector
 */
inline void set_gsl_vector(const double *A, gsl_vector *x, const int &n){
	for(int i = 0; i<n; ++i){
		gsl_vector_set(x,i,A[i]);
	}
};

/*! \fn set_gsl_vector
 * \brief	Write gsl vector to double array
 * \param	A	double array
 * \param	x	gsl vector
 * \param	n	Length of vector
 * \param	offset	Vector offset
 */
inline void set_gsl_vector(const double *A, gsl_vector *x, const int &n, const int &offset){
	for(int i = 0; i<n; ++i){
		gsl_vector_set(x,i+offset,A[i]);
	}
};

/*! \fn zero
 * \brief Set gsl vector to zeros
 * \param V gsl vector 
 * \param n length
 */
inline void zero(gsl_vector *V){
	gsl_vector_set_zero(V);
};

/*! \fn increment_gsl_vector
 * \brief	Increment gsl vector element
 * \param	x	gsl vector
 * \param	n	index
 * \param 	e	incremental value
 */
inline void increment_gsl_vector(gsl_vector *x, const int &n, const double &e){
	gsl_vector_set(x,n,gsl_vector_get(x,n) + e);
};

/*
 * \fn cartesian2spherical
 * \brief Cartesian to Spherical coordinate conversion
 * \param cartesian	Coordinate
 * \param r	Radius
 * \param theta	Azimuth
 * \param phi	Latitude 
 *
 */
inline void cartesian2spherical(gsl_vector *cartesian, double &r, double &theta, double &phi){
	double x = gsl_vector_get(cartesian,0);
	double y = gsl_vector_get(cartesian,1);
	double z = gsl_vector_get(cartesian,2);
	r = gsl_hypot3(x,y,z);
	if(r==0.0){
		phi = 0.0;
		theta = 0.0;
	}
	else{
		phi = atan2(y,x);
		theta = atan2(gsl_hypot(x,y),z);
	}
}

/*! \fn norm
 * \brief Distance
 * \param A vector 
 * \param B vector
 * Returns sqrt((A-B)*(A-B)')
 */
inline double norm(const gsl_vector *A, const gsl_vector *B){
	double result = 0.0;
	for(int i = 0; i<A->size; ++i){
		result += gsl_pow_2(gsl_vector_get(A,i) - gsl_vector_get(B,i));
	}
	return sqrt(result);
}

/*! \fn cross
 * \brief Vector cross product
 * \param A vector 
 * \param B vector
 * \param C vector
 * Returns A = B x C
 */
inline void cross(gsl_vector *A, const gsl_vector *B, const gsl_vector *C){
	double Bx = gsl_vector_get(B,0);
	double By = gsl_vector_get(B,1);
	double Bz = gsl_vector_get(B,2);
	double Cx = gsl_vector_get(C,0);
	double Cy = gsl_vector_get(C,1);
	double Cz = gsl_vector_get(C,2);
	gsl_vector_set(A,0,By*Cz - Bz*Cy);
	gsl_vector_set(A,1,Bz*Cx - Bx*Cz);
	gsl_vector_set(A,2,Bx*Cy - By*Cx);
}

/*! \fn flip_sign
 * \brief Flip sign of vector
 * \param A vector
 */
inline void flip_sign(gsl_vector *A){
	for(int i = 0; i<A->size; ++i) {
		gsl_vector_set(A,i,-gsl_vector_get(A,i));
	}
}

/*
 *
 * MATRIX OPERATIONS
 *
 */

/*! \fn get_gsl_matrix
 * \brief	Write gsl vector to double array
 * \param	M	gsl matrix
 * \param	A	double matrix
 * \param	n	number of rows
 * \param	m	number of cols
 */
inline void get_gsl_matrix(const gsl_matrix *M, double *A, const int &n, const int &m){
	for(int i = 0; i<n; ++i){
		for(int j = 0; j<m; ++j){
			A[i*m+j] = gsl_matrix_get(M,i,j);
		}
	}
};

/*! \fn set_gsl_matrix
 * \brief	Write double matrix to gsl matrix
 * \param	A	double matrix
 * \param	M	gsl matrix
 * \param	n	number of rows
 * \param	m	number of cols
 */
inline void set_gsl_matrix(const double *A, gsl_matrix *M, const int &n, const int &m){
	for(int i = 0; i<n; ++i){
		for(int j = 0; j<m; ++j){
			gsl_matrix_set(M,i,j,A[i*m+j]);
		}
	}
};

#endif
