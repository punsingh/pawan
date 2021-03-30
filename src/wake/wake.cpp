/*! PArticle Wake ANalysis
 * \file wake.cpp
 * \brief Routines for Wake calculations
 *
 * @author Puneet Singh
 * @date 03/28/2021
 */
#include "wake.h"

__wake::__wake(){
	_position = gsl_matrix_alloc(10,3);
	_vorticity = gsl_matrix_alloc(10,3);
	for(int i = 0; i<10; ++i){
		for(int j = 0; j<3; ++j){
			gsl_matrix_set(_position,i,j,2.0*i+7*j-0.2*i*j);
			gsl_matrix_set(_vorticity,i,j,i-5*j+M_PI*i*j);
		}
	}

}

__wake::__wake(const __wake &w){
	_position = gsl_matrix_alloc(10,3);
	_vorticity = gsl_matrix_alloc(10,3);
	for(int i = 0; i<w._position->size1; ++i){
		for(int j = 0; j<w._position->size2; ++j){
			gsl_matrix_set(_position,i,j,gsl_matrix_get(w._position,i,j));
			gsl_matrix_set(_vorticity,i,j,gsl_matrix_get(w._vorticity,i,j));
		}
	}

}

__wake::~__wake(){
	gsl_matrix_free(_position);
	gsl_matrix_free(_vorticity);
	gsl_vector_free(_radius);
	gsl_vector_free(_volume);
	gsl_vector_free(_birthstrength);
}

void __wake::print(){
	DISP("_position",_position);
	DISP("_vorticity",_vorticity);
}

