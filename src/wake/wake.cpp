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

__wake::~__wake(){
	gsl_matrix_free(_position);
	gsl_matrix_free(_vorticity);
}

void __wake::print(){
	for(int i = 0; i<10; ++i){
		for(int j = 0; j<3; ++j){
			std::cout << "\t" << gsl_matrix_get(_position,i,j); 
		}
		std::cout << std::endl;
	}
	for(int i = 0; i<10; ++i){
		for(int j = 0; j<3; ++j){
			std::cout << "\t" << gsl_matrix_get(_vorticity,i,j); 
		}
		std::cout << std::endl;
	}

}

