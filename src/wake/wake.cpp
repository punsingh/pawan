/*! PArticle Wake ANalysis
 * \file wake.cpp
 * \brief Routines for Wake calculations
 *
 * @author Puneet Singh
 * @date 03/28/2021
 */
#include "wake.h"

pawan::__wake::__wake(__io *op){
	_numDimensions = 3;
	_numParticles = 10;
	_position = gsl_matrix_alloc(_numParticles,_numDimensions);
	_vorticity = gsl_matrix_alloc(_numParticles,_numDimensions);
	for(int i = 0; i<_numParticles; ++i){
		for(int j = 0; j<_numDimensions; ++j){
			gsl_matrix_set(_position,i,j,2.0*i+7*j-0.2*i*j);
			gsl_matrix_set(_vorticity,i,j,i-5*j+M_PI*i*j);
		}
	}
	OUT("_position",_position);
	OUT("_vorticity",_vorticity);
	std::string filename = op->getFile() + ".wk";
	FILE *f = fopen(filename.c_str(),"wb");
	gsl_matrix_fwrite(f,_position);
	gsl_matrix_fwrite(f,_vorticity);
	fclose(f);
	FILE *f2 = fopen(filename.c_str(),"rb");
	gsl_matrix_fread(f2,_vorticity);
	gsl_matrix_fread(f2,_position);
	fclose(f2);

}

//pawan::__wake::__wake(const __wake &w){
	//_position = gsl_matrix_alloc(w._position->size1,_numParticles);
	//_vorticity = gsl_matrix_alloc(w._position->size1,_numDimensions);
	//for(int i = 0; i<w._position->size1; ++i){
		//for(int j = 0; j<w._position->size2; ++j){
			//gsl_matrix_set(_position,i,j,gsl_matrix_get(w._position,i,j));
			//gsl_matrix_set(_vorticity,i,j,gsl_matrix_get(w._vorticity,i,j));
		//}
	//}

//}

pawan::__wake::~__wake(){
	gsl_matrix_free(_position);
	gsl_matrix_free(_vorticity);
	gsl_vector_free(_radius);
	gsl_vector_free(_volume);
	gsl_vector_free(_birthstrength);
}

void pawan::__wake::print(){
	OUT("_position",_position);
	OUT("_vorticity",_vorticity);
}

