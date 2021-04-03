/*! PArticle Wake ANalysis
 * \file wake.cpp
 * \brief Routines for Wake calculations
 *
 * @author Puneet Singh
 * @date 03/28/2021
 */
#include "wake.h"

void pawan::__wake::create_particles(const int &n){
	_numDimensions = 3;
	_numParticles = n;
	_position = gsl_matrix_alloc(_numParticles,_numDimensions);
	_vorticity = gsl_matrix_alloc(_numParticles,_numDimensions);
	_radius = gsl_vector_alloc(_numParticles);
	_volume = gsl_vector_alloc(_numParticles);
	_birthstrength = gsl_vector_alloc(_numParticles);
}

pawan::__wake::__wake(){
	create_particles(1);
	for(int i = 0; i<_numParticles; ++i){
		for(int j = 0; j<_numDimensions; ++j){
			gsl_matrix_set(_position,i,j,0.0);
			gsl_matrix_set(_vorticity,i,j,0.0);
		}
		gsl_vector_set(_radius,i,1.0);
		gsl_vector_set(_volume,i,1.0);
		gsl_vector_set(_birthstrength,i,1.0);
	}
}

pawan::__wake::__wake(const double &gamma, const double &radius, const double &core, const int &nRadial){
	create_particles(nRadial);
	double dPsi = 2.0*M_PI/nRadial; 		/*!< Azimuthal step size */ 
	double psi = 0.0; 				/*!< Initial azimuth location */ 
	double strength = gamma*radius*dPsi;		/*!< Vorticity magnitude */ 
	double volume = 4.0*M_PI*gsl_pow_3(core)/3.0;	/*!< Volume of annular segment */ 
	for(int i = 0; i<_numParticles; ++i){
		psi = i*dPsi;
		gsl_matrix_set(_position,i,0,radius*cos(psi));
		gsl_matrix_set(_position,i,1,radius*sin(psi));
		gsl_matrix_set(_position,i,2,0.0);
		gsl_matrix_set(_vorticity,i,0,-strength*sin(psi));
		gsl_matrix_set(_vorticity,i,1,strength*cos(psi));
		gsl_matrix_set(_vorticity,i,2,0.0);
		gsl_vector_set(_radius,i,core);
		gsl_vector_set(_volume,i,volume);
		gsl_vector_set(_birthstrength,i,strength);
	}
}

pawan::__wake::~__wake(){
	gsl_matrix_free(_position);
	gsl_matrix_free(_vorticity);
	gsl_vector_free(_radius);
	gsl_vector_free(_volume);
	gsl_vector_free(_birthstrength);
}

void pawan::__wake::print(){
	OUT("_numParticles",_numParticles);
	OUT("_position",_position);
	OUT("_vorticity",_vorticity);
	OUT("_radius",_radius);
	OUT("_volume",_volume);
	OUT("_birthstrength",_birthstrength);
}

void pawan::__wake::write(__io *IO){
	FILE *f = IO->create_binary_file(".wake");
	fwrite(&_numParticles,sizeof(size_t),1,f);	
	gsl_matrix_fwrite(f,_position);
	gsl_matrix_fwrite(f,_vorticity);
	gsl_vector_fwrite(f,_radius);
	gsl_vector_fwrite(f,_volume);
	gsl_vector_fwrite(f,_birthstrength);
	fclose(f);
}

void pawan::__wake::read(__io *IO){
	FILE *f = IO->open_binary_file(".wake");
	fread(&_numParticles,sizeof(size_t),1,f);	
	gsl_matrix_fread(f,_position);
	gsl_matrix_fread(f,_vorticity);
	gsl_vector_fread(f,_radius);
	gsl_vector_fread(f,_volume);
	gsl_vector_fread(f,_birthstrength);
	fclose(f);
}

