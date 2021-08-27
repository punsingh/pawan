/*! PArticle Wake ANalysis
 * \file vring.cpp
 * \brief Routines for creating Vortex rings
 *
 * @author Puneet Singh
 * @date 08/12/2021
 */
#include "vring.h"

pawan::__vring::__vring(const double &radius, const double &core, const int &nRadial, const int &nLayers, const bool &shift){
	create_particles(nRadial*gsl_pow_2(2*nLayers-1));
	double dPsi = 2.0*M_PI/nRadial; 					/*!< Azimuthal step size */ 
	double K = -log(0.01);
	double core2 = gsl_pow_2(core);
	double r0 = core/(2*nLayers - 1);					/*!< Central area */ 
	double sigma = pow(2*gsl_pow_2(M_PI*core)*radius/_numParticles,1./3.);	/*!< Smoothing radius */ 
	double volume = 4.0*M_PI*gsl_pow_3(sigma)/3.0;				/*!< Volume of annular segment */ 
	size_t index = 0;	
	double delpsi = shift?0.5*dPsi:0.; 
	for(size_t i = 0; i<nRadial; ++i){
		double psi = i*dPsi + delpsi;
		double sPsi = sin(psi);
		double cPsi = cos(psi);

		gsl_matrix_set(_position,index,0,radius*cPsi);
		gsl_matrix_set(_position,index,1,radius*sPsi);
		gsl_matrix_set(_position,index,2,0.0);
		gsl_matrix_set(_vorticity,index,0,-sPsi);
		gsl_matrix_set(_vorticity,index,1,cPsi);
		gsl_matrix_set(_vorticity,index,2,0.0);
		gsl_vector_set(_radius,index,sigma);
		gsl_vector_set(_volume,index,volume);
		gsl_vector_set(_birthstrength,index,1.0);
		index++;
		for(size_t j = 1; j<nLayers; ++j){
			double r = r0*(12.*gsl_pow_2(j)+1.0)/(6.*j);
			double dr = r0*(2. - 1./(6.*j*(j+1.)));
			size_t nCore = 8*j;
			double dTheta = 2.*M_PI/nCore;
			for(size_t k = 0; k<nCore; ++k){
				double theta = k*dTheta;
				double x = (radius + r*cos(theta))*cPsi;
				double y = (radius + r*cos(theta))*sPsi;
				double z = r*sin(theta);
				double ax = -sPsi;
				double ay = cPsi;
				gsl_matrix_set(_position,index,0,x);
				gsl_matrix_set(_position,index,1,y);
				gsl_matrix_set(_position,index,2,z);
				gsl_matrix_set(_vorticity,index,0,ax);
				gsl_matrix_set(_vorticity,index,1,ay);
				gsl_matrix_set(_vorticity,index,2,0.0);
				gsl_vector_set(_radius,index,sigma);
				gsl_vector_set(_volume,index,volume);
				gsl_vector_set(_birthstrength,index,1.0);
				index++;
			}
		}
	}
}


pawan::__vring::~__vring(){
}
