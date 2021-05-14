/*! PArticle Wake ANalysis
 * \file square.cpp
 * \brief Routines for creating Vortex squares
 *
 * @author Puneet Singh
 * @date 05/13/2021
 */
#include "square.h"

pawan::__square::__square(const double &gamma, const double &side, const double &core, const int &nSide){
	create_particles(4*nSide);
	double volume = 4.0*M_PI*gsl_pow_3(core)/3.0;	/*!< Volume of particle */ 
	double length = side/nSide;			/*!< Length of particle */ 
	double strength = gamma*length;			/*!< Strength of particle */ 
	double dist[4] = {0.5*side,0.5*side,-0.5*side,-0.5*side};
	double vort[4] = {strength,-strength,-strength,strength};
	for(size_t i = 0; i<4; ++i){
		size_t axisPerp = i%2;
		size_t axisParl = 1 - axisPerp;
		double sign = pow(-1,(i/2)%2); 
		for(size_t j = 0; j<nSide; ++j){
			size_t index = i*nSide + j;
			gsl_matrix_set(_position,index,axisPerp,dist[i]);
			gsl_matrix_set(_position,index,axisParl,sign*(-0.5*side + length*(j + 0.5)));
			gsl_matrix_set(_position,index,2,0.0);
			gsl_matrix_set(_vorticity,index,axisPerp,0.0);
			gsl_matrix_set(_vorticity,index,axisParl,vort[i]);
			gsl_matrix_set(_vorticity,index,2,0.0);
			gsl_vector_set(_radius,index,core);
			gsl_vector_set(_volume,index,volume);
			gsl_vector_set(_birthstrength,index,strength);
		}
	}
}

pawan::__square::__square(const double &gamma, const double &radius, const double &core, const int &nRadial, const int &nLayers, const bool &shift){
	//create_particles(nRadial*gsl_pow_2(2*nLayers-1));
	//double dPsi = 2.0*M_PI/nRadial; 					[>!< Azimuthal step size <] 
	//double K = -log(0.01);
	//double core2 = gsl_pow_2(core);
	//double strength = gamma*core2/K;					[>!< Vorticity magnitude <] 
	//double r0 = core/(2*nLayers - 1);					[>!< Central area <] 
	//double sigma = pow(2*gsl_pow_2(M_PI*core)*radius/_numParticles,1./3.);	[>!< Smoothing radius <] 
	//double volume = 4.0*M_PI*gsl_pow_3(sigma)/3.0;				[>!< Volume of annular segment <] 
	//size_t index = 0;	
	//double mag0 = strength*(1.0-exp(-K*gsl_pow_2(r0/core)))/(r0*r0);
	//double delpsi = shift?0.5*dPsi:0.; 
	//for(size_t i = 0; i<nRadial; ++i){
		//double psi = i*dPsi + delpsi;
		//double sPsi = sin(psi);
		//double cPsi = cos(psi);

		//gsl_matrix_set(_position,index,0,radius*cPsi);
		//gsl_matrix_set(_position,index,1,radius*sPsi);
		//gsl_matrix_set(_position,index,2,0.0);
		//gsl_matrix_set(_vorticity,index,0,-mag0*sPsi);
		//gsl_matrix_set(_vorticity,index,1,mag0*cPsi);
		//gsl_matrix_set(_vorticity,index,2,0.0);
		//gsl_vector_set(_radius,index,sigma);
		//gsl_vector_set(_volume,index,volume);
		//gsl_vector_set(_birthstrength,index,mag0);
		//index++;
		//for(size_t j = 1; j<nLayers; ++j){
			//double r = r0*(12.*gsl_pow_2(j)+1.0)/(6.*j);
			//double dr = r0*(2. - 1./(6.*j*(j+1.)));
			//size_t nCore = 8*j;
			//double dTheta = 2.*M_PI/nCore;
			//for(size_t k = 0; k<nCore; ++k){
				//double theta = k*dTheta;
				//double x = (radius + r*cos(theta))*cPsi;
				//double y = (radius + r*cos(theta))*sPsi;
				//double z = r*sin(theta);
				//double mag = strength*exp(-K*gsl_pow_2(r/core))*sinh(K*r*dr/core2)/core/dr;
				//double ax = -mag*sPsi;
				//double ay = mag*cPsi;
				//gsl_matrix_set(_position,index,0,x);
				//gsl_matrix_set(_position,index,1,y);
				//gsl_matrix_set(_position,index,2,z);
				//gsl_matrix_set(_vorticity,index,0,ax);
				//gsl_matrix_set(_vorticity,index,1,ay);
				//gsl_matrix_set(_vorticity,index,2,0.0);
				//gsl_vector_set(_radius,index,sigma);
				//gsl_vector_set(_volume,index,volume);
				//gsl_vector_set(_birthstrength,index,mag);
				//index++;
			//}
		//}
	//}
}


pawan::__square::~__square(){
}
