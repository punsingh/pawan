/*! PArticle Wake ANalysis
 * \file integration.cpp
 * \brief Routines for Integrations
 *
 * @author Puneet Singh
 * @date 04/21/2021
 */
#include "integration.h"

pawan::__integration::__integration(const double &t, const size_t &n){
	_dt = t/n;
	_t = t;
	_n = n;
}

void pawan::__integration::integrate(__interaction *S, __io *IO){
	gsl_vector *states = gsl_vector_calloc(S->_size);
	gsl_vector *rates = gsl_vector_calloc(S->_size);
	FILE *f = IO->create_binary_file(".wake");
	double t = 0.0;
	fwrite(&t,sizeof(double),1,f);	
	S->write(f);
	S->getStates(states);
	for(size_t i = 1; i<=_n; ++i){
		t = i*_dt;
		S->interact();
		S->getRates(rates);
		gsl_vector_scale(rates,_dt);
		gsl_vector_add(states,rates);
		S->setStates(states);
		fwrite(&t,sizeof(double),1,f);	
		S->write(f);
	}
	fclose(f);
	
	gsl_vector_free(states);
	gsl_vector_free(rates);
}
