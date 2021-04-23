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

void pawan::__integration::integrate(__wake *W, __io *IO){
	gsl_vector *states = gsl_vector_calloc(W->_size);
	gsl_vector *rates = gsl_vector_calloc(W->_size);
	FILE *f = IO->create_binary_file(".wake");
	double t = 0.0;
	fwrite(&t,sizeof(double),1,f);	
	W->write(f);
	W->getStates(states);
	for(size_t i = 1; i<=_n; ++i){
		t = i*_dt;
		W->calculateInteraction();
		W->getRates(rates);
		gsl_vector_scale(rates,_dt);
		gsl_vector_add(states,rates);
		W->setStates(states);
		W->print();
		fwrite(&t,sizeof(double),1,f);	
		W->write(f);
	}
	fclose(f);
	
	gsl_vector_free(states);
	gsl_vector_free(rates);
}
