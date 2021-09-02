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
    std::cout << "integration.cpp---------------------------------------------------real execution begins here------------------------------------------------------------ entry into __integration::integrate" << std::endl;
	gsl_vector *states = gsl_vector_calloc(S->_size);   /* allocates memory and initializes all elements to 0*/
	FILE *f = IO->create_binary_file(".wake");
	double t = 0.0;
	fwrite(&t,sizeof(double),1,f);	//writing out initial time
	S->write(f);
	S->getStates(states);
	double tStart = TIME();
	for(size_t i = 1; i<=_n; ++i){
		OUT("\tStep",i);
		t = i*_dt;
		step(_dt,S,states);
		fwrite(&t,sizeof(double),1,f);  //writing out each time step
//        std::cout << "integration.cpp------writing to .wake file" << std::endl;
		S->write(f);
	}
	fclose(f);
	double tEnd = TIME();
	OUT("Total Time (s)",tEnd - tStart);
	gsl_vector_free(states);
    std::cout << "integration.cpp------------------------------------------------------------------all execution ends------------------------------------------------------ exit from __integration::integrate" << std::endl;
}

void pawan::__integration::step(const double &dt, __interaction *S, gsl_vector* states){
//    std::cout << "---------------------------step() in integration.cpp being executed" << std::endl;
	gsl_vector *rates = gsl_vector_calloc(S->_size);
	S->interact();
	S->getRates(rates);
	gsl_vector_scale(rates,dt);
	gsl_vector_add(states,rates);
	S->setStates(states);
	gsl_vector_free(rates);
}
