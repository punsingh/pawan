/*! PArticle Wake ANalysis
 * \file integration.cpp
 * \brief Routines for Integrations
 *
 * @author Puneet Singh
 * @date 04/21/2021
 */
#include "integration.h"
#include <unistd.h>

pawan::__integration::__integration(const double &t, const size_t &n){
    DOUT("-----------------------in pawan::__integration::__integration()");
	_dt = t/n;
	_t = t;
	_n = n;
}

void pawan::__integration::integrate(__system *S, __io *IO, NetworkInterfaceTCP<OPawanRecvData,OPawanSendData> *networkCommunicatorTest){
    gsl_vector *states = gsl_vector_calloc(S->_maxsize); //ip: have the *rates pointer here outside too for speed up
    FILE *f = IO->create_binary_file(".wake");
    double t = 0.0;
    fwrite(&t,sizeof(double),1,f);
    S->write(f);  //write particles info as is

    double tStart = TIME();
    for(size_t i = 1; i<=_n; ++i){
        S->getStates(states);
        OUT("\tStep",i);
        t = i*_dt;
        step(_dt,S,states);
        fwrite(&t,sizeof(double),1,f);
        S->write(f);  //write particles info after interaction in this time step
        //S->updateVinfEffect(_dt,states);
        S->updateVinfEffect(_dt);
        OPawanSendData opawansenddata;
        networkCommunicatorTest->send_data(opawansenddata);
        //sleep(5);
        if (i!=_n) {
            OPawanRecvData opawanrecvdata;
            networkCommunicatorTest->recieve_data(opawanrecvdata);
            S->addParticles(&opawanrecvdata);
        }


    }
    fclose(f);
    double tEnd = TIME();
    OUT("Total Time (s)",tEnd - tStart);
    gsl_vector_free(states);
}

void pawan::__integration::step(const double &dt, __system *S, gsl_vector* states){
    gsl_vector *rates = gsl_vector_calloc(S->_maxsize);  //ip: large memory allocations would make code slow here
    S->solve();
    S->getRates(rates);
    gsl_vector_scale(rates,dt);
    gsl_vector_add(states,rates);
    S->setStates(states);
    gsl_vector_free(rates);
}

void pawan::__integration::integrate(__system *S, __io *IO, __rotor_ll *R){
    DOUT("-----------------------in pawan::__integration::integrate()");
	gsl_vector *states = gsl_vector_calloc(S->_size);
	FILE *f = IO->create_binary_file(".wake");
	double t = 0.0;
    fwrite(&_n,sizeof(size_t),1,R->_rotorf);  //write number of time steps to rotor motion file
	fwrite(&t,sizeof(double),1,f);
	S->write(f);
	S->getStates(states);
	double tStart = TIME();
	for(size_t i = 1; i<=_n; ++i){
		OUT("\tStep",i);
		t = i*_dt;
		step(_dt,S,states,R);
		fwrite(&t,sizeof(double),1,f);	
		S->write(f);
	}
	fclose(f);
	double tEnd = TIME();
	OUT("Total Time (s)",tEnd - tStart);
	gsl_vector_free(states);
}

void pawan::__integration::step(const double &dt,
                                __system *S,
                                gsl_vector* states,
                                __rotor_ll *R){
    DOUT("-----------------------in pawan::__integration::step()");
	gsl_vector *rates = gsl_vector_calloc(S->_size);
	R->update_motion(dt);
	S->solve();
	S->getRates(rates);
	gsl_vector_scale(rates,dt);
	gsl_vector_add(states,rates);
	S->setStates(states);
	gsl_vector_free(rates);
    R->update_wake();
}
