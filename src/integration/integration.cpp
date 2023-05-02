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

void pawan::__integration::integrate(__system *S,
                                     __io *IO,
                                     NetworkInterfaceTCP<OPawanRecvData,OPawanSendData> *networkCommunicatorTest,
                                     bool coupling,
                                     bool diagnose){
    gsl_vector *states = gsl_vector_calloc(S->_totalmaxsize); //ip: have the *rates pointer here outside too for speed up
    FILE *f = IO->create_binary_file(".wake");
    FILE *fdiag = IO->create_binary_file(".diagnosis");
    double t = 0.0;
    fwrite(&t,sizeof(double),1,f);
    S->write(f);  //write particles info as is
    if(diagnose) {
        S->writenu(fdiag);
        fwrite(&t,sizeof(double),1,fdiag);
        S->diagnose();
        S->writediagnosis(fdiag);
    }
    double tStart = TIME();
    for(size_t i = 1; i<=_n; ++i){
        S->getStates(states);
        OUT("\tStep",i);
        t = i*_dt;
        step(_dt,S,states);
        //S->relax();

        fwrite(&t,sizeof(double),1,f);
        S->write(f);  //write particles info after interaction in this time step

        if(diagnose){
            S->diagnose();
            fwrite(&t,sizeof(double),1,fdiag);
            S->writediagnosis(fdiag);
        }

        if (coupling) {
            OPawanRecvData opawanrecvdata;//create it once outside the loop
            networkCommunicatorTest->getrecieveBuffer(opawanrecvdata);
            S->updateVinfEffect(opawanrecvdata.Vinf,_dt);
            //S->updateBoundVorEffect(&opawanrecvdata,_dt);

            OPawanSendData opawansenddata;//create it once outside the loop and should be good
            S->getInflow(&opawanrecvdata,&opawansenddata);
            networkCommunicatorTest->send_data(opawansenddata);
            if (i != _n) {
                //S->diagnose();
                OPawanRecvData opawanrecvdata;
                networkCommunicatorTest->recieve_data(opawanrecvdata);
                S->addParticles(&opawanrecvdata);
            }
        }
    }
    fclose(f);
    double tEnd = TIME();
    OUT("Total Time (s)",tEnd - tStart);
    gsl_vector_free(states);
}

void pawan::__integration::step(const double &dt, __system *S, gsl_vector* states){
    gsl_vector *rates = gsl_vector_calloc(S->_totalmaxsize);  //ip: large memory allocations would make code slow here
    S->solve();
    S->getRates(rates);
    gsl_vector_scale(rates,dt);
    printf("Position change of 1st particle     __integration::step(): %+8.3e, %+8.3e, %+8.3e\n",
           gsl_vector_get(rates, 0),gsl_vector_get(rates, 1),gsl_vector_get(rates, 2));
    gsl_vector_add(states,rates);
    S->setStates(states);
    gsl_vector_free(rates);
}

