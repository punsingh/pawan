/*! PArticle Wake ANalysis
 * \file system.h
 * \brief Header for System class
 *
 * @author Puneet Singh
 * @date 06/24/2021
 */
#ifndef SYSTEM_H_
#define SYSTEM_H_

#include <stdio.h>
#include <iostream>
#include <vector>

#include <gsl/gsl_vector.h>
#include "src/utils/print_utils.h"
#include "src/networkinterface/networkdatastructures.h"

namespace pawan{
class __system{

	public:
		size_t _size;		/*!< Size of state vector */
		size_t _maxsize;
        size_t _maxnumParticles;
		//! Constructor
		/*
		 * Creates empty system object
		 */
		__system();
		
		//! Destructor
		/*
		 * Delete system
		 */
		~__system() = default;

		//! Solvel system 
		/*
		 * Solve system
		 */
		virtual void solve(){};

		//! Write all system data
		/*
		 * Write binary file with all system data
		 * \param f	Binary file
		 */
		virtual void write(FILE *f){};

		//! Set states
		/*
		 * Sets all system states
		 * \param state		System state
		 */
		virtual void setStates(const gsl_vector *state){};
		
		//! Get rates
		/*
		 * Get system rates
		 * \param rate		System rate
		 */
		virtual void getRates(gsl_vector *rate){};
		
		//! Get states
		/*
		 * Get current System states
		 * \param state		System state
		 */
		virtual void getStates(gsl_vector *state){};

        //! add new particles
        virtual void addParticles(PawanRecvData pawanrecvdata){};
        //!
        virtual void updateVinfEffect(double &dt, gsl_vector* states){};
        virtual void updateVinfEffect(double &dt){};

};
}
#endif
