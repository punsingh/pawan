/*! PArticle Wake ANalysis
 * \file interaction.h
 * \brief Header for Interaction class
 *
 * @author Puneet Singh
 * @date 04/24/2021
 */
#ifndef INTERACTION_H_
#define INTERACTION_H_

#include <stdio.h>
#include <iostream>
#include <gsl/gsl_vector.h>
#include "src/utils/print_utils.h"
#include "src/wake/wake.h"

namespace pawan{
class __interaction{

	protected:
		double _nu;		/*!< Kinematic viscosity */

	private:
		pawan::__wake *_W;	/*!< Pointer to wake object */
		
		//! Interact
		/*
		 * Compute interaction between particles of a single wake object
		 */
		virtual void interact(__wake *W);

	public:
		size_t _size;		/*!< Size of state vector */
		
		//! Constructor
		/*
		 * Creates empty integral
		 * \param	W	Wake object pointer
		 */
		__interaction(__wake *W);
		
		//! Destructor
		/*
		 * Deletes particles
		 */
		~__interaction() = default;
		
		//! Interact
		/*
		 * Compute interaction between wake particles
		 */
		void interact();

		//! Write all wake data
		/*
		 * Write binary file with all wake particle data
		 * \param f	Binary file
		 */
		virtual void write(FILE *f);

		//! Set states
		/*
		 * Sets all wake states
		 * \param state		Wake state
		 */
		virtual void setStates(const gsl_vector *state);
		
		//! Get rates
		/*
		 * Get wake state rates
		 * \param rate		Wake rate
		 */
		virtual void getRates(gsl_vector *rate);
		
		//! Get states
		/*
		 * Get current wake states
		 * \param state		Wake state
		 */
		virtual void getStates(gsl_vector *state);
};
}
#endif
