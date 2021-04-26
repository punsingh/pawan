/*! PArticle Wake ANalysis
 * \file parallel.h
 * \brief Header for OpenMP Parallel Interaction class
 *
 * @author Puneet Singh
 * @date 04/24/2021
 */
#ifndef PARALLEL_H_
#define PARALLEL_H_

#include <stdio.h>
#include <iostream>
#include <gsl/gsl_vector.h>
#include "src/utils/print_utils.h"
#include "src/wake/wake.h"
#include "src/interaction/interaction.h"
#include "src/interaction/interaction_utils.h"

namespace pawan{
class __parallel : public __interaction{

	private:
		//! Interact
		/*
		 * Compute interaction between particles of a single wake object
		 */
		virtual void interact(__wake *W);

	public:
		
		//! Constructor
		/*
		 * Creates empty integral
		 * \param	W	Wake object pointer
		 */
		__parallel(__wake *W);
		
		//! Destructor
		/*
		 * Deletes particles
		 */
		~__parallel() = default;
};
}
#endif
