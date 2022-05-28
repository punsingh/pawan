/*! PArticle Wake ANalysis
 * \file resolve.h
 * \brief Header for Resolve class
 *
 * @author Puneet Singh
 * @date 09/12/2021
 */
#ifndef RESOLVE_H_
#define RESOLVE_H_

#include <stdio.h>
#include <iostream>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include "src/utils/print_utils.h"
#include "src/utils/timing_utils.h"
#include "src/interaction/interaction.h"
#include "src/io/io.h"

namespace pawan{
class __resolve{

	protected:

	public:
		//! Constructor
		/*
		 * Default constructor
		 */
		__resolve() = default;
		
		//! Destructor
		/*
		 * Deletes particles
		 */
		~__resolve() = default;
		
		//! Calculates vorticity field and modifies states
		/*
		 * Integrates wake
		 * \param	IN	Interaction solver
		 * \param	IO	Input/Output file writing
		 */
		void rebuild(__interaction *IN, __io *IO);
};
}
#endif
