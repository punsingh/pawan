/*! PArticle Wake ANalysis
 * \file integration.h
 * \brief Header for Integration class
 *
 * @author Puneet Singh
 * @date 04/21/2021
 */
#ifndef INTEGRATION_H_
#define INTEGRATION_H_

#include <stdio.h>
#include <iostream>
#include <gsl/gsl_vector.h>
#include "src/utils/print_utils.h"
#include "src/wake/wake.h"
#include "src/io/io.h"

namespace pawan{
class __integration{

	private:
		double _dt;	/*!< Time step size */
		double _t;	/*!< Total time */
		size_t _n;	/*!< Number of time steps */

	public:
		//! Constructor
		/*
		 * Creates empty integral
		 * \param	t	Total time
		 * \param	n	Number of time steps
		 */
		__integration(const double &t, const size_t &n);
		
		//! Destructor
		/*
		 * Deletes particles
		 */
		~__integration() = default;
		
		//! Integrate
		/*
		 * Integrates wake
		 * \param	W	Wake object
		 * \param	IO	Input/Output file writing
		 */
		void integrate(__wake *W, __io *IO);
};
}
#endif
