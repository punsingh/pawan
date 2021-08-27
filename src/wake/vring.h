/*! PArticle Wake ANalysis
 * \file vring.h
 * \brief Header for Vortex Ring example class
 *
 * @author Puneet Singh
 * @date 08/12/2021
 */
#ifndef VRING_H_
#define VRING_H_

#include <stdio.h>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>

#include "src/wake/wake.h"
#include "src/wake/ring.h"

namespace pawan{

class __vring : public __ring{

	private:

	public:
		//! Constructor for vogel model vortex ring
		/*
		 * Creates vortex ring
		 * \param radius	Ring radius
		 * \param core		Core radius
		 * \param nRadial	Number of radial sections
		 * \param nCore		Number of core layers
		 * \param shift		Index from zero or halfpoint
		 */
		__vring(const double &radius, const double &core, const int &nRadial, const int &nCore, const bool &shift = false);
		
		//! Destructor
		/*
		 * Deletes particles
		 */
		~__vring();
};
}
#endif
