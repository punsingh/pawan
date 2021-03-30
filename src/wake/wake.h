/*! PArticle Wake ANalysis
 * \file wake.h
 * \brief Header for WAKE lass
 *
 * @author Puneet Singh
 * @date 03/28/2021
 */
#ifndef WAKE_H_
#define WAKE_H_

#include <stdio.h>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include "src/utils/gsl_utils.h"

class __wake{

	private:
		size_t _numParticles;
		gsl_matrix *_position;
		gsl_matrix *_vorticity;
		gsl_vector *_radius;
		gsl_vector *_volume;
		gsl_vector *_birthstrength;

	public:
		//! Constructor
		/*
		 * Creates random particles
		 */
		__wake();
		
		//! Copy Constructor
		/*
		 * Copies particles
		 */
		__wake(const __wake &w);
		
		//! Destructor
		/*
		 * Creates random particles
		 */
		~__wake();

		//! Print all wake particles
		/*
		 * Creates random particles
		 */
		virtual void print();

};

#endif
