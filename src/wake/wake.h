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

#include "src/io/io.h"
#include "src/utils/gsl_utils.h"

namespace pawan{

class __wake{

	private:
		std::string _file;		/*!< Name and path of wake data file*/
		
		size_t _numParticles;		/*!< Number of vortex particles */
		size_t _numDimensions;		/*!< Number of dimensions */
		gsl_matrix *_position;		/*!< Particle positions */
		gsl_matrix *_vorticity;		/*!< Particle vorticities */
		gsl_vector *_radius;		/*!< Particle smoothing radii */
		gsl_vector *_volume;		/*!< Particle volumes */
		gsl_vector *_birthstrength;	/*!< Strengths of particles at birth */

	public:
		//! Constructor
		/*
		 * Creates random particles
		 * \param op	Input/Output object
		 */
		__wake(__io *op);
		
		//! Destructor
		/*
		 * Deletes particles
		 */
		~__wake();

		//! Print all wake particles
		/*
		 * Print wake particle information
		 */
		virtual void print();

		//! Write wake data file
		/*
		 * Write binary file with all wake particle data
		 */
		virtual void write();

		//! Read wake data file
		/*
		 * Read binary file with all wake particle data
		 */
		virtual void read();

};
}
#endif
