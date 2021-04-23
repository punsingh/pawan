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
		size_t _numParticles;		/*!< Number of vortex particles */
		size_t _numDimensions;		/*!< Number of dimensions */
		gsl_matrix *_position;		/*!< Particle positions */
		gsl_matrix *_velocity;		/*!< Particle velocity */
		gsl_matrix *_vorticity;		/*!< Particle vorticities */
		gsl_matrix *_retvorcity;	/*!< Particle rate of change of vorticity*/
		gsl_vector *_radius;		/*!< Particle smoothing radii */
		gsl_vector *_volume;		/*!< Particle volumes */
		gsl_vector *_birthstrength;	/*!< Strengths of particles at birth */
		
		//! Create particles
		/*
		 * Creates empty particles
		 * \param n	Number of particles
		 */
		void create_particles(const int &n);

	public:
		size_t _size;			/*!< Size of state vector */
		
		//! Constructor
		/*
		 * Creates default
		 */
		__wake();
		
		//! Constructor for thin vortex ring
		/*
		 * Creates thin vortex ring
		 * \param gamma		Strength of ring	
		 * \param radius	Ring radius
		 * \param core		Core radius
		 * \param nRadial	Number of particles
		 */
		__wake(const double &gamma, const double &radius, const double &core, const int &nRadial);
		
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
		 * \param f	Binary file
		 */
		virtual void write(FILE *f);

		//! Read wake data file
		/*
		 * Read binary file with all wake particle data
		 * \param op	Input/Output object for file
		 */
		virtual void read(__io *IO);

		//! Calculate wake interactions
		/*
		 * Calculate velocities and rate of change of vorticities for all wake particles
		 */
		virtual void calculateInteraction();
		
		//! Set states
		/*
		 * Sets positions and vorticities of particles
		 * \param state		Wake state
		 */
		virtual void setStates(const gsl_vector *state);
		
		//! Get rates
		/*
		 * Get velocities and rate of change of vorticities of particles
		 * \param rate		Wake rate
		 */
		virtual void getRates(gsl_vector *rate);
		
		//! Get states
		/*
		 * Get position and vorticities of particles
		 * \param rate		Wake state
		 */
		virtual void getStates(gsl_vector *state);

};
}
#endif
