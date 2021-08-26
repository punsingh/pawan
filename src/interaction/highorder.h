/*! PArticle Wake ANalysis
 * \file highorder.h
 * \brief Inline functions for High Order Algebraic smoothing
 *
 * @author 	Puneet Singh
 * @date	08/20/2021
 *
 */

#ifndef HIGHORDER_H_
#define HIGHORDER_H_

#include "src/utils/gsl_utils.h"

/*! \fn inline double ZETASIG(const double &rho, const double &sigma, double &q, double &F, double &Z)
 * \brief Compute vortex particle kernel
 * \param	rho		double distance
 * \param	sigma		double radius
 */
inline double ZETASIG(	const double &rho, const double &sigma){
	double rho_bar = rho/sigma;
	return 1.875*M_1_PI/pow(rho_bar*rho_bar + 1.0,7.5)/pow(sigma,3);
};

/*! \fn inline void KERNEL(const double &rho, const double &sigma, double &q, double &F, double &Z)
 * \brief Compute velocity induced by vortex particle kernel
 * \param	rho		double distance
 * \param	sigma		double radius
 * \param	q		double Q
 * \param	F		double F
 * \param	Z		double Z
 */
inline void KERNEL(	const double &rho, 
			const double &sigma, 
			double &q, 
			double &F, 
			double &Z){
	double rho_bar = rho/sigma;
	Z = ZETASIG(rho,sigma);
	double rho_bar2 = gsl_pow_2(rho_bar);
	double sig3 = sigma*sigma*sigma;
	double phi = 0.25*M_1_PI*(rho_bar2 + 1.5)/pow(rho_bar2 + 1.0,2.0)/sig3;
	q = (phi/rho_bar - Z)/gsl_pow_2(rho_bar);
	F = (Z - 3*q)/gsl_pow_2(rho);
};

#endif
