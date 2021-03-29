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
		gsl_matrix *_position;
		gsl_matrix *_vorticity;

	public:

		__wake();
		~__wake();

		virtual void print();

};

#endif
