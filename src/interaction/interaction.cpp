/*! PArticle Wake ANalysis
 * \file interaction.cpp
 * \brief Routines for interactions
 *
 * @author Puneet Singh
 * @date 04/24/2021
 */
#include "interaction.h"
#include "interaction_utils.h"

pawan::__interaction::__interaction(__wake *W){
	_W = W;
	_size = _W->_size;
	_nu = 2.0e-5;
}

void pawan::__interaction::interact(){
	interact(_W);
}

void pawan::__interaction::interact(__wake *W){
	gsl_matrix_set_zero(W->_velocity);
	gsl_matrix_set_zero(W->_retvorcity);
	for(size_t i_src = 0; i_src < W->_numParticles; ++i_src){
		gsl_vector_const_view r_src = gsl_matrix_const_row(W->_position,i_src);
		gsl_vector_const_view a_src = gsl_matrix_const_row(W->_vorticity,i_src);
		gsl_vector_view dr_src = gsl_matrix_row(W->_velocity,i_src);
		gsl_vector_view da_src = gsl_matrix_row(W->_retvorcity,i_src);
		double s_src = gsl_vector_get(W->_radius,i_src);
		double v_src = gsl_vector_get(W->_volume,i_src);
		for(size_t i_trg = i_src + 1; i_trg < W->_numParticles; ++i_trg){
			gsl_vector_const_view r_trg = gsl_matrix_const_row(W->_position,i_trg);
			gsl_vector_const_view a_trg= gsl_matrix_const_row(W->_vorticity,i_trg);
			gsl_vector_view dr_trg = gsl_matrix_row(W->_velocity,i_trg);
			gsl_vector_view da_trg = gsl_matrix_row(W->_retvorcity,i_trg);
			double s_trg = gsl_vector_get(W->_radius,i_trg);
			double v_trg = gsl_vector_get(W->_volume,i_trg);
			INTERACT(_nu,s_src,s_trg,&r_src.vector,&r_trg.vector,&a_src.vector,&a_trg.vector,v_src,v_trg,&dr_src.vector,&dr_trg.vector,&da_src.vector,&da_trg.vector);
		}
	}
}

void pawan::__interaction::write(FILE *f){
	_W->write(f);
}

void pawan::__interaction::setStates(const gsl_vector *state){
	_W->setStates(state);
}

void pawan::__interaction::getRates(gsl_vector *rate){
	_W->getRates(rate);
}

void pawan::__interaction::getStates(gsl_vector *state){
	_W->getStates(state);
}
