/*! PArticle Wake ANalysis
 * \file interaction.cpp
 * \brief Routines for interactions
 *
 * @author Puneet Singh
 * @date 04/24/2021
 */
#include "interaction.h"
#include "interaction_utils.h"

pawan::__interaction::__interaction(){
	//_nu = 2.0e-2;
	_nu = 2.5e-3;
	_nWake = 0;
	_totalVorticity = gsl_vector_calloc(3);
	_linearImpulse = gsl_vector_calloc(3);
	_angularImpulse = gsl_vector_calloc(3);
}

pawan::__interaction::__interaction(__wake *W):__interaction(){
	addWake(W);
}

pawan::__interaction::__interaction(__wake *W1, __wake *W2):__interaction(){
	addWake(W1);
	addWake(W2);
}

pawan::__interaction::~__interaction(){
	gsl_vector_free(_totalVorticity);
	gsl_vector_free(_linearImpulse);
	gsl_vector_free(_angularImpulse);
}

void pawan::__interaction::addWake(__wake *W){
	_W.push_back(W);
	_size += W->_size;
	_nWake++;
}

void pawan::__interaction::solve(){
	interact();
}

void pawan::__interaction::resolve(){
	influence();
}

void pawan::__interaction::diagnose(){
	// Linear diagnostics
	gsl_vector_set_zero(_totalVorticity);
	gsl_vector_set_zero(_linearImpulse);
	gsl_vector_set_zero(_angularImpulse);
	
	gsl_vector *O = gsl_vector_calloc(3);
	gsl_vector *I = gsl_vector_calloc(3);
	gsl_vector *A = gsl_vector_calloc(3);
	for(auto &w : _W){
		// Total vorticity
		calculateTotalVorticity(w,O);
		gsl_vector_add(_totalVorticity,O);
		// Linear impulse
		calculateLinearImpulse(w,I);
		gsl_vector_add(_linearImpulse,I);
		// Angular impulse
		calculateAngularImpulse(w,A);
		gsl_vector_add(_angularImpulse,A);
	}
	OUTT("Total Vorticity",_totalVorticity);
	OUTT("Linear Impulse",_linearImpulse);
	OUTT("Angular Impulse",_angularImpulse);
	gsl_vector_free(O);
	gsl_vector_free(I);
	gsl_vector_free(A);

	// Quadratic diagnostics
	_helicity = 0.0;
	for(size_t i = 0; i<_nWake; ++i){
		_helicity += calculateHelicity(_W[i]);
		for(size_t j = i + 1; j < _nWake; ++j){
			_helicity += calculateHelicity(_W[i],_W[j]);
		}
	}
	OUT("Helicity",_helicity);
	_enstrophy = 0.0;
	for(size_t i = 0; i<_nWake; ++i){
		_enstrophy += calculateEnstrophy(_W[i]);
		for(size_t j = i + 1; j < _nWake; ++j){
			_enstrophy += calculateEnstrophy(_W[i],_W[j]);
		}
	}
	OUT("Enstrophy",_enstrophy);
}

void pawan::__interaction::interact(){
	for(auto &w : _W){
		gsl_matrix_set_zero(w->_velocity);
		gsl_matrix_set_zero(w->_retvorcity);
	}
	for(auto &w : _W){
		interact(w);
	}
	for(size_t i = 0; i<_nWake; ++i){
		for(size_t j = i + 1; j < _nWake; ++j){
			interact(_W[i],_W[j]);
		}
	}
}

void pawan::__interaction::influence(){
	for(auto &w : _W){
		gsl_matrix_set_zero(w->_vorticityfield);
	}
	for(auto &w : _W){
		influence(w);
	}
	for(size_t i = 0; i<_nWake; ++i){
		for(size_t j = i + 1; j < _nWake; ++j){
			influence(_W[i],_W[j]);
		}
	}
}

void pawan::__interaction::influence(__wake *W){
	for(size_t i_src = 0; i_src < W->_numParticles; ++i_src){
		gsl_vector_const_view r_src = gsl_matrix_const_row(W->_position,i_src);
		gsl_vector_const_view a_src = gsl_matrix_const_row(W->_vorticity,i_src);
		gsl_vector_view k_src = gsl_matrix_row(W->_vorticityfield,i_src);
		double s_src = gsl_vector_get(W->_radius,i_src);
		SELFINFLUENCE(s_src,&r_src.vector,&a_src.vector,&k_src.vector);
		for(size_t i_trg = i_src + 1; i_trg < W->_numParticles; ++i_trg){
			gsl_vector_const_view r_trg = gsl_matrix_const_row(W->_position,i_trg);
			gsl_vector_const_view a_trg= gsl_matrix_const_row(W->_vorticity,i_trg);
			gsl_vector_view k_trg = gsl_matrix_row(W->_vorticityfield,i_trg);
			double s_trg = gsl_vector_get(W->_radius,i_trg);
			INFLUENCE(s_src,s_trg,&r_src.vector,&r_trg.vector,&a_src.vector,&a_trg.vector,&k_src.vector,&k_trg.vector);
		}
	}
}

void pawan::__interaction::influence(__wake *W1, __wake *W2){
	for(size_t i_src = 0; i_src < W1->_numParticles; ++i_src){
		gsl_vector_const_view r_src = gsl_matrix_const_row(W1->_position,i_src);
		gsl_vector_const_view a_src = gsl_matrix_const_row(W1->_vorticity,i_src);
		gsl_vector_view k_src = gsl_matrix_row(W1->_vorticityfield,i_src);
		double s_src = gsl_vector_get(W1->_radius,i_src);
		for(size_t i_trg = 0; i_trg < W2->_numParticles; ++i_trg){
			gsl_vector_const_view r_trg = gsl_matrix_const_row(W2->_position,i_trg);
			gsl_vector_const_view a_trg = gsl_matrix_const_row(W2->_vorticity,i_trg);
			gsl_vector_view k_trg = gsl_matrix_row(W2->_vorticityfield,i_trg);
			double s_trg = gsl_vector_get(W2->_radius,i_trg);
			INFLUENCE(s_src,s_trg,&r_src.vector,&r_trg.vector,&a_src.vector,&a_trg.vector,&k_src.vector,&k_trg.vector);
		}
	}
}

void pawan::__interaction::interact(__wake *W){
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

void pawan::__interaction::interact(__wake *W1, __wake *W2){
	for(size_t i_src = 0; i_src < W1->_numParticles; ++i_src){
		gsl_vector_const_view r_src = gsl_matrix_const_row(W1->_position,i_src);
		gsl_vector_const_view a_src = gsl_matrix_const_row(W1->_vorticity,i_src);
		gsl_vector_view dr_src = gsl_matrix_row(W1->_velocity,i_src);
		gsl_vector_view da_src = gsl_matrix_row(W1->_retvorcity,i_src);
		double s_src = gsl_vector_get(W1->_radius,i_src);
		double v_src = gsl_vector_get(W1->_volume,i_src);
		for(size_t i_trg = 0; i_trg < W2->_numParticles; ++i_trg){
			gsl_vector_const_view r_trg = gsl_matrix_const_row(W2->_position,i_trg);
			gsl_vector_const_view a_trg = gsl_matrix_const_row(W2->_vorticity,i_trg);
			gsl_vector_view dr_trg = gsl_matrix_row(W2->_velocity,i_trg);
			gsl_vector_view da_trg = gsl_matrix_row(W2->_retvorcity,i_trg);
			double s_trg = gsl_vector_get(W2->_radius,i_trg);
			double v_trg = gsl_vector_get(W2->_volume,i_trg);
			INTERACT(_nu,s_src,s_trg,&r_src.vector,&r_trg.vector,&a_src.vector,&a_trg.vector,v_src,v_trg,&dr_src.vector,&dr_trg.vector,&da_src.vector,&da_trg.vector);
		}
	}
}

void pawan::__interaction::write(FILE *f){
	fwrite(&_nWake,sizeof(size_t),1,f);	
	for(auto &w: _W){
		w->write(f);
	}
}

void pawan::__interaction::setStates(const gsl_vector *state){
	size_t offset = 0;
	for(auto &w: _W){
		gsl_vector_const_view st = gsl_vector_const_subvector(state,offset,w->_size);
		w->setStates(&st.vector);
		offset += w->_size;
	}
}

void pawan::__interaction::getRates(gsl_vector *rate){
	size_t offset = 0;
	for(auto &w: _W){
		gsl_vector_view rt = gsl_vector_subvector(rate,offset,w->_size);
		w->getRates(&rt.vector);
		offset += w->_size;
	}
}

void pawan::__interaction::getStates(gsl_vector *state){
	size_t offset = 0;
	for(auto &w: _W){
		gsl_vector_view st = gsl_vector_subvector(state,offset,w->_size);
		w->getStates(&st.vector);
		offset += w->_size;
	}
}

void pawan::__interaction::getIdealRates(gsl_vector *rate){
	size_t offset = 0;
	for(auto &w: _W){
		gsl_vector_view rt = gsl_vector_subvector(rate,offset,w->_size);
		w->getIdealRates(&rt.vector);
		offset += w->_size;
	}
}

void pawan::__interaction::calculateTotalVorticity(__wake *W, gsl_vector *O){
	gsl_vector_set_zero(O);
	for(size_t i = 0; i < W->_numParticles; ++i){
		gsl_vector_const_view a = gsl_matrix_const_row(W->_vorticity,i);
		gsl_vector_add(O,&a.vector);
	}
}

void pawan::__interaction::calculateLinearImpulse(__wake *W, gsl_vector *I){
	gsl_vector_set_zero(I);
	gsl_vector *rxa = gsl_vector_calloc(3);
	for(size_t i = 0; i < W->_numParticles; ++i){
		gsl_vector_const_view r = gsl_matrix_const_row(W->_position,i);
		gsl_vector_const_view a = gsl_matrix_const_row(W->_vorticity,i);
		gsl_cross(&r.vector,&a.vector,rxa);
		gsl_vector_add(I,rxa);
	}
	gsl_vector_scale(I,0.5);
	gsl_vector_free(rxa);
}

void pawan::__interaction::calculateAngularImpulse(__wake *W, gsl_vector *A){
	gsl_vector_set_zero(A);
	gsl_vector *rxa = gsl_vector_calloc(3);
	gsl_vector *rxrxa = gsl_vector_calloc(3);
	for(size_t i = 0; i < W->_numParticles; ++i){
		gsl_vector_const_view r = gsl_matrix_const_row(W->_position,i);
		gsl_vector_const_view a = gsl_matrix_const_row(W->_vorticity,i);
		gsl_cross(&r.vector,&a.vector,rxa);
		gsl_cross(&r.vector,rxa,rxrxa);
		gsl_vector_add(A,rxrxa);
	}
	gsl_vector_scale(A,1.0/3.0);
	gsl_vector_free(rxrxa);
	gsl_vector_free(rxa);
}

double pawan::__interaction::calculateEnstrophy(__wake *W){
	double s = 0.0;
	for(size_t I = 0; I < W->_numParticles; ++I){
		gsl_vector_const_view rI = gsl_matrix_const_row(W->_position,I);
		gsl_vector_const_view aI = gsl_matrix_const_row(W->_vorticity,I);
		double sI = gsl_vector_get(W->_radius,I);
		s += ENSTROPHY(sI,&aI.vector);
		//double ens = ENSTROPHY(sI,&aI.vector);
		//s += ens;
		//OUT("e",ens);
		for(size_t J = I + 1; J < W->_numParticles; ++J){
			//OUT("\tJ",J);
			gsl_vector_const_view rJ = gsl_matrix_const_row(W->_position,J);
			gsl_vector_const_view aJ = gsl_matrix_const_row(W->_vorticity,J);
			double sJ = gsl_vector_get(W->_radius,J);
			// doubling S because S(i,j) = S(j,i)
			s += 2*ENSTROPHY(sI,sJ,&rI.vector,&rJ.vector,&aI.vector,&aJ.vector);
			//OUT("\ts",s);
		}
	}
	return s;
}

double pawan::__interaction::calculateEnstrophy(__wake *W1, __wake *W2){
	double s = 0.0;
	for(size_t I = 0; I < W1->_numParticles; ++I){
		gsl_vector_const_view rI = gsl_matrix_const_row(W1->_position,I);
		gsl_vector_const_view aI = gsl_matrix_const_row(W1->_vorticity,I);
		double sI = gsl_vector_get(W1->_radius,I);
		for(size_t J = 0; J < W2->_numParticles; ++J){
			gsl_vector_const_view rJ = gsl_matrix_const_row(W2->_position,J);
			gsl_vector_const_view aJ = gsl_matrix_const_row(W2->_vorticity,J);
			double sJ = gsl_vector_get(W2->_radius,J);
			s += ENSTROPHY(sI,sJ,&rI.vector,&rJ.vector,&aI.vector,&aJ.vector);
		}
	}
	// doubling S because S(i,j) = S(j,i)
	return 2.0*s;
}

double pawan::__interaction::calculateHelicity(__wake *W){
	double h = 0.0;
	for(size_t I = 0; I < W->_numParticles; ++I){
		gsl_vector_const_view rI = gsl_matrix_const_row(W->_position,I);
		gsl_vector_const_view aI = gsl_matrix_const_row(W->_vorticity,I);
		double sI = gsl_vector_get(W->_radius,I);
		for(size_t J = I + 1; J < W->_numParticles; ++J){
			gsl_vector_const_view rJ = gsl_matrix_const_row(W->_position,J);
			gsl_vector_const_view aJ = gsl_matrix_const_row(W->_vorticity,J);
			double sJ = gsl_vector_get(W->_radius,J);
			h += HELICITY(sI,sJ,&rI.vector,&rJ.vector,&aI.vector,&aJ.vector);
		}
	}
	// doubling H because H(i,j) = H(j,i)
	return 2.0*h;
}

double pawan::__interaction::calculateHelicity(__wake *W1, __wake *W2){
	double h = 0.0;
	for(size_t I = 0; I < W1->_numParticles; ++I){
		gsl_vector_const_view rI = gsl_matrix_const_row(W1->_position,I);
		gsl_vector_const_view aI = gsl_matrix_const_row(W1->_vorticity,I);
		double sI = gsl_vector_get(W1->_radius,I);
		for(size_t J = 0; J < W2->_numParticles; ++J){
			gsl_vector_const_view rJ = gsl_matrix_const_row(W2->_position,J);
			gsl_vector_const_view aJ = gsl_matrix_const_row(W2->_vorticity,J);
			double sJ = gsl_vector_get(W2->_radius,J);
			h += HELICITY(sI,sJ,&rI.vector,&rJ.vector,&aI.vector,&aJ.vector);
		}
	}
	// doubling H because H(i,j) = H(j,i)
	return 2.0*h;
}

double pawan::__interaction::calculateKineticEnergy(__wake *W){
	double ke = 0.0;
	for(size_t I = 0; I < W->_numParticles; ++I){
		gsl_vector_const_view rI = gsl_matrix_const_row(W->_position,I);
		gsl_vector_const_view aI = gsl_matrix_const_row(W->_vorticity,I);
		double sI = gsl_vector_get(W->_radius,I);
		ke += KINETICENERGY(sI,&aI.vector);
		for(size_t J = I + 1; J < W->_numParticles; ++J){
			gsl_vector_const_view rJ = gsl_matrix_const_row(W->_position,J);
			gsl_vector_const_view aJ = gsl_matrix_const_row(W->_vorticity,J);
			double sJ = gsl_vector_get(W->_radius,J);
			ke += 2.0*KINETICENERGY(sI,sJ,&rI.vector,&rJ.vector,&aI.vector,&aJ.vector);
		}
	}
	return 2.0*ke;
}

double pawan::__interaction::calculateKineticEnergy(__wake *W1, __wake *W2){
	double ke = 0.0;
	for(size_t I = 0; I < W1->_numParticles; ++I){
		gsl_vector_const_view rI = gsl_matrix_const_row(W1->_position,I);
		gsl_vector_const_view aI = gsl_matrix_const_row(W1->_vorticity,I);
		double sI = gsl_vector_get(W1->_radius,I);
		for(size_t J = 0; J < W2->_numParticles; ++J){
			gsl_vector_const_view rJ = gsl_matrix_const_row(W2->_position,J);
			gsl_vector_const_view aJ = gsl_matrix_const_row(W2->_vorticity,J);
			double sJ = gsl_vector_get(W2->_radius,J);
			ke += KINETICENERGY(sI,sJ,&rI.vector,&rJ.vector,&aI.vector,&aJ.vector);
		}
	}
	// doubling KE because KE(i,j) = KE(j,i)
	return 2.0*ke;
}
