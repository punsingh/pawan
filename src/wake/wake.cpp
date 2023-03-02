/*! PArticle Wake ANalysis
 * \file wake.cpp
 * \brief Routines for Wake calculations
 *
 * @author Puneet Singh
 * @date 03/28/2021
 */
#include "wake.h"
#define SHEDVOR 1 //model only trail vortices if 0

void pawan::__wake::create_particles(const int &n){
    DOUT("------------------------------------------------------in pawan::__wake::create_particles()");
	_numParticles    = n;
	_size            = 2*_numParticles*3;
    _maxnumParticles = _numParticles;
    _maxsize         = 2*_maxnumParticles*3;
	_position        = gsl_matrix_calloc(_numParticles,3);
	_velocity        = gsl_matrix_calloc(_numParticles,3);
	_vorticity       = gsl_matrix_calloc(_numParticles,3);
	_retvorcity      = gsl_matrix_calloc(_numParticles,3);
	_radius          = gsl_vector_calloc(_numParticles);
	_volume          = gsl_vector_calloc(_numParticles);
	_birthstrength   = gsl_vector_calloc(_numParticles);
	_vorticityfield  = gsl_matrix_calloc(_numParticles,3);
    _active          = gsl_vector_calloc(_numParticles);
}

void pawan::__wake::initialise_memory(){
    DOUT("------------------------------------------------------in pawan::__wake::create_particles()");
    _numParticles = 0; //no particle info in arrays yet
    _size = 2*_numParticles*3; //this becomes irrelevant in the case of large memory allocation
    _maxnumParticles = 10000;
    _maxsize         = 2*_maxnumParticles*3;
    _position        = gsl_matrix_calloc(_maxnumParticles,3);
    _velocity        = gsl_matrix_calloc(_maxnumParticles,3);
    _vorticity       = gsl_matrix_calloc(_maxnumParticles,3);
    _retvorcity      = gsl_matrix_calloc(_maxnumParticles,3);
    _radius          = gsl_vector_calloc(_maxnumParticles);
    _volume          = gsl_vector_calloc(_maxnumParticles);
    _birthstrength   = gsl_vector_calloc(_maxnumParticles);
    _vorticityfield  = gsl_matrix_calloc(_maxnumParticles,3);
    _active          = gsl_vector_calloc(_maxnumParticles);
}

pawan::__wake::__wake(){
    DOUT("------------------------------------------------------in pawan::__wake::__wake()");
	create_particles(1);
	for(int i = 0; i<_numParticles; ++i){
		for(int j = 0; j<3; ++j){
			gsl_matrix_set(_position,i,j,0.0);
			gsl_matrix_set(_vorticity,i,j,0.0);
			gsl_matrix_set(_vorticityfield,i,j,0.0);
		}
		gsl_vector_set(_radius,i,1.0);
		gsl_vector_set(_volume,i,1.0);
		gsl_vector_set(_birthstrength,i,1.0);
	}
}

void pawan::__wake::addParticles(PawanRecvData pawanrecvdata){
    double *Vinf = pawanrecvdata->Vinf;
    int NbOfLfnLines = pawanrecvdata->NbOfLfnLines;
    int *NbOfAst = pawanrecvdata->NbOfAst;
    double hres = pawanrecvdata->hres;
    double acrossa = pawanrecvdata->acrossa;
    double *TEpos = pawanrecvdata->TEpos;
    double *circ = pawanrecvdata->circ;
    double *TEpos_prev = pawanrecvdata->TEpos_prev;
    double *circ_prev = pawanrecvdata->circ_prev;
    double *astpos = pawanrecvdata->astpos;

    double Vinfmag = sqrt(Vinf[0]*Vinf[0] + Vinf[1]*Vinf[1] + Vinf[2]*Vinf[2]);
    int npidx = _numParticles; //new particle index start
    int astidx = 0;
    for(size_t ilfn = 0; ilfn<NbOfLfnLines; ++ilfn){
        for(size_t iast = 0; iast<NbOfAst[ilfn]; ++iast) {

            //evaluate trailing vortex parameters corresponding to a airsta
            double trailvorpos[3];
            trailvorpos[0] = TEpos[3 * astidx];
            trailvorpos[1] = TEpos[3 * astidx + 1];
            trailvorpos[2] = TEpos[3 * astidx + 2];

            double localflowdir[3];  //need more astkin details to get this accurately for blades
            localflowdir[0] = Vinf[0] / Vinfmag;
            localflowdir[1] = Vinf[1] / Vinfmag;
            localflowdir[2] = Vinf[2] / Vinfmag;

            double trailvormag; //magnitude and sign of trail vortices
            if (iast == 0)                   //if beginning (from root) of a new lfnline
                trailvormag = -circ[astidx];//clockwise vortex when looking from behind is -ve
            else if (iast == NbOfAst[ilfn] - 1)//last airstation gives anticlockwise vortex
                trailvormag = circ[astidx];      //anticlockwise vortex is +ve
            else
                trailvormag = circ[astidx] - circ[astidx + 1];

            for (int k = 0; k < 3; ++k) {
                gsl_matrix_set(_position, npidx, k, trailvorpos[k]);
                gsl_matrix_set(_vorticity, npidx, k, localflowdir[k] * trailvormag);
                gsl_matrix_set(_vorticityfield, npidx, k, 0.0); //not sure about this
            }
            gsl_vector_set(_radius, npidx, 0.01);
            gsl_vector_set(_volume, npidx, 0.01);
            gsl_vector_set(_birthstrength, npidx, trailvormag);
            gsl_vector_set(_active, npidx, 1.0);
            npidx++;

            if (SHEDVOR) {
                //evaluate shed vortex parameters corresponding to a airsta
                double shedvorpos[3];
                double shedvordir[3];
                double shedvormag; //magnitude of shed vortices
                double shedvordirmag;

                if (iast != NbOfAst[ilfn] - 1) {//shed vortex only between airsta
                    shedvordir[0] = TEpos[3 * astidx] - TEpos[3 * (astidx + 1)];
                    shedvordir[1] = TEpos[3 * astidx + 1] - TEpos[3 * (astidx + 1) + 1];
                    shedvordir[2] = TEpos[3 * astidx + 2] - TEpos[3 * (astidx + 1) + 2];

                    shedvordirmag = sqrt(shedvordir[0]*shedvordir[0] + shedvordir[1]*shedvordir[1] + shedvordir[2]*shedvordir[2]);

                    shedvorpos[0] = (TEpos[3 * astidx] + TEpos[3 * (astidx + 1)]) / 2;
                    shedvorpos[1] = (TEpos[3 * astidx + 1] + TEpos[3 * (astidx + 1) + 1]) / 2;
                    shedvorpos[2] = (TEpos[3 * astidx + 2] + TEpos[3 * (astidx + 1) + 2]) / 2;

                    shedvormag = circ[astidx] - circ_prev[astidx];

                    for (int k = 0; k < 3; ++k) {
                        gsl_matrix_set(_position, npidx, k, shedvorpos[k]);
                        gsl_matrix_set(_vorticity, npidx, k, shedvordir[k] * shedvormag/shedvordirmag);
                        gsl_matrix_set(_vorticityfield, npidx, k, 0.0); //not sure about this
                    }
                    gsl_vector_set(_radius, npidx, 0.01);
                    gsl_vector_set(_volume, npidx, 0.01);
                    gsl_vector_set(_birthstrength, npidx, shedvormag);
                    gsl_vector_set(_active, npidx, 1.0);

                    npidx++;
                }
            }
            astidx++;
        }
    }
    _numParticles = npidx;
    _size = 2*_numParticles*3;
}
pawan::__wake::__wake(PawanRecvData pawanrecvdata){

    initialise_memory();   //large enough memory allocation
    addParticles(pawanrecvdata);
}
void pawan::__wake::updateVinfEffect(double &dt, gsl_vector* states){
    double Vinf[3];
    Vinf[0] = -1000;
    Vinf[1] = 0;
    Vinf[2] = 0;
    //printf("Vinf: %+8.3e, %+8.3e, %+8.3e\n", Vinf[0], Vinf[1], Vinf[2]);
    int idx;
    for(int i = 0; i<_numParticles; ++i) {
        for (int j = 0; j < 3; ++j) {
            idx = i*3+j;
            if(idx==0)
                printf("Position of 1st particle before update: %+8.3e\n",
                       gsl_vector_get(states, 0));
            gsl_vector_set(states, idx, Vinf[j]*dt + gsl_vector_get(states, idx));
            if(idx==0)
                printf("Position of 1st particle after update: %+8.3e\n",
                       gsl_vector_get(states, 0));
        }
    }
}

void pawan::__wake::updateVinfEffect(double &dt){
    double Vinf[3];
    Vinf[0] = -10;
    Vinf[1] = 0;
    Vinf[2] = 0;
    //printf("Vinf: %+8.3e, %+8.3e, %+8.3e\n", Vinf[0], Vinf[1], Vinf[2]);
    for(int i = 0; i<_numParticles; ++i) {
        for (int j = 0; j < 3; ++j) {
            if(i==0 && j==0)
                printf("Position of 1st particle     before update: %+8.3e\n",
                       gsl_matrix_get(_position, 0, 0));
            gsl_matrix_set(_position, i, j, Vinf[j]*dt + gsl_matrix_get(_position, i, j));
            if(i==0 && j==0)
                printf("Position of 1st particle     after update: %+8.3e\n",
                       gsl_matrix_get(_position, 0, 0));
        }
    }
}

pawan::__wake::~__wake(){
    DOUT("------------------------------------------------------in pawan::__wake::~__wake()");
	gsl_matrix_free(_position);
	gsl_matrix_free(_velocity);
	gsl_matrix_free(_vorticity);
	gsl_matrix_free(_retvorcity);
	gsl_vector_free(_radius);
	gsl_vector_free(_volume);
	gsl_vector_free(_birthstrength);
	gsl_matrix_free(_vorticityfield);
    gsl_vector_free(_active);
}

void pawan::__wake::print(){
    DOUT("------------------------------------------------------in pawan::__wake::print()");
	OUT("_numParticles",_numParticles);
	OUT("_position",_position);
	OUT("_vorticity",_vorticity);
	OUT("_radius",_radius);
	OUT("_volume",_volume);
	OUT("_birthstrength",_birthstrength);
	OUT("_vorticityfield",_vorticityfield);
    OUT("_active",_active);
}

void pawan::__wake::save(FILE *f){
    DOUT("------------------------------------------------------in pawan::__wake::save()");
	fwrite(&_numParticles,sizeof(size_t),1,f);
	gsl_matrix_fwrite(f,_position);
	gsl_matrix_fwrite(f,_vorticity);
	gsl_vector_fwrite(f,_radius);
	gsl_vector_fwrite(f,_volume);
	gsl_vector_fwrite(f,_birthstrength);
	//gsl_matrix_fwrite(f,_vorticityfield);
}

void pawan::__wake::write(FILE *f){
    DOUT("------------------------------------------------------in pawan::__wake::write()");
	fwrite(&_numParticles,sizeof(size_t),1,f);
	gsl_matrix_fwrite(f,_position);
	gsl_matrix_fwrite(f,_vorticity);
	gsl_vector_fwrite(f,_radius);
    gsl_vector_fwrite(f,_active);
}

void pawan::__wake::read(__io *IO){
    DOUT("------------------------------------------------------in pawan::__wake::read()");
	FILE *f = IO->open_binary_file(".wake");
	fread(&_numParticles,sizeof(size_t),1,f);	
	gsl_matrix_fread(f,_position);
	gsl_matrix_fread(f,_vorticity);
	gsl_vector_fread(f,_radius);
	gsl_vector_fread(f,_volume);
	gsl_vector_fread(f,_birthstrength);
	//gsl_matrix_fread(f,_vorticityfield);
	fclose(f);
}

void pawan::__wake::setStates(const gsl_vector *state){
    DOUT("------------------------------------------------------in pawan::__wake::setStates()");
	//OUT("setStates");
	//size_t np = state->size/2/3;
	size_t offset = _maxnumParticles*3;
	size_t vecsize = _maxnumParticles*3;
	gsl_matrix_const_view pos = gsl_matrix_const_view_vector(state,_maxnumParticles,3);
	gsl_matrix_memcpy(_position,&pos.matrix);
	gsl_vector_const_view vor = gsl_vector_const_subvector(state,offset,vecsize);
	gsl_matrix_const_view vrx = gsl_matrix_const_view_vector(&vor.vector,_maxnumParticles,3);
	gsl_matrix_memcpy(_vorticity,&vrx.matrix);
}

void pawan::__wake::getRates(gsl_vector *rate){
    DOUT("------------------------------------------------------in pawan::__wake::getRates()");
	//OUT("getRates");
	for(size_t i = 0; i<_numParticles; ++i){
		for(size_t j = 0; j<3; ++j){
			size_t ind = i*3 + j;
			gsl_vector_set(rate,ind,gsl_matrix_get(_velocity,i,j));
			ind += _maxsize/2;
			gsl_vector_set(rate,ind,gsl_matrix_get(_retvorcity,i,j));
		}
	}
}

void pawan::__wake::getStates(gsl_vector *state){
    DOUT("------------------------------------------------------in pawan::__wake::getStates()");
	//OUT("getStates");
	for(size_t i = 0; i<_numParticles; ++i){
		for(size_t j = 0; j<3; ++j){
			size_t ind = i*3 + j;
			gsl_vector_set(state,ind,gsl_matrix_get(_position,i,j));
			ind += _maxsize/2;
			gsl_vector_set(state,ind,gsl_matrix_get(_vorticity,i,j));
		}
	}
}

void pawan::__wake::getIdealRates(gsl_vector *rate){
    DOUT("------------------------------------------------------in pawan::__wake::getIdealRates()");
	//OUT("getIdealRates");
	gsl_vector_set_zero(rate);
}

void pawan::__wake::translate(const size_t &n, const double &x){
    DOUT("------------------------------------------------------in pawan::__wake::translate()");
	for(int i = 0; i<_numParticles; ++i){
		gsl_matrix_set(_position,i,n,x + gsl_matrix_get(_position,i,n));
	}
}
