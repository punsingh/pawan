/*! PArticle Wake ANalysis
 * \file wake.cpp
 * \brief Routines for Wake calculations
 *
 * @author Puneet Singh
 * @date 03/28/2021
 */
#include "wake.h"

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
    _maxnumParticles = MAXNUMPARTICLES;
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

void pawan::__wake::split(){

    size_t npsplit = 0;//number of particles added by splitting
    for(size_t i = 0; i<_numParticles; ++i){
        double ibirthstrength = gsl_vector_get(_birthstrength,i);
        gsl_vector_view ipos = gsl_matrix_row(_position,i);
        gsl_vector_view ivor = gsl_matrix_row(_vorticity,i);
        double ivormag = gsl_blas_dnrm2(&ivor.vector);
        double iradius = gsl_vector_get(_radius,i);
        double ivol = gsl_vector_get(_volume,i);

        //can this be parallelised when all threads could be changing
        // _vorticity and _position at the same time?
        if(ivormag>=2*ibirthstrength && ibirthstrength>0.0001) {//ip: hard-coded, a more appropriate condition needs to be found
            printf("splitting!!!--------- %d --------- ", i);
            gsl_vector *deltapos = gsl_vector_calloc(3);
            gsl_vector_memcpy(deltapos,&ivor.vector);
            gsl_blas_dscal(0.25*iradius/ivormag,deltapos);

            //add a particle at the end of array
            size_t splitidx = _numParticles+npsplit;
            gsl_vector_view ipos_split = gsl_matrix_row(_position,splitidx);
            gsl_vector_memcpy(&ipos_split.vector,&ipos.vector);
            gsl_vector_add(&ipos_split.vector,deltapos);
            gsl_vector_view ivor_split = gsl_matrix_row(_vorticity,splitidx);
            gsl_vector_memcpy(&ivor_split.vector,&ivor.vector);
            gsl_blas_dscal(0.5,&ivor_split.vector);
            gsl_vector_set(_radius,splitidx,iradius);
            gsl_vector_set(_volume,splitidx,ivol);  //need to confirm this
            gsl_vector_set(_birthstrength,splitidx,0.5*ivormag);

            //replace particle parameters in-place
            gsl_vector_sub(&ipos.vector,deltapos);
            gsl_blas_dscal(0.5,&ivor.vector);
            gsl_vector_set(_birthstrength,i,0.5*ivormag);

            gsl_vector_free(deltapos);
            npsplit++;
        }
    }
    _numParticles = _numParticles + npsplit;
    _size = 2*_numParticles*3;
}

void pawan::__wake::merge(){

    int npmerge = 0;//number of particles 'lost' due to merging
    for(size_t i = 0; i<_numParticles; ++i){
        gsl_vector_view ipos = gsl_matrix_row(_position,i);
        gsl_vector_view ivor = gsl_matrix_row(_vorticity,i);
        double ivormag = gsl_blas_dnrm2(&ivor.vector);
        double iradius = gsl_vector_get(_radius,i);

        for(size_t j = i + 1; j < _numParticles; ++j) {
            gsl_vector_view jpos = gsl_matrix_row(_position,j);
            gsl_vector_view jvor = gsl_matrix_row(_vorticity,j);
            double jvormag = gsl_blas_dnrm2(&jvor.vector);
            double jradius = gsl_vector_get(_radius,j);

            gsl_vector *displacement = gsl_vector_calloc(3);
            gsl_vector_memcpy(displacement,&ipos.vector);
            gsl_vector_sub(displacement,&jpos.vector);
            double rho = gsl_blas_dnrm2(displacement);
            if (rho<0.125*(iradius+jradius) && jvormag!=0.0 && ivormag!=0.0) {
                double alignment=0.0;
                gsl_blas_ddot(&ivor.vector,&jvor.vector,&alignment);
                alignment = 1 - alignment/ivormag/jvormag;
                if(alignment>=0.0001){//DBL_EPSILON??
                    printf("merging!!!--------- %d ----and---- %d--------- ", i,j);
                    gsl_vector *ipos_m = gsl_vector_calloc(3);
                    gsl_vector *tmp = gsl_vector_calloc(3);
                    gsl_vector_memcpy(ipos_m,&ipos.vector);
                    gsl_blas_dscal(ivormag,ipos_m);
                    gsl_vector_memcpy(tmp,&jpos.vector);
                    gsl_blas_dscal(jvormag,tmp);
                    gsl_vector_add(ipos_m,tmp);
                    gsl_blas_dscal(ivormag+jvormag,ipos_m);
                    gsl_vector_memcpy(&ipos.vector,ipos_m);    //in-place replacement
                    gsl_vector_add(&ivor.vector,&jvor.vector);

                    gsl_vector_set_zero(&jvor.vector);     //"deleting" the other particle
                    gsl_vector_set(_volume,j,0.0);

                    gsl_vector_free(ipos_m);
                    gsl_vector_free(tmp);
                }
            }
            gsl_vector_free(displacement);
        }
    }
}

void pawan::__wake::relax(){

    split();
    merge();
}

void pawan::__wake::getlfnvec(double* vec,
                             const double* mat,
                             const int rowsize,
                             const int axis,
                             const int offsetrows,
                             const int size){
    int offset = offsetrows*rowsize;
    for(size_t i = 0; i<size; ++i){
        vec[i] = mat[offset + rowsize * i + axis];
    }
}

void pawan::__wake::addParticles(PawanRecvData pawanrecvdata){
    OUT("_numParticles",_numParticles);
    double *Vinf = pawanrecvdata->Vinf;
    int NbOfLfnLines = pawanrecvdata->NbOfLfnLines;
    int *NbOfAst = pawanrecvdata->NbOfAst;
    double *lfnlen = pawanrecvdata->lfnlen;
    //int spanres = pawanrecvdata->spanres;
    int spanres = 15;  //number of particles over the span of a lfnline
    double acrossa = pawanrecvdata->acrossa;
    double deltat = pawanrecvdata->deltat;
    double *spandisc = pawanrecvdata->span_disc;
    double *TEpos = pawanrecvdata->TEpos;
    double *circ = pawanrecvdata->circ;
    double *TEpos_prev = pawanrecvdata->TEpos_prev;
    double *circ_prev = pawanrecvdata->circ_prev;
    double *astpos = pawanrecvdata->astpos;
    double *astvel = pawanrecvdata->astvel;

    double Vinfmag = sqrt(Vinf[0]*Vinf[0] + Vinf[1]*Vinf[1] + Vinf[2]*Vinf[2]);
    int npidx = _numParticles; //new particle index start
    int astidx = 0;
    int offsetrows=0;
    for(size_t ilfn = 0; ilfn<NbOfLfnLines; ++ilfn) {

        double *TEposx_lfn = (double *) calloc(NbOfAst[ilfn], sizeof(double));
        double *TEposy_lfn = (double *) calloc(NbOfAst[ilfn], sizeof(double));
        double *TEposz_lfn = (double *) calloc(NbOfAst[ilfn], sizeof(double));
        getlfnvec(TEposx_lfn, TEpos, 3, 0, offsetrows, NbOfAst[ilfn]);
        getlfnvec(TEposy_lfn, TEpos, 3, 1, offsetrows, NbOfAst[ilfn]);
        getlfnvec(TEposz_lfn, TEpos, 3, 2, offsetrows, NbOfAst[ilfn]);
        double *TEposxprev_lfn = (double *) calloc(NbOfAst[ilfn], sizeof(double));
        double *TEposyprev_lfn = (double *) calloc(NbOfAst[ilfn], sizeof(double));
        double *TEposzprev_lfn = (double *) calloc(NbOfAst[ilfn], sizeof(double));
        getlfnvec(TEposxprev_lfn, TEpos_prev, 3, 0, offsetrows, NbOfAst[ilfn]);
        getlfnvec(TEposyprev_lfn, TEpos_prev, 3, 1, offsetrows, NbOfAst[ilfn]);
        getlfnvec(TEposzprev_lfn, TEpos_prev, 3, 2, offsetrows, NbOfAst[ilfn]);
        double *circ_lfn = (double *) calloc(NbOfAst[ilfn], sizeof(double));
        getlfnvec(circ_lfn, circ, 1, 0, offsetrows, NbOfAst[ilfn]);
        double *circprev_lfn = (double *) calloc(NbOfAst[ilfn], sizeof(double));
        getlfnvec(circprev_lfn, circ_prev, 1, 0, offsetrows, NbOfAst[ilfn]);
        double *astposx_lfn = (double *) calloc(NbOfAst[ilfn], sizeof(double));
        double *astposy_lfn = (double *) calloc(NbOfAst[ilfn], sizeof(double));
        double *astposz_lfn = (double *) calloc(NbOfAst[ilfn], sizeof(double));
        getlfnvec(astposx_lfn, astpos, 3, 0, offsetrows, NbOfAst[ilfn]);
        getlfnvec(astposy_lfn, astpos, 3, 1, offsetrows, NbOfAst[ilfn]);
        getlfnvec(astposz_lfn, astpos, 3, 2, offsetrows, NbOfAst[ilfn]);
        double *astvelx_lfn = (double *) calloc(NbOfAst[ilfn], sizeof(double));
        double *astvely_lfn = (double *) calloc(NbOfAst[ilfn], sizeof(double));
        double *astvelz_lfn = (double *) calloc(NbOfAst[ilfn], sizeof(double));
        getlfnvec(astvelx_lfn, astvel, 3, 0, offsetrows, NbOfAst[ilfn]);
        getlfnvec(astvely_lfn, astvel, 3, 1, offsetrows, NbOfAst[ilfn]);
        getlfnvec(astvelz_lfn, astvel, 3, 2, offsetrows, NbOfAst[ilfn]);
        double *spandisc_lfn = (double *) calloc(NbOfAst[ilfn], sizeof(double));
        getlfnvec(spandisc_lfn, spandisc, 1, 0, offsetrows, NbOfAst[ilfn]);
        offsetrows += NbOfAst[ilfn];//also just astidx can be used

        gsl_interp_accel *accel_ptr;
        gsl_spline *circ_spline, *TEposx_spline, *TEposy_spline, *TEposz_spline;
        gsl_spline *circprev_spline, *TEposxprev_spline, *TEposyprev_spline, *TEposzprev_spline;
        gsl_spline *astposx_spline, *astposy_spline, *astposz_spline;
        gsl_spline *astvelx_spline, *astvely_spline, *astvelz_spline;
        accel_ptr = gsl_interp_accel_alloc();
        circ_spline = gsl_spline_alloc(gsl_interp_cspline, NbOfAst[ilfn]);
        TEposx_spline = gsl_spline_alloc(gsl_interp_cspline, NbOfAst[ilfn]);
        TEposy_spline = gsl_spline_alloc(gsl_interp_cspline, NbOfAst[ilfn]);
        TEposz_spline = gsl_spline_alloc(gsl_interp_cspline, NbOfAst[ilfn]);
        circprev_spline = gsl_spline_alloc(gsl_interp_cspline, NbOfAst[ilfn]);
        TEposxprev_spline = gsl_spline_alloc(gsl_interp_cspline, NbOfAst[ilfn]);
        TEposyprev_spline = gsl_spline_alloc(gsl_interp_cspline, NbOfAst[ilfn]);
        TEposzprev_spline = gsl_spline_alloc(gsl_interp_cspline, NbOfAst[ilfn]);
        astposx_spline = gsl_spline_alloc(gsl_interp_cspline, NbOfAst[ilfn]);
        astposy_spline = gsl_spline_alloc(gsl_interp_cspline, NbOfAst[ilfn]);
        astposz_spline = gsl_spline_alloc(gsl_interp_cspline, NbOfAst[ilfn]);
        astvelx_spline = gsl_spline_alloc(gsl_interp_cspline, NbOfAst[ilfn]);
        astvely_spline = gsl_spline_alloc(gsl_interp_cspline, NbOfAst[ilfn]);
        astvelz_spline = gsl_spline_alloc(gsl_interp_cspline, NbOfAst[ilfn]);
        gsl_spline_init(circ_spline, spandisc_lfn, circ_lfn, NbOfAst[ilfn]);
        gsl_spline_init(TEposx_spline, spandisc_lfn, TEposx_lfn, NbOfAst[ilfn]);
        gsl_spline_init(TEposy_spline, spandisc_lfn, TEposy_lfn, NbOfAst[ilfn]);
        gsl_spline_init(TEposz_spline, spandisc_lfn, TEposz_lfn, NbOfAst[ilfn]);
        gsl_spline_init(circprev_spline, spandisc_lfn, circprev_lfn, NbOfAst[ilfn]);
        gsl_spline_init(TEposxprev_spline, spandisc_lfn, TEposxprev_lfn, NbOfAst[ilfn]);
        gsl_spline_init(TEposyprev_spline, spandisc_lfn, TEposyprev_lfn, NbOfAst[ilfn]);
        gsl_spline_init(TEposzprev_spline, spandisc_lfn, TEposzprev_lfn, NbOfAst[ilfn]);
        gsl_spline_init(astposx_spline, spandisc_lfn, astposx_lfn, NbOfAst[ilfn]);
        gsl_spline_init(astposy_spline, spandisc_lfn, astposy_lfn, NbOfAst[ilfn]);
        gsl_spline_init(astposz_spline, spandisc_lfn, astposz_lfn, NbOfAst[ilfn]);
        gsl_spline_init(astvelx_spline, spandisc_lfn, astvelx_lfn, NbOfAst[ilfn]);
        gsl_spline_init(astvely_spline, spandisc_lfn, astvely_lfn, NbOfAst[ilfn]);
        gsl_spline_init(astvelz_spline, spandisc_lfn, astvelz_lfn, NbOfAst[ilfn]);

        //get peak circulation and index
        int iast_peakcirc_lfn = 0;          //airsta index where peak circulation occurs
        double peakcirc_lfn = circ_lfn[0];  //peak lfnline circulation
        double spandisc_peakcirc_lfn = 0.0; //non-dim radius where peak circ occurs
        for (int i = 1; i < NbOfAst[ilfn]; i++) {
            if (circ[i] > peakcirc_lfn) {
                peakcirc_lfn = circ_lfn[i];
                iast_peakcirc_lfn = i;
                spandisc_peakcirc_lfn = spandisc_lfn[i];
            }
        }

        double sigma = 0.0;
        double hres = lfnlen[ilfn]/spanres;
        if(SHEDVOR)
            hres =  lfnlen[ilfn] / (2 * spanres - 1);
        else
            hres = lfnlen[ilfn] / spanres;
        sigma = 1.5 * hres;

        for (int ista = 0; ista < spanres + 1; ++ista) {//new discretized station index
            double tepos[3];
            tepos[0] = gsl_spline_eval(TEposx_spline, ista * (1. / spanres), accel_ptr);
            tepos[1] = gsl_spline_eval(TEposy_spline, ista * (1. / spanres), accel_ptr);
            tepos[2] = gsl_spline_eval(TEposz_spline, ista * (1. / spanres), accel_ptr);
            double teposprev[3];
            teposprev[0] = gsl_spline_eval(TEposxprev_spline, ista * (1. / spanres), accel_ptr);
            teposprev[1] = gsl_spline_eval(TEposyprev_spline, ista * (1. / spanres), accel_ptr);
            teposprev[2] = gsl_spline_eval(TEposzprev_spline, ista * (1. / spanres), accel_ptr);
            double teconv[3];//convection velocity of vparticles at TE
            /*teconv[0]=gsl_spline_eval(astvelx_spline, ista * (1. / spanres), accel_ptr);
            teconv[1]=gsl_spline_eval(astvely_spline, ista * (1. / spanres), accel_ptr);
            teconv[2]=gsl_spline_eval(astvelz_spline, ista * (1. / spanres), accel_ptr);
            */
            teconv[0] = Vinf[0]; teconv[1] = Vinf[1]; teconv[2] = Vinf[2];
            //teconv[0] = 0; teconv[1] = 0; teconv[2] = 0;
            double del_teposmag;
            //del_teposmag = sqrt(pow(tepos[0]-teposprev[0],2) + pow(tepos[1]-teposprev[1],2) + pow(tepos[2]-teposprev[2],2));
            del_teposmag = sqrt(pow(tepos[0]-(teposprev[0]+deltat*teconv[0]),2)
                                    + pow(tepos[1]-(teposprev[1]+deltat*teconv[1]),2)
                                    + pow(tepos[2]-(teposprev[2]+deltat*teconv[2]),2));

            //double QUANT = del_teposmag/hres;

            double vx = gsl_spline_eval(astvelx_spline, ista * (1. / spanres), accel_ptr);
            double vy = gsl_spline_eval(astvely_spline, ista * (1. / spanres), accel_ptr);
            double vz = gsl_spline_eval(astvelz_spline, ista * (1. / spanres), accel_ptr);
            double vmag = sqrt(vx * vx + vy * vy + vz * vz);
            double localflowdir[3];  //need more astkin details to get this accurately for blades
            localflowdir[0] = vx / vmag;
            localflowdir[1] = vy / vmag;
            localflowdir[2] = vz / vmag;

            double trailvormag; //magnitude and sign of trail vortices
            //assuming only a single peak in the circulation distribution for now
            if ((double) ista / spanres <= spandisc_peakcirc_lfn && ista != 0)
                trailvormag = -(gsl_spline_eval(circ_spline, ista * (1. / spanres), accel_ptr)
                                - gsl_spline_eval(circ_spline, (ista - 1) * (1. / spanres), accel_ptr));
            else if ((double) ista / spanres > spandisc_peakcirc_lfn && ista != spanres)
                trailvormag = gsl_spline_eval(circ_spline, ista * (1. / spanres), accel_ptr)
                              - gsl_spline_eval(circ_spline, (ista + 1) * (1. / spanres), accel_ptr);

            //trailvormag=0;
            if (ista == 0)
                trailvormag = -circ_lfn[0];//if beginning (from root) of a new lfnline
                                           //clockwise vortex when looking from behind is -ve
            if (ista == spanres)
                trailvormag = circ_lfn[NbOfAst[ilfn] - 1];//last airstation gives anticlockwise vortex
                                                          //anticlockwise vortex is +ve
            double alphai = trailvormag*hres;
            double voli = 0.05*hres;   //ip: not entirely sure about this (is this an empirical parameter??)

            double trailvorpos[3];
            for(int jsta = 0; jsta<=(int)(del_teposmag/hres); ++jsta ) {
                //printf("spline = %10.5e %10.5e %10.5e %10.5e \n", ista, ista*1, ista*2.0, ista*2.0/hres);
                //evaluate trailing vortex parameters corresponding to a airsta
                if (del_teposmag != 0.0) {//in case lfnline is moving in inertial frame
                    trailvorpos[0] = tepos[0] + jsta * ((teposprev[0]+deltat*teconv[0]) - tepos[0]) / (del_teposmag / hres);
                    trailvorpos[1] = tepos[1] + jsta * ((teposprev[1]+deltat*teconv[1]) - tepos[1]) / (del_teposmag / hres);
                    trailvorpos[2] = tepos[2] + jsta * ((teposprev[2]+deltat*teconv[2]) - tepos[2]) / (del_teposmag / hres);
                } else {//in case lfnline is stationary in inertial frame
                    trailvorpos[0] = tepos[0];
                    trailvorpos[1] = tepos[1];
                    trailvorpos[2] = tepos[2];
                }
                for (int k = 0; k < 3; ++k) {
                    gsl_matrix_set(_position, npidx, k, trailvorpos[k]);
                    gsl_matrix_set(_vorticity, npidx, k, localflowdir[k] * alphai);
                    gsl_matrix_set(_vorticityfield, npidx, k, 0.0); //later
                }
                gsl_vector_set(_radius, npidx, sigma);
                gsl_vector_set(_volume, npidx, voli);
                gsl_vector_set(_birthstrength, npidx, fabs(alphai));
                gsl_vector_set(_active, npidx, 1.0);//unsused for now
                if (npidx == 0)
                    printf("First particle created at : %8.3e, %8.3e, %8.3e \n", trailvorpos[0], trailvorpos[1],trailvorpos[2]);
                npidx++;
            }
            if (SHEDVOR) {
                //evaluate shed vortex parameters corresponding to a airsta
                double shedvorpos[3];
                double shedvordir[3];
                double shedvormag; //magnitude of shed vortices
                double shedvordirmag;

                if (ista != spanres) {//shed vortex only between airsta
                    double tepos_nxt[3];
                    tepos_nxt[0] = gsl_spline_eval(TEposx_spline, (ista+1) * (1. / spanres), accel_ptr);
                    tepos_nxt[1] = gsl_spline_eval(TEposy_spline, (ista+1) * (1. / spanres), accel_ptr);
                    tepos_nxt[2] = gsl_spline_eval(TEposz_spline, (ista+1) * (1. / spanres), accel_ptr);
                    double teposprev_nxt[3];
                    teposprev_nxt[0] = gsl_spline_eval(TEposxprev_spline, (ista+1) * (1. / spanres), accel_ptr);
                    teposprev_nxt[1] = gsl_spline_eval(TEposyprev_spline, (ista+1) * (1. / spanres), accel_ptr);
                    teposprev_nxt[2] = gsl_spline_eval(TEposzprev_spline, (ista+1) * (1. / spanres), accel_ptr);

                    shedvordir[0] =   gsl_spline_eval(astposx_spline, ista * (1. / spanres), accel_ptr)
                                    - gsl_spline_eval(astposx_spline, (ista+1) * (1. / spanres), accel_ptr);
                    shedvordir[1] =   gsl_spline_eval(astposy_spline, ista * (1. / spanres), accel_ptr)
                                    - gsl_spline_eval(astposy_spline, (ista+1) * (1. / spanres), accel_ptr);
                    shedvordir[2] =   gsl_spline_eval(astposz_spline, ista * (1. / spanres), accel_ptr)
                                    - gsl_spline_eval(astposz_spline, (ista+1) * (1. / spanres), accel_ptr);

                    shedvordirmag = sqrt(shedvordir[0]*shedvordir[0] + shedvordir[1]*shedvordir[1] + shedvordir[2]*shedvordir[2]);

                    //ip: not exactly correct (since last station not reached?)
//                    shedvormag =   gsl_spline_eval(circ_spline,     ista * (1. / hres), accel_ptr)
//                                 - gsl_spline_eval(circprev_spline, ista * (1. / hres), accel_ptr);
                    shedvormag =   0.5*(   gsl_spline_eval(circ_spline,     ista *     (1. / spanres),     accel_ptr)
                                         + gsl_spline_eval(circ_spline,     (ista+1) * (1. / spanres),     accel_ptr))
                                 - 0.5*(   gsl_spline_eval(circprev_spline, ista *     (1. / spanres),     accel_ptr)
                                         + gsl_spline_eval(circprev_spline, (ista+1) * (1. / spanres),     accel_ptr));
                    //placing shed vortex particle between two trailing vortex particles
                    for(int ksta = 0; ksta<=(int)(del_teposmag/hres); ++ksta ) {
                        if (del_teposmag != 0.0) {//in case lfnline is moving in inertial frame
                            shedvorpos[0] = 0.5*(   tepos[0]     + ksta * ((teposprev[0]    + deltat*teconv[0]) - tepos[0]    ) / (del_teposmag / hres)
                                                  + tepos_nxt[0] + ksta * ((teposprev_nxt[0]+ deltat*teconv[0]) - tepos_nxt[0]) / (del_teposmag / hres));
                            shedvorpos[1] = 0.5*(   tepos[1]     + ksta * ((teposprev[1]    + deltat*teconv[1]) - tepos[1]    ) / (del_teposmag / hres)
                                                    + tepos_nxt[1] + ksta * ((teposprev_nxt[1]+ deltat*teconv[1]) - tepos_nxt[1]) / (del_teposmag / hres));
                            shedvorpos[2] = 0.5*(   tepos[2]     + ksta * ((teposprev[2]    + deltat*teconv[2]) - tepos[2]    ) / (del_teposmag / hres)
                                                    + tepos_nxt[2] + ksta * ((teposprev_nxt[2]+ deltat*teconv[2]) - tepos_nxt[2]) / (del_teposmag / hres));

                        } else {//in case lfnline is stationary in inertial frame
                            shedvorpos[0] = 0.5*(  gsl_spline_eval(TEposx_spline, ista * (1. / spanres), accel_ptr)
                                                   + gsl_spline_eval(TEposx_spline, (ista+1) * (1. / spanres), accel_ptr));
                            shedvorpos[1] = 0.5*(  gsl_spline_eval(TEposy_spline, ista * (1. / spanres), accel_ptr)
                                                   + gsl_spline_eval(TEposy_spline, (ista+1) * (1. / spanres), accel_ptr));
                            shedvorpos[2] = 0.5*(  gsl_spline_eval(TEposz_spline, ista * (1. / spanres), accel_ptr)
                                                   + gsl_spline_eval(TEposz_spline, (ista+1) * (1. / spanres), accel_ptr));
                        }
                        double alphai_shed = shedvormag*hres;
                        for (int k = 0; k < 3; ++k) {
                            gsl_matrix_set(_position, npidx, k, shedvorpos[k]);
                            gsl_matrix_set(_vorticity, npidx, k, shedvordir[k] * alphai_shed/shedvordirmag);
                            gsl_matrix_set(_vorticityfield, npidx, k, 0.0); //later
                        }
                        gsl_vector_set(_radius, npidx, sigma);
                        gsl_vector_set(_volume, npidx, voli);
                        gsl_vector_set(_birthstrength, npidx, fabs(alphai_shed));
                        gsl_vector_set(_active, npidx, 1.0);//unsused for now

                        npidx++;
                    }
                }
            }
        }
        gsl_spline_free(circ_spline);
        gsl_spline_free(TEposx_spline);
        gsl_spline_free(TEposy_spline);
        gsl_spline_free(TEposz_spline);
        gsl_spline_free(circprev_spline);
        gsl_spline_free(TEposxprev_spline);
        gsl_spline_free(TEposyprev_spline);
        gsl_spline_free(TEposzprev_spline);
        gsl_spline_free(astposx_spline);
        gsl_spline_free(astposy_spline);
        gsl_spline_free(astposz_spline);
        gsl_interp_accel_free(accel_ptr);
        free(TEposx_lfn);free(TEposy_lfn);free(TEposz_lfn);
        free(TEposxprev_lfn);free(TEposyprev_lfn);free(TEposzprev_lfn);
        free(circ_lfn);free(circprev_lfn);free(spandisc_lfn);
        free(astposx_lfn);free(astposy_lfn);free(astposz_lfn);
        free(astvelx_lfn);free(astvely_lfn);free(astvelz_lfn);

    }
    _numParticles = npidx;
    _size = 2*_numParticles*3;
}
pawan::__wake::__wake(PawanRecvData pawanrecvdata){

    initialise_memory();   //large enough memory allocation
    addParticles(pawanrecvdata);
}

void pawan::__wake::addParticles(__wake *W){

    for(int i = 0; i<W->_numParticles; ++i){
        for(int j = 0; j<3; ++j){
            gsl_matrix_set(_position,i,j,gsl_matrix_get(W->_position,i,j));
            gsl_matrix_set(_vorticity,i,j,gsl_matrix_get(W->_vorticity,i,j));
            gsl_matrix_set(_vorticityfield,i,j,gsl_matrix_get(W->_vorticityfield,i,j));
        }
        gsl_vector_set(_radius,i,gsl_vector_get(W->_radius,i));
        gsl_vector_set(_volume,i,gsl_vector_get(W->_volume,i));
        gsl_vector_set(_birthstrength,i,gsl_vector_get(W->_birthstrength,i));
    }
    _numParticles = W->_numParticles;
    _size = 2*_numParticles*3;
}

pawan::__wake::__wake(__wake *W){

    initialise_memory();   //large enough memory allocation
    addParticles(W);
}

void pawan::__wake::updateVinfEffect(const double *Vinf, double &dt){

    for(int i = 0; i<_numParticles; ++i) {
        for (int j = 0; j < 3; ++j) {
            if(i==0 && j==0)
                printf("Position of 1st particle     before update: %+8.3e, %+8.3e, %+8.3e\n",
                       gsl_matrix_get(_position, 0, 0),gsl_matrix_get(_position, 0, 1),gsl_matrix_get(_position, 0, 2));
            gsl_matrix_set(_position, i, j, Vinf[j]*dt + gsl_matrix_get(_position, i, j));
            if(i==0 && j==0)
                printf("Position of 1st particle     after update: %+8.3e, %+8.3e, %+8.3e\n",
                       gsl_matrix_get(_position, 0, 0),gsl_matrix_get(_position, 0, 1),gsl_matrix_get(_position, 0, 2));
        }
    }
}

void pawan::__wake::updateBoundVorEffect(PawanRecvData pawanrecvdata,double &dt){
    int NbOfLfnLines = pawanrecvdata->NbOfLfnLines;
    int *NbOfAst = pawanrecvdata->NbOfAst;
    double *astpos = pawanrecvdata->astpos;
    double *circ = pawanrecvdata->circ;

    int astidx = 0;
    for(int i = 0; i<_numParticles; ++i) {
        gsl_vector *vb = gsl_vector_calloc(3);
        gsl_vector_const_view r = gsl_matrix_const_row(_position, i);
        for (size_t ilfn = 0; ilfn < NbOfLfnLines; ++ilfn) {
            for (size_t iast = 0; iast < NbOfAst[ilfn]-1; ++iast) {
                double t=0.0;
                gsl_vector *Si = gsl_vector_calloc(3);
                gsl_vector *rast1 = gsl_vector_calloc(3);
                gsl_vector *rast2 = gsl_vector_calloc(3);
                for(int j = 0; j<3; ++j) {
                    gsl_vector_set(rast1, j, astpos[3 * astidx + j]);
                    gsl_vector_set(rast2, j, astpos[3 * (astidx+1) + j]);
                }

                gsl_vector *rast12 = gsl_vector_calloc(3);
                gsl_vector *rast1r = gsl_vector_calloc(3);
                gsl_vector_memcpy(rast12,rast2);
                gsl_vector_sub(rast12,rast1);
                double span = gsl_blas_dnrm2(rast12);
                gsl_vector_memcpy(rast1r,&r.vector);
                gsl_vector_sub(rast1r,rast1);
                gsl_blas_ddot(rast1r,rast12,&t);
                t=t/(span*span);

                gsl_vector_memcpy(Si,rast12);
                gsl_blas_dscal(t,Si);
                gsl_vector_add(Si,rast1);

                gsl_vector *licap = gsl_vector_calloc(3);
                gsl_vector_memcpy(licap,rast12);
                gsl_blas_dscal(span,licap);

                gsl_vector *sicap = gsl_vector_calloc(3);
                gsl_vector_memcpy(sicap,&r.vector);
                gsl_vector_sub(sicap,Si);
                double simag = gsl_blas_dnrm2(sicap);
                gsl_blas_dscal(simag,sicap);

                gsl_vector *rast2r = gsl_vector_calloc(3);
                gsl_vector_memcpy(rast2r,&r.vector);
                gsl_vector_sub(rast2r,rast2);
                double alpha = 0.0;
                gsl_blas_ddot(rast2r,rast12,&alpha);
                alpha = acos(-alpha/gsl_blas_dnrm2(rast2r)/span);
                double beta = 0.0;
                gsl_blas_ddot(rast1r,rast12,&beta);
                beta = acos(-beta/gsl_blas_dnrm2(rast1r)/span);

                double circ1 = circ[astidx];
                double circ2 = circ[astidx+1];
                //t can be -ve
                double gamma0 = (circ1*fabs(span*(1+t)) + circ2*fabs(t*span))/span;
                double gamma1 = (circ2-circ1)/span;
                gsl_vector *vbi = gsl_vector_calloc(3);
                gsl_cross(licap,sicap,vbi);
                gsl_blas_dscal(0.25*M_1_PI/simag,vbi);
                double factor =gamma0*(cos(beta)+cos(alpha)) + simag*gamma1*(sin(beta)-sin(alpha));
                gsl_blas_dscal(factor,vbi);

                gsl_vector_add(vb,vbi);

                astidx++;
            }
        }

        for (int j = 0; j < 3; ++j) {
            gsl_matrix_set(_position, i, j,
                           gsl_vector_get(vb,j)*dt + gsl_matrix_get(_position, i, j));
        }
    }

/*
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
*/
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
    gsl_vector_fwrite(f,_volume);
    gsl_vector_fwrite(f,_birthstrength);
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

void pawan::__wake::translate(const double *x){
    for(size_t i = 0; i<_numParticles; ++i){
        for(size_t j = 0; j<3; ++j) {
            gsl_matrix_set(_position, i, j, x[j] + gsl_matrix_get(_position, i, j));
        }
    }
}

void pawan::__wake::rotate(const size_t &n,const double &a){
    for(size_t i = 0; i<_numParticles; ++i){
        if (n==0) { //rotation about x-axis
            gsl_matrix_set(_position, i, 1,
                           cos(a) * gsl_matrix_get(_position, i, 1) - sin(a) * gsl_matrix_get(_position, i, 2));
            gsl_matrix_set(_position, i, 2,
                           sin(a) * gsl_matrix_get(_position, i, 1) + cos(a) * gsl_matrix_get(_position, i, 2));
        }
        else if (n==1) { //rotation about y-axis
            gsl_matrix_set(_position, i, 0,
                           cos(a) * gsl_matrix_get(_position, i, 0) + sin(a) * gsl_matrix_get(_position, i, 2));
            gsl_matrix_set(_position, i, 2,
                           -sin(a) * gsl_matrix_get(_position, i, 0) + cos(a) * gsl_matrix_get(_position, i, 2));
        }
        else if (n==2) { //rotation about z-axis
            gsl_matrix_set(_position, i, 1,
                           cos(a) * gsl_matrix_get(_position, i, 0) - sin(a) * gsl_matrix_get(_position, i, 1));
            gsl_matrix_set(_position, i, 2,
                           sin(a) * gsl_matrix_get(_position, i, 0) + cos(a) * gsl_matrix_get(_position, i, 1));
        }
    }

}