/*! PArticle Wake ANalysis
 * \file parallel.cpp
 * \brief Routines for Open MP parallilized interactions
 *
 * @author Puneet Singh
 * @date 04/24/2021
 */
#include "parallel.h"
//#include "testcuda1.h"



pawan::__parallel::__parallel(__wake *W):__interaction(W){}
pawan::__parallel::__parallel(__wake *W1, __wake *W2):__interaction(W1,W2){}

void pawan::__parallel::interact(__wake *W){
    //testcuda_call();

    double** Wpos_arr = la_gslalloc(W->_position,W->_numParticles,W->_numDimensions);
    double** Wvor_arr = la_gslalloc(W->_vorticity,W->_numParticles,W->_numDimensions);
    double** Wretvor_arr = la_gslalloc(W->_retvorcity,W->_numParticles,W->_numDimensions);
    double** Wvel_arr = la_gslalloc(W->_velocity,W->_numParticles,W->_numDimensions);
    double* Wrad_vec = la_gslalloc(W->_radius,W->_numParticles);
    double* Wvol_vec = la_gslalloc(W->_volume,W->_numParticles);
    la_print(Wpos_arr, W->_numParticles, W->_numDimensions);
    /*
    la_print(Wpos_arr, W->_numParticles, W->_numDimensions);
    la_print(Wretvor_arr, W->_numParticles, W->_numDimensions);
    la_print(Wrad_vec, W->_numParticles);
    la_print(Wvol_vec, W->_numParticles);

    printf("printing one way---------------------------------- \n");
    printf("%f \n",*(*(Wpos_arr)));
    printf("%f \n",*(*(Wpos_arr+1)));
    printf("%f \n",*(*(Wpos_arr+2)));

    printf("printing another way---------------------------------- \n");
    printf("%f \n",*(*Wpos_arr));
    printf("%f \n",*(*Wpos_arr+1));
    printf("%f \n",*(*Wpos_arr+2));
*/
     for(size_t i_src = 0; i_src < W->_numParticles; ++i_src){
/*
//        const double *r_src = la_gsl_matrix_row(Wpos_arr,i_src,W->_numDimensions);//make the return itself const??
//        const double *a_src = la_gsl_matrix_row(Wvor_arr,i_src,W->_numDimensions);//make the return itself const??
         const double *r_src = la_gsl_matrix_row(Wpos_arr,i_src);
         const double *a_src = la_gsl_matrix_row(Wvor_arr,i_src);
        double s_src = la_gsl_vector_get(Wrad_vec,i_src);// these could be const as well ??
        double v_src = la_gsl_vector_get(Wvol_vec,i_src);// these could be const as well ??
        //set_gsl_drow(r_src_vec, Wpos_arr, i_src, W->_numParticles);
        //la_printVec(r_src_vec, W->_numParticles);
*/

		gsl_vector_const_view r_src = gsl_matrix_const_row(W->_position,i_src);
		gsl_vector_const_view a_src = gsl_matrix_const_row(W->_vorticity,i_src);
		double s_src = gsl_vector_get(W->_radius,i_src);
		double v_src = gsl_vector_get(W->_volume,i_src);
		double vx = 0.0, vy = 0.0, vz = 0.0;
		double qx = 0.0, qy = 0.0, qz = 0.0;
		#pragma omp parallel for reduction(+:vx,vy,vz,qx,qy,qz)
		for(size_t i_trg = i_src + 1; i_trg < W->_numParticles; ++i_trg){
/*            const double *r_trg = la_gsl_matrix_row(Wpos_arr,i_trg); //make the return itself const??
            const double *a_trg = la_gsl_matrix_row(Wvor_arr,i_trg); //make the return itself const??
            double *dr_trg = la_gsl_matrix_row(Wvel_arr,i_trg);
            double *da_trg= la_gsl_matrix_row(Wretvor_arr,i_trg);
            double s_trg = la_gsl_vector_get(Wrad_vec,i_trg);// these could be const as well ??
            double v_trg = la_gsl_vector_get(Wvol_vec,i_trg);// these could be const as well ??
*/			gsl_vector_const_view r_trg = gsl_matrix_const_row(W->_position,i_trg);
			gsl_vector_const_view a_trg= gsl_matrix_const_row(W->_vorticity,i_trg);
			gsl_vector_view dr_trg = gsl_matrix_row(W->_velocity,i_trg);
			gsl_vector_view da_trg = gsl_matrix_row(W->_retvorcity,i_trg);
			double s_trg = gsl_vector_get(W->_radius,i_trg);
			double v_trg = gsl_vector_get(W->_volume,i_trg);
			double vx_s = 0.0, vy_s = 0.0, vz_s = 0.0;
			double qx_s = 0.0, qy_s = 0.0, qz_s = 0.0;
			INTERACT(_nu,s_src,s_trg,&r_src.vector,&r_trg.vector,&a_src.vector,&a_trg.vector,v_src,v_trg,&dr_trg.vector,&da_trg.vector,vx_s,vy_s,vz_s,qx_s,qy_s,qz_s);
 //           INTERACT(_nu,s_src,s_trg,r_src,r_trg,a_src,a_trg,v_src,v_trg,dr_trg,da_trg,vx_s,vy_s,vz_s,qx_s,qy_s,qz_s);

            vx += vx_s;
			vy += vy_s;
			vz += vz_s;
			qx += qx_s;
			qy += qy_s;
			qz += qz_s;
		}
        //std::cout << omp_get_num_threads() << std::endl;
/*
        double *dr_src = la_gsl_matrix_row(Wvel_arr,i_src);
		la_gsl_vector_set(dr_src,0,vx + la_gsl_vector_get(dr_src,0));
         la_gsl_vector_set(dr_src,1,vy + la_gsl_vector_get(dr_src,1));
         la_gsl_vector_set(dr_src,2,vz + la_gsl_vector_get(dr_src,2));
         double *da_src = la_gsl_matrix_row(Wretvor_arr,i_src);
         la_gsl_vector_set(da_src,0,qx + la_gsl_vector_get(da_src,0));
         la_gsl_vector_set(da_src,1,qy + la_gsl_vector_get(da_src,1));
         la_gsl_vector_set(da_src,2,qz + la_gsl_vector_get(da_src,2));
         */
          gsl_vector_view dr_src = gsl_matrix_row(W->_velocity,i_src);
          gsl_vector_set(&dr_src.vector,0,vx + gsl_vector_get(&dr_src.vector,0));
          gsl_vector_set(&dr_src.vector,1,vy + gsl_vector_get(&dr_src.vector,1));
          gsl_vector_set(&dr_src.vector,2,vz + gsl_vector_get(&dr_src.vector,2));
          gsl_vector_view da_src = gsl_matrix_row(W->_retvorcity,i_src);
          gsl_vector_set(&da_src.vector,0,qx + gsl_vector_get(&da_src.vector,0));
          gsl_vector_set(&da_src.vector,1,qy + gsl_vector_get(&da_src.vector,1));
          gsl_vector_set(&da_src.vector,2,qz + gsl_vector_get(&da_src.vector,2));
 	}

    //la_gslalloc_reverse(Wpos_arr,W->_position,W->_numParticles,W->_numDimensions);
    //la_gslalloc_reverse(Wretvor_arr,W->_retvorcity,W->_numParticles,W->_numDimensions);

    la_dealloc(Wpos_arr,W->_numParticles);
    la_dealloc(Wvor_arr,W->_numParticles);
    la_dealloc(Wretvor_arr,W->_numParticles);
    la_dealloc(Wvel_arr,W->_numParticles);
    la_dealloc(Wrad_vec);
    la_dealloc(Wvol_vec);
}

/*!
 * void pawan::__parallel::interact(__wake *W){
    for(size_t i_src = 0; i_src < W->_numParticles; ++i_src){
        gsl_vector_const_view r_src = gsl_matrix_const_row(W->_position,i_src);
        gsl_vector_const_view a_src = gsl_matrix_const_row(W->_vorticity,i_src);
        double s_src = gsl_vector_get(W->_radius,i_src);
        double v_src = gsl_vector_get(W->_volume,i_src);
        double vx = 0.0, vy = 0.0, vz = 0.0;
        double qx = 0.0, qy = 0.0, qz = 0.0;
#pragma omp parallel for reduction(+:vx,vy,vz,qx,qy,qz)
        for(size_t i_trg = i_src + 1; i_trg < W->_numParticles; ++i_trg){
            gsl_vector_const_view r_trg = gsl_matrix_const_row(W->_position,i_trg);
            gsl_vector_const_view a_trg= gsl_matrix_const_row(W->_vorticity,i_trg);
            gsl_vector_view dr_trg = gsl_matrix_row(W->_velocity,i_trg);
            gsl_vector_view da_trg = gsl_matrix_row(W->_retvorcity,i_trg);
            double s_trg = gsl_vector_get(W->_radius,i_trg);
            double v_trg = gsl_vector_get(W->_volume,i_trg);
            double vx_s = 0.0, vy_s = 0.0, vz_s = 0.0;
            double qx_s = 0.0, qy_s = 0.0, qz_s = 0.0;
            INTERACT(_nu,s_src,s_trg,&r_src.vector,&r_trg.vector,&a_src.vector,&a_trg.vector,v_src,v_trg,&dr_trg.vector,&da_trg.vector,vx_s,vy_s,vz_s,qx_s,qy_s,qz_s);
            vx += vx_s;
            vy += vy_s;
            vz += vz_s;
            qx += qx_s;
            qy += qy_s;
            qz += qz_s;
        }
        //std::cout << omp_get_num_threads() << std::endl;
        gsl_vector_view dr_src = gsl_matrix_row(W->_velocity,i_src);
        gsl_vector_set(&dr_src.vector,0,vx + gsl_vector_get(&dr_src.vector,0));
        gsl_vector_set(&dr_src.vector,1,vy + gsl_vector_get(&dr_src.vector,1));
        gsl_vector_set(&dr_src.vector,2,vz + gsl_vector_get(&dr_src.vector,2));
        gsl_vector_view da_src = gsl_matrix_row(W->_retvorcity,i_src);
        gsl_vector_set(&da_src.vector,0,qx + gsl_vector_get(&da_src.vector,0));
        gsl_vector_set(&da_src.vector,1,qy + gsl_vector_get(&da_src.vector,1));
        gsl_vector_set(&da_src.vector,2,qz + gsl_vector_get(&da_src.vector,2));
    }
}
*/

void pawan::__parallel::interact(__wake *W1, __wake *W2){
/*
	for(size_t i_src = 0; i_src < W1->_numParticles; ++i_src){
		gsl_vector_const_view r_src = gsl_matrix_const_row(W1->_position,i_src);
		gsl_vector_const_view a_src = gsl_matrix_const_row(W1->_vorticity,i_src);
		double s_src = gsl_vector_get(W1->_radius,i_src);
		double v_src = gsl_vector_get(W1->_volume,i_src);
		double vx = 0.0, vy = 0.0, vz = 0.0;
		double qx = 0.0, qy = 0.0, qz = 0.0;
		#pragma omp parallel for reduction(+:vx,vy,vz,qx,qy,qz)
		for(size_t i_trg = 0; i_trg < W2->_numParticles; ++i_trg){
			gsl_vector_const_view r_trg = gsl_matrix_const_row(W2->_position,i_trg);
			gsl_vector_const_view a_trg = gsl_matrix_const_row(W2->_vorticity,i_trg);
			gsl_vector_view dr_trg = gsl_matrix_row(W2->_velocity,i_trg);
			gsl_vector_view da_trg = gsl_matrix_row(W2->_retvorcity,i_trg);
			double s_trg = gsl_vector_get(W2->_radius,i_trg);
			double v_trg = gsl_vector_get(W2->_volume,i_trg);
			double vx_s = 0.0, vy_s = 0.0, vz_s = 0.0;
			double qx_s = 0.0, qy_s = 0.0, qz_s = 0.0;
			INTERACT(_nu,s_src,s_trg,&r_src.vector,&r_trg.vector,&a_src.vector,&a_trg.vector,v_src,v_trg,&dr_trg.vector,&da_trg.vector,vx_s,vy_s,vz_s,qx_s,qy_s,qz_s);
			vx += vx_s;
			vy += vy_s;
			vz += vz_s;
			qx += qx_s;
			qy += qy_s;
			qz += qz_s;
		}
		gsl_vector_view dr_src = gsl_matrix_row(W1->_velocity,i_src);
		gsl_vector_view da_src = gsl_matrix_row(W1->_retvorcity,i_src);
		gsl_vector_set(&dr_src.vector,0,vx + gsl_vector_get(&dr_src.vector,0));
		gsl_vector_set(&dr_src.vector,1,vy + gsl_vector_get(&dr_src.vector,1));
		gsl_vector_set(&dr_src.vector,2,vz + gsl_vector_get(&dr_src.vector,2));
		gsl_vector_set(&da_src.vector,0,qx + gsl_vector_get(&da_src.vector,0));
		gsl_vector_set(&da_src.vector,1,qy + gsl_vector_get(&da_src.vector,1));
		gsl_vector_set(&da_src.vector,2,qz + gsl_vector_get(&da_src.vector,2));
	}
*/
 }
