/*! Particle Wake Analysis
 * \file pawan.cpp
 * \brief Main executable code
 * @author Puneet Singh
 * @date 03/28/2021
 */

#include <iostream>
#include <iomanip> // Required for set precision

#include "utils/print_utils.h"
#include "io/io.h"
#include "wake/wake.h"
#include "wake/ring.h"
#include "wake/square.h"
#include "wake/vring.h"
#include "src/interaction/interaction.h"
#include "src/interaction/parallel.h"
#include "src/integration/integration.h"
#include "src/resolve/resolve.h"
#include "src/integration/rk4.h"

int main(int argc, char* argv[]){

	std::cout << std::setprecision(3) << std::scientific;
	PAWAN();	// Print PAWAN
	// Ideal vortex ring validation
	pawan::__io *IO = new pawan::__io();
	pawan::__wake *W = new pawan::__vring(1.0,0.1,2);
	pawan::__interaction *S = new pawan::__interaction(W);
	//pawan::__interaction *S = new pawan::__parallel(W);
	pawan::__resolve *R = new pawan::__resolve();
	S->diagnose();
	R->rebuild(S,IO);
	W->print();
	S->diagnose();
	S->solve();
	W->print();
	//pawan::__wake *W = new pawan::__wake(1.0,1.0,1.0,1024);
	//pawan::__wake *W = new pawan::__ring(1.0,1.0,1.0,32);
	//pawan::__wake *W = new pawan::__ring(0.4,1.0,0.2,32,3);
	//pawan::__wake *W = new pawan::__ring(0.4,1.0,0.2,32,3);
	//pawan::__wake *W = new pawan::__square(3.,2.0,0.1,21);
	//pawan::__wake *W2 = new pawan::__ring(0.4,1.0,0.2,32,3,true);
	//FILE *f = IO->create_binary_file(".wake");
	//double t = 0.0;
	//fwrite(&t,sizeof(double),1,f);	
	////pawan::__wake *W = new pawan::__vring(1.0,0.1,32,3);
	//size_t nWake = 1;
	//fwrite(&nWake,sizeof(size_t),1,f);	
	//W->write(f);
	//fclose(f);
	//W->translate(2,-1);
	//W2->translate(2,-0.4);
	//W->print();
	//pawan::__interaction *S = new pawan::__interaction(W);
	//pawan::__interaction *S = new pawan::__parallel(W);
	//pawan::__interaction *S = new pawan::__interaction(W,W2);
	//pawan::__interaction *S = new pawan::__parallel(W,W2);
	//pawan::__integration *IN = new pawan::__integration(0.02,10);
	//pawan::__integration *IN = new pawan::__integration(0.025,1);
	//pawan::__integration *IN = new pawan::__rk4(4,64);
	//pawan::__integration *IN = new pawan::__rk4(0.1,2);
	//IN->integrate(S,IO);
	//delete IN;
	delete R;
	delete S;
	delete W;
	delete IO;

	// End
	return EXIT_SUCCESS;
}
