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
#include "src/interaction/interaction.h"
#include "src/interaction/parallel.h"
#include "src/integration/integration.h"
#include "src/integration/rk4.h"

int main(int argc, char* argv[]){
	std::cout << std::setprecision(16) << std::scientific;
	PAWAN();
	pawan::__io *IO = new pawan::__io();
//	pawan::__wake *W = new pawan::__wake(1.0,1.0,1.0,1024);
//	pawan::__wake *W = new pawan::__ring(1.0,1.0,1.0,32);
//	pawan::__wake *W = new pawan::__ring(0.4,1.0,0.2,32,1);
	pawan::__wake *W = new pawan::__square(3.,20.0,0.1,2);
	//pawan::__wake *W2 = new pawan::__ring(0.4,1.0,0.2,32,3,true);
    std::cout << "pawan.cpp------------translating created object" << std::endl;
    W->translate(2,-1);
    //	W2->translate(2,-0.4);
	//W->print();
	//pawan::__interaction *S = new pawan::__interaction(W);
	pawan::__interaction *S = new pawan::__parallel(W);          //is there a specific reason this a pointer of base class???
	//pawan::__interaction *S = new pawan::__interaction(W,W2);
	//pawan::__interaction *S = new pawan::__parallel(W,W2);
	//pawan::__integration *IN = new pawan::__integration(2,64);
	pawan::__integration *IN = new pawan::__rk4(0.1,1);
	IN->integrate(S,IO);
	delete IN;
	delete S;
	delete W;
	delete IO;

	// End
	printf("---------------------+++++++++++++++++++++++++!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EVERYTHING FINISHED SUCCESFULLY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!+++++++++++++++++++++++++++-------------------------------------");
	return EXIT_SUCCESS;
}
