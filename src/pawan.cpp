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
#include "src/interaction/interaction.h"
#include "src/interaction/parallel.h"
#include "src/integration/integration.h"

int main(int argc, char* argv[]){

	std::cout << std::setprecision(16) << std::scientific;
	std::cout << "PArticle Wake ANalysis" << std::endl;
	pawan::__io *IO = new pawan::__io();
	//pawan::__wake *W = new pawan::__wake(1.0,1.0,1.0,1024);
	pawan::__wake *W = new pawan::__ring(1.0,1.0,1.0,32);
	//pawan::__wake *W = new pawan::__ring(0.1,1.0,0.2,32,3);
	//W->print();
	//pawan::__interaction *S = new pawan::__interaction(W);
	pawan::__interaction *S = new pawan::__parallel(W);
	pawan::__integration *IN = new pawan::__integration(1,16);
	IN->integrate(S,IO);
	delete IN;
	delete S;
	delete W;
	delete IO;

	// End
	return EXIT_SUCCESS;
}
