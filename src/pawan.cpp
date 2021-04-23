/*! Particle Wake Analysis
 * \file orca.cpp
 * \brief Main executable code
 * @author Puneet Singh
 * @date 03/28/2021
 */

#include <iostream>
#include <iomanip> // Required for set precision

#include "utils/print_utils.h"
#include "io/io.h"
#include "wake/wake.h"
#include "src/integration/integration.h"

int main(int argc, char* argv[]){

	std::cout << std::setprecision(16) << std::scientific;
	std::cout << "PArticle Wake ANalysis" << std::endl;
	pawan::__io *IO = new pawan::__io();
	pawan::__wake *W = new pawan::__wake(1.0,1.0,1.0,16);
	pawan::__integration *IN = new pawan::__integration(1.0,5);
	IN->integrate(W,IO);
	//W->calculateInteraction();
	//W->print();
	//W->write(IO);
	delete IN;
	delete W;
	delete IO;

	// End
	return EXIT_SUCCESS;
}
