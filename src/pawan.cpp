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

int main(int argc, char* argv[]){

	std::cout << std::setprecision(16) << std::scientific;
	std::cout << "PArticle Wake ANalysis" << std::endl;
	pawan::__io *IO = new pawan::__io();
	IO->print();
	OUT(IO->getFile());
	pawan::__wake *W = new pawan::__wake(IO);
	W->print();
	delete W;
	delete IO;

	// End
	return EXIT_SUCCESS;
}
