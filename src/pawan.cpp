/*! Particle Wake Analysis
 * \file orca.cpp
 * \brief Main executable code
 * @author Puneet Singh
 * @date 03/28/2021
 */

#include <iostream>
#include <stdio.h>

#include "wake/wake.h"

int main(int argc, char* argv[]){

	std::cout << "PArticle Wake ANalysis" << std::endl;
	__wake *W = new __wake();
	W->print();
	delete W;

	// End
	return EXIT_SUCCESS;
}
