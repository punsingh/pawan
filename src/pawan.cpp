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
//#include <float.h>
int main(int argc, char* argv[]){
	std::cout << std::setprecision(16) << std::scientific;
	PAWAN();
	pawan::__io *IO = new pawan::__io();
    pawan::__wake *W = new pawan::__square(2.,5.0,0.1,80);
    //pawan::__wake *W = new pawan::__ring(1.0,5.0,0.1,80);
    W->translate(2,-1);
/*    pawan::__interaction *S = new pawan::__parallel(W);
    pawan::__integration *IN = new pawan::__rk4(5,100);
    IN->integrate(S,IO);
*/
    //for cuda parallelisation
    pawan::__interaction *S = new pawan::__interaction(W);
    pawan::__integration *IN = new pawan::__integration(5.0,100);
    IN->integrate_cuda(S,IO,W);

    delete IN;
    delete S;
    delete W;
    delete IO;
	// End
	printf("---------------------+++++++++++++++++++++++++!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EVERYTHING FINISHED SUCCESFULLY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!+++++++++++++++++++++++++++-------------------------------------\n");
    //std::cout << FLT_DIG << DBL_DIG << std::endl;
    return EXIT_SUCCESS;
}
