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
#include "src/networkinterface/networkdatastructures.h"
#include "src/networkinterface/networkinterface.h"
#include "src/networkinterface/networkinterface.cpp" //templates included this way

#define OUTPUTIP "127.0.0.1"
#define NETWORKBUFFERSIZE 50
#define PORT 8899


int main(int argc, char* argv[]){

    NetworkInterfaceTCP<OPawanRecvData,OPawanSendData>
            networkCommunicatorTest(PORT, OUTPUTIP, PORT, NETWORKBUFFERSIZE, true);
    networkCommunicatorTest.socket_init();
    OPawanRecvData opawanrecvdata;
    networkCommunicatorTest.recieve_data(opawanrecvdata);

    std::cout << std::setprecision(16) << std::scientific;
    PAWAN();
    pawan::__io *IO = new pawan::__io();
    PawanRecvData pawanrecvdata = &opawanrecvdata;
    pawan::__wake *W = new pawan::__wake(pawanrecvdata);
    pawan::__interaction *S = new pawan::__interaction(W);
    pawan::__integration *IN = new pawan::__integration(1.568659,100);
//    pawan::__integration *IN = new pawan::__rk4(0.1568659,10);
    IN->integrate(S,IO,&networkCommunicatorTest);

    delete IN;
    delete S;
    delete IO;
    // End
    printf("---------------------+++++++++++++++++++++++++!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EVERYTHING FINISHED SUCCESFULLY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!+++++++++++++++++++++++++++-------------------------------------\n");
    //std::cout << FLT_DIG << DBL_DIG << std::endl;
    return EXIT_SUCCESS;
}