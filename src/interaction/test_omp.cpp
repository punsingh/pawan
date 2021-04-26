#include <iostream>
#include <iomanip> // Required for set precision
#include <sstream>
#include <string>
int main(int argc, char* argv[]){
	#pragma omp parallel
	{
		std::cout << "Thread" << std::endl;
	}
	// End
	return EXIT_SUCCESS;
}
