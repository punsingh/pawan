/*! PArticle Wake ANalysis
 * \file io.cpp
 * \brief Routines for Input Output operations
 *
 * @author Puneet Singh
 * @date 04/01/2021
 */
#include "io.h"
#include <unistd.h>

pawan::__io::__io(){
	_root = "data/";
	_name = "temp";
}

void pawan::__io::print(){
	OUT(_root);
	OUT(_name);
}

std::string pawan::__io::getFile(){
	return (_root + _name);
}

FILE* pawan::__io::create_binary_file(std::string suffix = ".bin"){
	std::string filename = _root + _name + suffix;
	FILE *f = fopen(filename.c_str(),"wb");
	return f;
}

FILE* pawan::__io::open_binary_file(std::string suffix = ".bin"){
	std::string filename = _root + _name + suffix;
	FILE *f = fopen(filename.c_str(),"rb");
	return f;
}

