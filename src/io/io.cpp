/*! PArticle Wake ANalysis
 * \file io.cpp
 * \brief Routines for Input Output operations
 *
 * @author Puneet Singh
 * @date 04/01/2021
 */
#include "io.h"

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

