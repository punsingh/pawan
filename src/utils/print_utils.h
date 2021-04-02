/*! PArticle Wake ANalysis
 * @file print_utils.h
 * @brief Print utilities for PAWAN
 *
 * @author Puneet Singh
 * @date 04/01/2021
 *
 */

#ifndef PRINT_UTILS_H_
#define PRINT_UTILS_H_

#include <iostream>
#include <sstream>
#include <string>

/*! \fn inline void OUT(std::string s, std::ostream &os = std::cout)
 * \brief Print string
 * \param s	String
 * \param os	Output stream
 */
inline void OUT(std::string s, std::ostream &os = std::cout){
	os << s << std::endl;
};

/*! \fn inline void OUT(std::string s, T v, std::ostream &os = std::cout)
 * \brief Print string and value
 * \param s	String
 * \param v	Output quantity
 * \param os	Output stream
 */
//template <typename T>
//inline void OUT(std::string s, T v, std::ostream &os = std::cout){
	//os << s << " = "<< v << std::endl;
//};

#endif
