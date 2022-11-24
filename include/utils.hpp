// include guard
#ifndef __utils_hpp__
#define __utils_hpp__

#include <vector>

// Return a string with a double in scientific notation
std::string scientific_format(double d, const int &width, const int &prec);

// Return a string with a vector<double> in scientific notation
std::string scientific_format(const std::vector<double> &v, const int &width, const int &prec);

#endif // __utils_hpp__