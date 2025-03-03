#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "tableau.hh"
#include <iostream>
#include <cmath>
#include <valarray>
#include <vector>

using std::function, std::valarray, std::vector;
using std::find, std::fill, std::begin, std::end;
using std::abs, std::max, std::min;
using std::pow;
using std::cout, std::endl;

double norm(valarray<double>);
valarray<double> element_max(valarray<double>, valarray<double>);
void Dormand_Prince(function<valarray<double>(double, valarray<double>, void*)>, void*, vector<double>*, vector<valarray<double>>*, double, vector<double>*, vector<valarray<double>>*, double, double, int);

#endif
