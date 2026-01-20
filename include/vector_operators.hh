#pragma once
#include <vector>
#include <cmath>

using std::vector;
using std::max;
using std::abs;

inline void add_scalar(const double a, const vector<double>& x, vector<double>& out){
    for(int i=0; i<x.size(); i++){
        out[i] = a + x[i];
    }
}

inline void multiply_scalar(const double a, const vector<double>& x, vector<double>& out){
    for(int i=0; i<x.size(); i++){
        out[i] = a*x[i];
    }
}

inline void multiply_scalar_inplace(vector<double>& x, const double a){
    for(int i=0; i<x.size(); i++){
        x[i] *= a;
    }
}

inline void linear_combo(const double a, const vector<double>& x, const double b, const vector<double>& y, vector<double>& out){
    for(int i=0; i<x.size(); i++){
        out[i] = a*x[i] + b*y[i];
    }
}

inline void multiply_elements(const vector<double>& x, const vector<double>& y, vector<double>& out){
    for(int i=0; i<x.size(); i++){
        out[i] = x[i]*y[i];
    }
}

inline void divide_elements(const vector<double>& x, const vector<double>& y, vector<double>& out){
    for(int i=0; i<x.size(); i++){
        out[i] = x[i]/y[i];
    }
}

inline void divide_elements_inplace(vector<double>& x, const vector<double>& y){
    for(int i=0; i<x.size(); i++){
        x[i] /= y[i];
    }
}

inline void add_vector_inplace(vector<double>& x, const double a, const vector<double>& y){
    for(int i=0; i<x.size(); i++){
        x[i] += a*y[i];
    }
}

inline void add_abs_vector_inplace(vector<double>& x, const double a, const vector<double>& y){
    for(int i=0; i<x.size(); i++){
        x[i] += a*std::abs(y[i]);
    }
}

inline void element_max(const vector<double>& x, const vector<double>& y, vector<double>& out){
    for(int i=0; i<x.size(); i++){
        out[i] = max(std::abs(x[i]), std::abs(y[i]));
    }
}
