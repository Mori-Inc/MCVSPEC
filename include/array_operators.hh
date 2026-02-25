#pragma once
#include <cmath>
#include <algorithm>
#include "integration.hh" // for State

template <size_t n_dim>
inline void add_scalar(const double a, const State<n_dim>& x, State<n_dim>& out){
    for(size_t i=0; i<n_dim; i++){
        out[i] = a + x[i];
    }
}

template <size_t n_dim>
inline void multiply_scalar(const double a, const State<n_dim>& x, State<n_dim>& out){
    for(size_t i=0; i<n_dim; i++){
        out[i] = a*x[i];
    }
}

template <size_t n_dim>
inline void multiply_scalar_inplace(State<n_dim>& x, const double a){
    for(size_t i=0; i<n_dim; i++){
        x[i] *= a;
    }
}

template <size_t n_dim>
inline void linear_combo(const double a, const State<n_dim>& x, const double b, const State<n_dim>& y, State<n_dim>& out){
    for(size_t i=0; i<n_dim; i++){
        out[i] = a*x[i] + b*y[i];
    }
}

template <size_t n_dim>
inline void multiply_elements(const State<n_dim>& x, const State<n_dim>& y, State<n_dim>& out){
    for(size_t i=0; i<n_dim; i++){
        out[i] = x[i]*y[i];
    }
}

template <size_t n_dim>
inline void divide_elements(const State<n_dim>& x, const State<n_dim>& y, State<n_dim>& out){
    for(size_t i=0; i<n_dim; i++){
        out[i] = x[i]/y[i];
    }
}

template <size_t n_dim>
inline void divide_elements_inplace(State<n_dim>& x, const State<n_dim>& y){
    for(size_t i=0; i<n_dim; i++){
        x[i] /= y[i];
    }
}

template <size_t n_dim>
inline void add_vector_inplace(State<n_dim>& x, const double a, const State<n_dim>& y){
    for(size_t i=0; i<n_dim; i++){
        x[i] += a*y[i];
    }
}

template <size_t n_dim>
inline void add_abs_vector_inplace(State<n_dim>& x, const double a, const State<n_dim>& y){
    for(size_t i=0; i<n_dim; i++){
        x[i] += a*std::abs(y[i]);
    }
}

template <size_t n_dim>
inline void element_max(const State<n_dim>& x, const State<n_dim>& y, State<n_dim>& out){
    for(size_t i=0; i<n_dim; i++){
        out[i] = std::max(std::abs(x[i]), std::abs(y[i]));
    }
}
