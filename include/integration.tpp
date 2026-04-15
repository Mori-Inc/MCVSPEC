#pragma once
#include "integration.hh"
#include "tableau.hh"
#include "array_operators.hh"
#include <algorithm>
#include <numeric>
#include <array>
#include <cmath>

template <size_t n_dim>
inline double norm(const State<n_dim>& x){
    double norm = 0.;
    for(size_t i=0; i<n_dim; i++){
        norm += x[i]*x[i];
    }
    return sqrt(norm/double(n_dim));
}

template <size_t n_dim, class RHS>
Integrator<n_dim,RHS>::Integrator(RHS rhs, const double absolute_err, const double relative_err):
    func(std::move(rhs)), abs_err(absolute_err), rel_err(relative_err)
{}

template <size_t n_dim, class RHS>
void Integrator<n_dim,RHS>::Set_Initial_Step(const double& t0, const State<n_dim>& y0){
    tol.fill(abs_err);
    add_abs_vector_inplace(tol, rel_err, y0);
    double h_0 = 1e-2*norm(y0)/norm(k[0]);
    linear_combo(1.,y0, dir*h_0, k[0],y_internal);
    t_internal = t0+dir*h_0;
    func(t_internal, y_internal, k[1]);
    linear_combo(1., k[1], -1., k[0], buffer);
    divide_elements_inplace(buffer, tol);
    double delta = norm(buffer)/h_0;
    divide_elements(k[0], tol, buffer);
    double h_1 = std::pow(1e-2/std::max(delta,norm(buffer)),1./tableau::order);
    h = std::min(1e2*h_0,h_1);
}

template <size_t n_dim, class RHS>
void Integrator<n_dim,RHS>::Initialize(double& t, const double t_end, const State<n_dim>& y){
    dir = (0. < (t_end-t)) - ((t_end-t) < 0.);
    func(t, y, k[0]);
    Set_Initial_Step(t, y);
}

template <size_t n_dim, class RHS>
void Integrator<n_dim,RHS>::Prepare_Step(const double& t,const State<n_dim>& y, double& h_new){
    using tableau::a;
    using tableau::b;
    using tableau::c;
    using tableau::e;
    bool step_succeded = false;
    bool step_failed = false;
    double err_norm;
    while(!step_succeded){
        multiply_scalar(e[0], k[0], err_arr);
        for(size_t i=1; i<tableau::n_stages; i++){
            dy.fill(0.0);
            for(size_t j = 0; j<i; j++){
                add_vector_inplace(dy, a[i][j], k[j]);
            }
            t_internal = t+c[i]*dir*h;
            linear_combo(1.,y,dir*h,dy,y_internal);
            func(t_internal, y_internal, k[i]);
            add_vector_inplace(err_arr, e[i],k[i]);
        }
        dy.fill(0.0);
        for(size_t i=0; i<tableau::n_stages; i++){
            add_vector_inplace(dy, b[i],k[i]);
        }
        multiply_scalar_inplace(dy, dir*h);
        dt = dir*h;
        linear_combo(1.,y,1.,dy,y_internal);
        func(t+dt, y_internal, k[tableau::n_stages]);
        add_vector_inplace(err_arr, e[tableau::n_stages],k[tableau::n_stages]);
        double sqr_err = 0;
        for(size_t i=0; i<y.size(); i++){
            const double tol_i = abs_err + rel_err*std::max(std::abs(y[i]), std::abs(y_internal[i]));
            const double err_i = h*err_arr[i]/tol_i;
            sqr_err += err_i*err_i;
        }
        err_norm = sqrt(sqr_err/n_dim);
        if(err_norm < 1.){
            step_succeded = true;
            if (err_norm == 0){
                h_new = 5.*h;
            }
            else if(step_failed){
                h_new = std::min(0.9*std::pow(err_norm,-1./tableau::order), 1.)*h;
            }
            else{
                h_new = std::min(0.9*std::pow(err_norm,-1./tableau::order), 5.)*h;
            }
        }
        else if(err_norm >= 1.){
            step_failed = true;
            h *= std::max(0.9*std::pow(err_norm,-1./tableau::order), 0.2);
        }
        else{
            // if err is nan
            step_failed = true;
            h *= 0.2;
        }
    }
}

template <size_t n_dim, class RHS>
void Integrator<n_dim,RHS>::Step(double& t, State<n_dim>& y){
    double h_new = h;
    dt = h;
    Prepare_Step(t, y, h_new);
    add_vector_inplace(y, 1., dy);
    t += dt;
    h = h_new;
    swap(k[0], k[tableau::n_stages]);
}

template <size_t n_dim, class RHS>
void Integrator<n_dim,RHS>::Dense_Step(double& t, State<n_dim>& y){
    using tableau::p;
    double h_new=h;
    dt = h;
    t_old = t;
    q[0]=y;
    Prepare_Step(t, y, h_new);
    for(size_t i = 1; i<tableau::order; i++){
        q[i].fill(0.0);
        for(size_t j = 0; j<tableau::n_stages+1; j++){
            add_vector_inplace(q[i], p[j][i-1], k[j]);
        }
        multiply_scalar_inplace(q[i], std::pow(h,1-static_cast<int>(i)));
    }
    add_vector_inplace(y, 1., dy);
    t += dt;
    h = h_new;
    swap(k[0], k[tableau::n_stages]);
}

template <size_t n_dim, class RHS>
void Integrator<n_dim,RHS>::Interpolate(double t, State<n_dim>& y){
    double dt = dir*(t-t_old);
    y = q[0];
    for(size_t i=1; i<tableau::order; i++){
        add_vector_inplace(y, dir*std::pow(dt,i), q[i]);
    }
}
