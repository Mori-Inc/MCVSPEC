#pragma once
#include "integration.hh"
#include <algorithm>
#include <numeric>
#include <vector>
using std::fill;
using std::begin;
using std::end;
using std::max;
using std::min;
using std::pow;
using std::inner_product;

template <typename model>
using exec = void (model::*)(double, const vector<double>&, vector<double>&) const;

inline double norm(const vector<double>& x){
    return sqrt(inner_product(x.begin(), x.end(), x.begin(), 0)/x.size());
}

template <typename model, exec<model> call>
Integrator<model,call>::Integrator(const model& function, const int n_dim, const double absolute_err, const double relative_err):
    func(&function), abs_err(absolute_err), rel_err(relative_err), buffer(n_dim), y_internal(n_dim), t_internal(0), tol(n_dim, abs_err)
{}

template <typename model, exec<model> call>
void Integrator<model,call>::Set_Initial_Step(const double& t0, const vector<double>& y0){
    add_abs_vector_inplace(tol, rel_err, y0);
    double h_0 = 1e-2*norm(y0)/norm(k[0]);
    linear_combo(1.,y0, dir*h_0, k[0],y_internal);
    t_internal = t0+dir*h_0;
    (func->*call)(t_internal, y_internal, k[1]);
    linear_combo(1., k[1], -1., k[2], buffer);
    divide_elements_inplace(buffer, tol);
    double delta = norm(buffer)/h_0;
    divide_elements(k[0], tol, buffer);
    double h_1 = pow(1e-2/max(delta,norm(buffer)),1./order);
    h = min(1e2*h_0,h_1);
}

template <typename model, exec<model> call>
void Integrator<model,call>::Initialize(double& t, const double t_end, const vector<double>& y){
    fill(begin(k), end(k), vector<double>(0.,y.size()));
    fill(begin(q), end(q), vector<double>(0.,y.size()));
    dir = (0. < (t_end-t)) - ((t_end-t) < 0.);
    (func->*call)(t, y, k[0]);
    Set_Initial_Step(t, y);
}

template <typename model, exec<model> call>
int Integrator<model,call>::Integrate(double& t, const double t_end, vector<double>& y){
    Initialize(t, t_end, y);
    int n_steps=0;

    while(dir*(t_end-t) > 0 && n_steps < max_itter){
        h = min(h,dir*(t_end-t));
        Step(t, y);
        n_steps++;
    }
    if(n_steps>=max_itter){
        return 0;
    }
    return 1;
}

template <typename model, exec<model> call>
void Integrator<model,call>::Prepare_Step(const double& t, double& dt,const vector<double>& y, vector<double>& dy, double& h_new){
    using tableau::a;
    using tableau::b;
    using tableau::c;
    using tableau::e;
    bool step_succeded = false;
    bool step_failed = false;
    vector<double> err_arr(y.size());
    double err_norm;
    while(!step_succeded){
        multiply_scalar(e[0], k[0], err_arr);
        for(uint i=1; i<n_stages; i++){
            fill(dy.begin(), dy.end(), 0.0);
            for(int j = 0; j<i; j++){
                add_vector_inplace(dy, a[i][j], k[j]);
            }
            t_internal = t+c[i]*dir*h;
            linear_combo(1.,y,dir*h,dy,y_internal);
            (func->*call)(t_internal, y_internal, k[i]);
            add_vector_inplace(err_arr, e[i],k[i]);
        }
        fill(dy.begin(), dy.end(), 0.0);
        for(uint i=0; i<n_stages; i++){
            add_vector_inplace(dy, b[i],k[i]);
        }
        multiply_scalar_inplace(dy, dir*h);
        dt = dir*h;
        linear_combo(1.,y,1.,dy,y_internal);
        (func->*call)(t+dt, y_internal, k[n_stages]);
        fill(tol.begin(), tol.end(), abs_err);
        element_max(y, y_internal,y_internal);
        add_abs_vector_inplace(tol, rel_err, y_internal);
        add_vector_inplace(err_arr, h*e[n_stages],k[n_stages]);
        divide_elements_inplace(err_arr,tol);
        err_norm = norm(err_arr);
        if(err_norm < 1.){
            step_succeded = true;
            if (err_norm == 0){
                h_new = 5.*h;
            }
            else if(step_failed){
                h_new = min(0.9*pow(err_norm,-1./order), 1.)*h;
            }
            else{
                h_new = min(0.9*pow(err_norm,-1./order), 5.)*h;
            }
        }
        else if(err_norm >= 1.){
            step_failed = true;
            h *= max(0.9*pow(err_norm,-1./order), 0.2);
        }
        else{
            // if err is nan
            step_failed = true;
            h *= 0.2;
        }
    }
}

template <typename model, exec<model> call>
void Integrator<model,call>::Step(double& t, vector<double>& y){
    double dt=h, h_new=h;
    vector<double> dy(y.size());
    Prepare_Step(t, dt, y, dy, h_new);
    add_vector_inplace(y, 1., dy);
    t += dt;
    h = h_new;
    k[0] = k[n_stages];
}

template <typename model, exec<model> call>
void Integrator<model,call>::Dense_Step(double& t, vector<double>& y){
    using tableau::p;
    double dt=h, h_new=h;
    vector<double> dy(y.size());
    t_old = t;
    q[0]=y;
    Prepare_Step(t, dt, y, dy, h_new);
    for(int i = 1; i<order; i++){
        fill(q[i].begin(), q[i].end(), 0.);
        for(int j = 0; j<n_stages+1; j++){
            add_vector_inplace(q[i], p[j][i-1], k[j]);
        }
        multiply_scalar_inplace(q[i], pow(h,1-i));
    }
    add_vector_inplace(y, 1., dy);
    t += dt;
    h = h_new;
    k[0] = k[n_stages];
}

template <typename model, exec<model> call>
void Integrator<model,call>::Interpolate(double t, vector<double>& y){
    double dt = dir*(t-t_old);
    y = q[0];
    for(uint i=1; i<order; i++){
        add_vector_inplace(y, dir*pow(dt,i), q[i]);
    }
}
