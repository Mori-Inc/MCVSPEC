#pragma once
#include "integration.hh"
#include <valarray>
using std::fill;
using std::begin;
using std::end;
using std::abs;
using std::max;
using std::min;
using std::pow;

template <typename model>
using exec = void (model::*)(double, const valarray<double>&, valarray<double>&) const;

inline double norm(valarray<double> x){
    return sqrt((x*x).sum()/x.size());
}
inline valarray<double> element_max(valarray<double> x, valarray<double> y){
    valarray<double> max_arr(x.size());
    for(int i = 0; i<x.size(); i++){
        max_arr[i] = max(x[i],y[i]);
    }
    return max_arr;
}

template <typename model, exec<model> call>
Integrator<model,call>::Integrator(const model& function, const double absolute_err, const double relative_err):
    func(&function), abs_err(absolute_err), rel_err(relative_err)
{}

template <typename model, exec<model> call>
void Integrator<model,call>::Set_Initial_Step(const double t0, const valarray<double> y0){
    valarray<double> tol = abs_err + rel_err*abs(y0);
    double h_0 = 1e-2*norm(y0)/norm(k[0]);
    (func->*call)(t0+dir*h_0, y0+dir*h_0*k[0], k[1]);
    double delta = norm((k[1]-k[0])/tol)/h_0;
    double h_1 = pow(1e-2/max(delta,norm(k[0]/tol)),1./order);
    h = min(1e2*h_0,h_1);
}

template <typename model, exec<model> call>
void Integrator<model,call>::Initialize(double& t, const double t_end, const valarray<double>& y){
    fill(begin(k), end(k), valarray<double>(0.,y.size()));
    fill(begin(q), end(q), valarray<double>(0.,y.size()));
    dir = (0. < (t_end-t)) - ((t_end-t) < 0.);
    (func->*call)(t, y, k[0]);
    Set_Initial_Step(t, y);
}

template <typename model, exec<model> call>
void Integrator<model,call>::Integrate(double& t, const double t_end, valarray<double>& y){
    Initialize(t, t_end, y);
    int n_steps=0;

    while(dir*(t_end-t) > 0 && n_steps < max_itter){
        h = min(h,dir*(t_end-t));
        Step(t, y);
    }
}

template <typename model, exec<model> call>
void Integrator<model,call>::Prepare_Step(const double& t, double& dt,const valarray<double>& y, valarray<double>& dy, double& h_new){
    using tableau::a;
    using tableau::b;
    using tableau::c;
    using tableau::e;
    bool step_succeded = false;
    bool step_failed = false;
    valarray<double> err_arr(y.size());
    double err_norm;
    while(!step_succeded){
        err_arr = e[0]*k[0];
        for(uint i=1; i<n_stages; i++){
            dy = 0;
            for(int j = 0; j<i; j++){
                dy += a[i][j]*k[j];
            }
            (func->*call)(t+c[i]*dir*h, y+dir*h*dy, k[i]);
            err_arr += e[i]*k[i];
        }
        dy = 0;
        for(uint i=0; i<n_stages; i++){
            dy += b[i]*k[i];
        }
        dy *= dir*h;
        dt = dir*h;

        (func->*call)(t+dt, y+dy, k[n_stages]);
        tol = abs_err + rel_err*element_max(abs(y), abs(y+dy));
        err_norm = norm((err_arr+e[n_stages]*k[n_stages])*h/tol);
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
void Integrator<model,call>::Step(double& t, valarray<double>& y){
    double dt=h, h_new=h;
    valarray<double> dy(y.size());
    Prepare_Step(t, dt, y, dy, h_new);
    y += dy;
    t += dt;
    h = h_new;
    k[0] = k[n_stages];
}

template <typename model, exec<model> call>
void Integrator<model,call>::Dense_Step(double& t, valarray<double>& y){
    using tableau::p;
    double dt=h, h_new=h;
    valarray<double> dy(y.size());
    t_old = t;
    q[0]=y;
    Prepare_Step(t, dt, y, dy, h_new);
    for(int i = 1; i<order; i++){
        q[i] = 0.;
        for(int j = 0; j<n_stages+1; j++){
            q[i] += k[j]*p[j][i-1];
        }
        q[i] /= pow(h,i-1);
    }
    y += dy;
    t += dt;
    h = h_new;
    k[0] = k[n_stages];
}

template <typename model, exec<model> call>
void Integrator<model,call>::Interpolate(double t, valarray<double>& y){
    double dt = dir*(t-t_old);
    y = q[0];
    for(uint i=1; i<order; i++){
        y += dir*q[i]*pow(dt,i);
    }
}
