#include "integration.hh"
using std::find;
using std::fill;
using std::begin;
using std::end;
using std::abs;
using std::max;
using std::min;
using std::pow;

using namespace tableau;

double norm(valarray<double> x){
    return sqrt((x*x).sum()/x.size());
}
valarray<double> element_max(valarray<double> x, valarray<double> y){
    valarray<double> max_arr(x.size());
    for(int i = 0; i<x.size(); i++){
        max_arr[i] = max(x[i],y[i]);
    }
    return max_arr;
}

Integrator::Integrator(){}

Integrator::Integrator(valarray<double> (*function)(double, valarray<double>, void*), const int n_dims){
    double pars = 1;
    func = {function, &pars};
    n_dim = n_dims;
    abs_err = absolute_err;
    rel_err = relative_err;
}

Integrator::Integrator(valarray<double> (*function)(double, valarray<double>, void*), void* pars, const int n_dims){
    func = {function, pars};
    n_dim = n_dims;
    abs_err = absolute_err;
    rel_err = relative_err;
}

Integrator::Integrator(valarray<double> (*function)(double, valarray<double>, void*), const int n_dims, const double absolute_err, const double relative_err){
    double pars = 1;
    func = {function, &pars};
    n_dim = n_dims;
    abs_err = absolute_err;
    rel_err = relative_err;
}

void Integrator::Set_Initial_Step(const double dir, const double t0, const valarray<double> y0){
    valarray<double> tol = abs_err + rel_err*abs(y0);
    double h_0 = 1e-2*norm(y0)/norm(k[0]);
    valarray<double> f_1 = func(t0+dir*h_0, y0+dir*h_0*k[0]);
    double delta = norm((f_1-k[0])/tol)/h_0;
    double h_1 = pow(1e-2/max(delta,norm(k[0]/tol)),1./order);
    h = min(1e2*h_0,h_1);
}

void Integrator::Integrate(void* parameters, const double t_start, const double t_end, const valarray<double> y_start, bool interpolate, bool y_bounds){
    func.pars = parameters;
    fill(begin(k), end(k), valarray<double>(0.,n_dim));
    fill(begin(q), end(q), valarray<double>(0.,n_dim));
    t.resize(1);
    y.resize(1);
    t[0] = t_start;
    y[0] = y_start;

    double dir = (0. < (t_end-t[0])) - ((t_end-t[0]) < 0.);
    k[0] = func(t[0], y[0]);
    Set_Initial_Step(dir, t[0], y[0]);

    valarray<double> y_new(n_dim);
    double t_new;

    bool before_bound = true;
    double bound_dir;
    if(y_bounds){
        bound_dir = (y_boundary-y_start[boundary_index]>0)-(y_boundary-y_start[boundary_index]<0);
    }

    while(dir*(t_end-t.back()) > 0 && t.size() < max_itter){
        h = min(h,dir*(t_end-t.back()));
        h = Step(dir, t.back(), y.back(), &t_new, &y_new);

        if(interpolate){
            Dense_Output(dir, t.back(), y.back(), t_new);
        }
        if(y_bounds){
            before_bound = bound_dir*(y_boundary-y_new[boundary_index]) > absolute_err;
        }

        y.push_back(y_new);
        t.push_back(t_new);
        k[0] = k[n_stages];
        if(!before_bound){
            return;
        }
    }
}

void Integrator::Integrate(void* parameters, const double t_start, const double t_end, const valarray<double> y_start){
    Integrate(parameters, t_start, t_end, y_start, false, false);
}

void Integrator::Integrate(void* parameters, const double t_start, const double t_end, const valarray<double> y_start, const vector<double> t_evals){

    t_eval = valarray<double>(t_evals.data(), t_evals.size());
    y_eval.resize(t_eval.size());
    Integrate(parameters, t_start, t_end, y_start, true, false);
}

void Integrator::Integrate(void* parameters, const double t_start, const double t_end, const valarray<double> y_start, const double y_bound, const int bound_ind){
    y_boundary = y_bound;
    boundary_index = bound_ind;
    Integrate(parameters, t_start, t_end, y_start, false, true);
}

double Integrator::Step(double dir, const double t_old, const valarray<double> y_old, double* t_new, valarray<double>* y_new){
    bool step_succeded = false;
    bool step_failed = false;
    valarray<double> dy(n_dim), err_arr(n_dim);
    double h_new, err_norm;
    while(!step_succeded){
        *t_new = t_old+dir*h;
        for(int i = 1; i<n_stages; i++){
            dy = 0;
            for(int j = 0; j<i; j++){
                dy += a[i][j]*k[j];
            }
            k[i] = func(t_old+c[i]*dir*h, y_old+dir*h*dy);
        }
        *y_new = 0;
        err_arr = 0;
        for(int i = 0; i<n_stages; i++){
            *y_new += b[i]*k[i];
            err_arr += e[i]*k[i];
        }
        *y_new = y_old + dir*h*(*y_new);
        k[n_stages] = func(*t_new, *y_new);
        tol = abs_err + rel_err*element_max(abs(y_old), abs(*y_new));
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
    return h_new;
}

void Integrator::Dense_Output(const double dir, const double t_old, const valarray<double> y_old, const double t_new){
    double sigma, h;
    h = abs(t_new-t_old);

    valarray<double> dy(n_dim);
    for(int i = 0; i<order-1; i++){
        q[i] = 0.;
        for(int j = 0; j<n_stages+1; j++){
            q[i] += k[j]*p[j][i];
        }
    }
    for(double t_interp:t_eval){
        if((dir*(t_new-t_interp)<0)|(dir*(t_interp-t_old)<0)){
            continue;
        }
        sigma = (t_interp-t_old)/(dir*h);;
        dy = 0.;
        for(int i = 0; i < order-1; i++){
            dy += q[i]*pow(sigma,i+1);
        }
        dy *= dir*h;
        int ind = find(begin(t_eval), end(t_eval), t_interp) - begin(t_eval);
        y_eval[ind] = y_old + dy;
    }
}
