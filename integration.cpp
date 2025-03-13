#include "integration.hh"
using std::find, std::fill, std::begin, std::end;
using std::abs, std::max, std::min;
using std::pow;

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
Integrator::Integrator(function<valarray<double>(double, valarray<double>, void*)> function, const int n_dims){
    double pars = 1;
    func = {function, &pars};
    n_dim = n_dims;
    abs_err = absolute_err;
    rel_err = relative_err;
}

Integrator::Integrator(function<valarray<double>(double, valarray<double>, void*)> function, const int n_dims, const double absolute_err, const double relative_err){
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

void Integrator::Integrate(void* parameters, const double t_start, const double t_end, const valarray<double> y_start){
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

    while(dir*(t_end-t.back()) > 0 && t.size() < max_itter){
        h = min(h,dir*(t_end-t.back()));
        h = Step(dir, t.back(), y.back(), &t_new, &y_new);
        y.push_back(y_new);
        t.push_back(t_new);
        k[0] = k[n_stages];
    }
}

void Integrator::Integrate(void* parameters, const double t_start, const double t_end, const valarray<double> y_start, const valarray<double> y_bound){
    func.pars = parameters;

    fill(begin(k), end(k), valarray<double>(0.,n_dim));
    fill(begin(q), end(q), valarray<double>(0.,n_dim));
    t.resize(1);
    y.resize(1);
    t[0] = t_start;
    y[0] = y_start;

    bool before_bound = true;
    valarray<double> bound_dir(n_dim);
    for(int i = 0; i < n_dim; i++){
        bound_dir[i] = (0. < (y_bound[i]-y_start[i])) - ((y_bound[i]-y_start[i]) < 0.);
    }

    double dir = (0. < (t_end-t[0])) - ((t_end-t[0]) < 0.);
    k[0] = func(t[0], y[0]);
    Set_Initial_Step(dir, t[0], y[0]);
    valarray<double> y_new(n_dim);
    double t_new;
    while(dir*(t_end-t.back()) > 0 && t.size() < max_itter && before_bound){
        h = min(h,dir*(t_end-t.back()));
        h = Step(dir, t.back(), y.back(), &t_new, &y_new);

        for(int i = 0; i < n_dim; i++){
            if(bound_dir[i]*(y_bound[i]-y_new[i]) < 0.){
                before_bound = false;
            }
        }

        y.push_back(y_new);
        t.push_back(t_new);
        k[0] = k[n_stages];
        deriv.push_back(k[0]);
    }
}

void Integrator::Integrate(void* parameters, const double t_start, const double t_end, const valarray<double> y_start, const vector<double> t_evals){
    func.pars = parameters;
    t = t_evals;
    y.resize(t_evals.size());

    fill(begin(k), end(k), valarray<double>(0.,n_dim));
    fill(begin(q), end(q), valarray<double>(0.,n_dim));

    double dir = (0. < (t_end-t[0])) - ((t_end-t[0]) < 0.);
    k[0] = func(t_start, y_start);
    Set_Initial_Step(dir, t_start, y_start);
    valarray<double> y_old(n_dim), y_new(n_dim);
    double t_old, t_new, h_new;
    y_old = y_start;
    t_old = t_start;

    while(dir*(t_end-t_old) > 0 && t.size() < max_itter){
        h = min(h,dir*(t_end-t_old));
        h_new = Step(dir, t_old, y_old, &t_new, &y_new);
        Dense_Output(dir, t_old, y_old, t_new);
        h = h_new;
        y_old = y_new;
        t_old = t_new;
        k[0] = k[n_stages];
    }
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
    double sigma;
    valarray<double> dy(n_dim);
    for(int i = 0; i<order-1; i++){
        q[i] = 0.;
        for(int j = 0; j<n_stages+1; j++){
            q[i] += k[j]*p[j][i];
        }
    }
    for(double t_interp:t){
        if((dir*(t_new-t_interp)<0)|(dir*(t_interp-t_old)<0)){
            continue;
        }
        sigma = (t_interp-t_old)/(dir*h);;
        dy = 0.;
        for(int i = 0; i < order-1; i++){
            dy += q[i]*pow(sigma,i+1);
        }
        dy *= dir*h;
        int ind = find(t.begin(), t.end(), t_interp) - t.begin();
        y[ind] = y_old + dy;
    }
}
