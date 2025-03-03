#include "integration.hh"

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

void Dormand_Prince(function<valarray<double>(double, valarray<double>, void*)> func, void* args, vector<double>* t, vector<valarray<double>>* y, double t_bound, vector<double>* t_eval, vector<valarray<double>>* y_eval, double abs_err, double rel_err, int max_itter){
    int n_dim = y[0].size();
    valarray<double> k[n_stages+1];
    fill(begin(k), end(k), valarray<double>(0.,n_dim));
    valarray<double> q[order-1];
    fill(begin(q), end(q), valarray<double>(0.,n_dim));
    y_eval->resize(t_eval->size());

    double dir = (0. < (t_bound-(*t)[0])) - ((t_bound-(*t)[0]) < 0.); // direction of integration
    k[0] = func((*t)[0], (*y)[0],args);
    valarray<double> tol = abs_err + rel_err*abs((*y)[0]);
    double h_0 = 1e-2*norm((*y)[0])/norm(k[0]);
    valarray<double> f_1 = func((*t)[0]+dir*h_0, (*y)[0]+dir*h_0*k[0], args);
    double delta = norm((f_1-k[0])/tol)/h_0;
    double h_1 = pow(1e-2/max(delta,norm(k[0]/tol)),1./order);
    double h = min(1e2*h_0,h_1);

    valarray<double> y_new(n_dim), dy(n_dim), err(n_dim);
    double t_new, err_norm, factor, sigma;
    bool step_failed = false;

    while(dir*(t_bound-t->back()) > 0 && t->size() < max_itter){
        h = min(h,dir*(t_bound-t->back()));
        t_new = t->back()+dir*h;

        for(int i = 1; i<n_stages; i++){
            dy = 0;
            for(int j = 0; j<i; j++){
                dy += a[i][j]*k[j];
            }
            k[i] = func(t->back()+c[i]*dir*h, y->back()+dir*h*dy, args);
        }
        y_new = 0;
        err = 0;
        for(int i = 0; i<n_stages; i++){
            y_new += b[i]*k[i];
            err += e[i]*k[i];
        }
        y_new = y->back() + dir*h*y_new;
        k[n_stages] = func(t_new, y_new, args);
        tol = abs_err + rel_err*element_max(abs(y->back()), abs(y_new));
        err_norm = norm((err+e[n_stages]*k[n_stages])*h/tol);

        if(err_norm < 1.){
            if (err_norm == 0){
                factor = 10.;
            }
            else if(step_failed){
                factor = min(0.9*pow(err_norm,-1./order), 1.);
            }
            else{
                factor = min(0.9*pow(err_norm,-1./order), 5.);
            }
            step_failed = false;
        }
        else if(err_norm >= 1.){
            step_failed = true;
            h *= max(0.9*pow(err_norm,-1./order), 0.2);
            continue;
        }
        else{
            // if err is nan
            h *= 0.2;
            step_failed = true;
            continue;
        }

        if(t_eval->size() > 0){
            for(int i = 0; i<order-1; i++){
                q[i] = 0.;
                for(int j = 0; j<n_stages+1; j++){
                    q[i] += k[j]*p[j][i];
                }
            }
            for(double t_interp:*t_eval){
                if((dir*(t_new-t_interp)<0)|(dir*(t_interp-t->back())<0)){
                    continue;
                }
                sigma = (t_interp-t->back())/(dir*h);;
                dy = 0.;
                for(int i = 0; i < order-1; i++){
                    dy += q[i]*pow(sigma,i+1);
                }
                dy *= dir*h;
                int ind = find(t_eval->begin(), t_eval->end(), t_interp) - t_eval->begin();
                (*y_eval)[ind] = y->back() + dy;
            }
        }

        h *= factor;
        y->push_back(y_new);
        t->push_back(t_new);
        k[0] = k[n_stages];
    }

    if (t->size() == max_itter){
        cout << "Warning steps exceded max_itter, completion is not guaranteed. final (t,y) = (" << t->back() << ", " << y->back()[0] << ')' << endl;
    }
}
