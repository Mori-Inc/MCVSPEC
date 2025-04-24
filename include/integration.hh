#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "tableau.hh"
#include <cmath>
#include <valarray>
#include <vector>

using std::valarray;
using std::vector;

const double absolute_err = 1e-8;
const double relative_err = 1e-6;
const double max_itter = 100000;

double norm(valarray<double>); // computes norm of array
valarray<double> element_max(valarray<double>, valarray<double>); // returns max of to arrays element wise (i.e. return[i] is the max of (array_a[i], array_b[i]))

struct equation{
    valarray<double> (*func)(double, valarray<double>, void*);
    void* pars;
    valarray<double> operator()(double t, valarray<double> y) const {
        return func(t,y,pars);
    }
};

class Integrator{
    public:
        vector<double> t;
        vector<valarray<double>> y;
        vector<valarray<double>> deriv;

        valarray<double> t_eval;
        valarray<valarray<double>> y_eval;
    private:
        equation func;
        double abs_err, rel_err;
        int n_dim;
        valarray<double> k[n_stages+1];
        valarray<double> q[order-1];
        valarray<double> tol;
        double h;

    public:
        Integrator();
        Integrator(valarray<double> (*func)(double, valarray<double>, void*), const int);
        Integrator(valarray<double> (*func)(double, valarray<double>, void*), void*, const int);
        Integrator(valarray<double> (*func)(double, valarray<double>, void*), const int, const double, const double);
        void Integrate(void*, const double, const double, const valarray<double>);
        void Integrate(void*, const double, const double, const valarray<double>, const vector<double>);
        void Integrate(void*, const double, const double, const valarray<double>, const valarray<double>);
    private:
        double Step(double, const double, const valarray<double>, double*, valarray<double>*);
        void Set_Initial_Step(const double, const double, const valarray<double>);
        void Dense_Output(const double, const double, const valarray<double>, const double);
};

#endif
