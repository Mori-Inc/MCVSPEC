#pragma once

#include "tableau.hh"
#include <valarray>
#include <vector>

using std::valarray;
using std::vector;

using tableau::n_stages;
using tableau::order;

const int max_itter = 1000000;

template <typename model>
using exec = void (model::*)(double, const valarray<double>&, valarray<double>&) const;
template <typename model, exec<model> call>
class Integrator{
    public:
        vector<double> t;
        vector<valarray<double>> y;
        valarray<double> t_eval;
        valarray<valarray<double>> y_eval;
    private:
        const model* func;
        double abs_err, rel_err;
        valarray<double> k[n_stages+1];
        valarray<double> q[order-1];
        valarray<double> tol;
        double h;
        double y_boundary;
        int boundary_index;

    public:
        explicit Integrator(const model&, const double absolute_err=1e-8, const double relative_err=1e-6);
        void Integrate(const double, const double, const valarray<double>);
        void Integrate(const double, const double, const valarray<double>, const vector<double>);
        void Integrate(const double, const double, const valarray<double>, const double, const int);
    private:
        void Integrate(const double, const double, const valarray<double>, bool, bool);
        double Step(double, const double, const valarray<double>, double*, valarray<double>*);
        void Set_Initial_Step(const double, const double, const valarray<double>);
        void Dense_Output(const double, const double, const valarray<double>, const double);
};

#include "integration.tpp"
