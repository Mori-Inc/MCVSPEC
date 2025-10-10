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
    private:
        const model* func;
        double abs_err, rel_err;
        valarray<double> k[n_stages+1];
        valarray<double> q[order];
        valarray<double> tol;
        double dir;
        double h;
        double t_old,h_old;

    public:
        explicit Integrator(const model&, const double absolute_err=1e-8, const double relative_err=1e-6);
        void Initialize(double&, const double, const valarray<double>&);
        void Integrate(double&, const double, valarray<double>&);
        void Step(double&, valarray<double>&);
        void Dense_Step(double&, valarray<double>&);
        void Interpolate(double, valarray<double>&);
    private:
        void Prepare_Step(const double&, double&, const valarray<double>&, valarray<double>&, double&);
        void Set_Initial_Step(const double, const valarray<double>);
        void Dense_Output(const double, const double, const valarray<double>, const double);
};

#include "integration.tpp"
