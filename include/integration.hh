#pragma once

#include "vector_operators.hh"
#include "tableau.hh"
#include <vector>

using std::vector;

using tableau::n_stages;
using tableau::order;

const int max_itter = 1000000;

template <typename model>
using exec = void (model::*)(double, const vector<double>&, vector<double>&) const;
template <typename model, exec<model> call>
class Integrator{
    private:
        const model* func;
        double abs_err, rel_err;
        vector<double> k[n_stages+1];
        vector<double> q[order];
        vector<double> buffer;
        vector<double> y_internal;
        double t_internal;
        vector<double> tol;
        double dir;
        double h;
        double t_old;

    public:
        explicit Integrator(const model&, const int, const double absolute_err=1e-8, const double relative_err=1e-6);
        void Initialize(double&, const double, const vector<double>&);
        int Integrate(double&, const double, vector<double>&);
        void Step(double&, vector<double>&);
        void Dense_Step(double&, vector<double>&);
        void Interpolate(double, vector<double>&);
    private:
        void Prepare_Step(const double&, double&, const vector<double>&, vector<double>&, double&);
        void Set_Initial_Step(const double&, const vector<double>&);
        void Dense_Output(const double, const double, const vector<double>, const double);
};

#include "integration.tpp"
