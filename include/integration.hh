#pragma once

#include <array>
#include "tableau.hh"

template <size_t n_dim>
using State = std::array<double, n_dim>;

template <size_t n_dim, class RHS>
class Integrator{
    private:
        RHS func;
        double abs_err = 1e-8;
        double rel_err = 1e-6;

        std::array<State<n_dim>, tableau::n_stages+1> k{};
        std::array<State<n_dim>, tableau::order> q{};
        double t_internal = 0.0;
        State<n_dim> y_internal{};
        double dt = 0.0;
        State<n_dim> dy{};
        State<n_dim> err_arr{};
        State<n_dim> tol{};
        double dir = 1.;
        double h = 0.0;
        double t_old = 0.0;

        State<n_dim> buffer{};

    public:
        explicit Integrator(RHS, const double absolute_err = 1e-8, const double relative_err = 1e-6);
        void Initialize(double&, const double, const State<n_dim>&);
        void Step(double&, State<n_dim>&);
        void Dense_Step(double&, State<n_dim>&);
        void Interpolate(double, State<n_dim>&);
    private:
        void Prepare_Step(const double&, const State<n_dim>&, double&);
        void Set_Initial_Step(const double&, const State<n_dim>&);
        void Dense_Output(const double, const double, const State<n_dim>, const double);
};

#include "integration.tpp"
