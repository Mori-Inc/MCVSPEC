#pragma once

#include <array>
#include "tableau.hh"

#ifndef MCVSPEC_MAX_STEPS
#define MCVSPEC_MAX_STEPS 1000000
#endif

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
        long n_steps = 0; // steps taken since last Initialize()

        State<n_dim> buffer{};

    public:
        // Hard cap on accepted steps per integration. Healthy MCVSPEC solutions need far fewer;
        // solutions creeping toward a sonic point can otherwise take ~1e9 tiny steps (an apparent hang).
        static constexpr long max_steps = MCVSPEC_MAX_STEPS;
        explicit Integrator(RHS, const double absolute_err = 1e-8, const double relative_err = 1e-6);
        void Initialize(double&, const double, const State<n_dim>&);
        int Step(double&, State<n_dim>&);
        int Dense_Step(double&, State<n_dim>&);
        void Interpolate(double, State<n_dim>&);
    private:
        int Prepare_Step(const double&, const State<n_dim>&, double&);
        void Set_Initial_Step(const double&, const State<n_dim>&);
        void Dense_Output(const double, const double, const State<n_dim>, const double);
};

#include "integration.tpp"
