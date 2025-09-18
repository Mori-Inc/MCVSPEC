#pragma once

class Dipole{
    public:
        void Solve_Coordinates(double);
        double Get_Radial_Distance();
        double Get_Projection();
        double Get_Metric(int);
        double Get_Convergence();
    private:
        double u;
        double r;
        double proj_r_w;
        double metric[3];
        double cannalle_vars[3];
        double convergance;
};
