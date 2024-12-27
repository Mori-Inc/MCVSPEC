#include "mcvspec_core.hh"
#include "mcvspec_namespaces.hh"

using namespace mcvspec;
using std::valarray;
using std::pow;
using std::abs;
using std::sqrt;
using std::cbrt;
using std::cerr;
using std::endl;

namespace mcvspec{
    const int max_num_grid_points = 200;
    valarray<double> density(max_num_grid_points);
    valarray<double> pressure(max_num_grid_points);
    valarray<double> temperature(max_num_grid_points);
    valarray<double> altitude(max_num_grid_points);
    valarray<double> velocity(max_num_grid_points);
    valarray<double> electron_dens(max_num_grid_points);
    valarray<double> offset_altitude(max_num_grid_points);
}

double Calculate_Magnetic_Field(){
    return pow(2.,5./2.)*pow(mag_radius,7./4.)*pow(wd_radius,3)*pow(grav_const*mass*accretion_rate*accretion_rate,1./4.);
}
double Calculate_White_Dwarf_Radius(){
    return solar_radius*(0.0225/wd_mol_mass)*sqrt(1.0-pow(mass/mass_limit,4./3.))/cbrt(mass/mass_limit);
}
double Calculate_Accretion_Rate(){
    return luminosity/(grav_const*mass*(1./wd_radius - 1./mag_radius)*4*pi*wd_radius*wd_radius*accretion_area);
}
double Estimate_Shock_Height(){
    double integral = (39.*sqrt(3.) - 20*pi)/96.;
    return pow(free_fall_velocity,3.)*integral/(2.*bremss_const*accretion_rate);
}
double Calculate_B_Free_Shock_Height(){
    return (pow(free_fall_velocity,2)/(2*(accretion_rate/free_fall_velocity)*bremss_const))*((39.*sqrt(3.)-20.*pi)/96.);
}
double Calculate_Shock_Temperature(){
    return (3./8.)*grav_const*mass*col_mol_mass*hydrg_mass*(1/(wd_radius+shock_height) - 1/mag_radius)/boltz_const;
}
double Calculate_Electron_Density(){
    return (accretion_rate/free_fall_velocity)*shock_velocity*electron_ion_ratio/(hydrg_mass*col_mol_mass);
}
double Calculate_Epsilon_Zero(){
    return 9.1e-3*pow(b_field/10.,2.85)*pow(shock_temperature/1e8,2.)*pow(shock_electron_dens/1e16,-1.85)*pow(b_free_shock_height/1e7,-0.85);
}

double Epsilon_Diff(double epsilon_s, void* eps_zero){
    double epsilon_zero = *(double *)eps_zero;
    return epsilon_zero*pow(1. + epsilon_s, 0.425) - epsilon_s;
}
double Epsilon_Diff_Derivative(double epsilon_s, void* eps_zero){
    double epsilon_zero = *(double *)eps_zero;
    return 0.425*epsilon_zero*pow(1. + epsilon_s, -0.575) - 1.;
}
void Epsilon_Diff_and_Derivative(double epsilon_s, void* eps_zero, double* diff, double* diff_deriv){
    *diff = Epsilon_Diff(epsilon_s, eps_zero);
    *diff_deriv = Epsilon_Diff_Derivative(epsilon_s, eps_zero);
}

gsl_function_fdf gsl_func = {&Epsilon_Diff, &Epsilon_Diff_Derivative, &Epsilon_Diff_and_Derivative, &epsilon_zero};
gsl_root_fdfsolver *newton_solver = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_newton);

double Root_Finder(double initial_guess, int max_itter, double tolerance){
    shock_electron_dens = Calculate_Electron_Density(); // density times n_electrons/gram
    shock_temperature = Calculate_Shock_Temperature();
    b_free_shock_height = Calculate_B_Free_Shock_Height(); // equation 7 from Wu 1994 (https://ui.adsabs.harvard.edu/link_gateway/1994ApJ...426..664W/doi:10.1086/174103)
    epsilon_zero = Calculate_Epsilon_Zero();

    gsl_func.params = &epsilon_zero;
    gsl_root_fdfsolver_set(newton_solver, &gsl_func, initial_guess);
    int i = 0;
    double epsilon_s = initial_guess;
    while (gsl_root_test_residual(Epsilon_Diff(epsilon_s, gsl_func.params), tolerance) == GSL_CONTINUE && i < max_itter){
        i += 1;
        gsl_root_fdfsolver_iterate(newton_solver);
        epsilon_s = gsl_root_fdfsolver_root(newton_solver);
    }
    if (gsl_root_test_residual(Epsilon_Diff(epsilon_s, gsl_func.params), tolerance) != GSL_SUCCESS){
        cerr << "Newton solver failed to converge after " << max_itter << " iterations." << endl;
    }
    gsl_root_fdfsolver_free(newton_solver);
    return epsilon_s;
}

int Normalized_Position_Derivative(double velocity_normed, const double position_normed[], double pos_derivative[], void* eps_s){
    double epsilon_s = *(double *)eps_s;
    double prefactor = pow(free_fall_velocity,3)/(2.*bremss_const*accretion_rate);
    pos_derivative[0] = prefactor*(pow(velocity_normed,1.5)*(5.0-8.*velocity_normed)/sqrt(1-velocity_normed)/
                        (1.+epsilon_s*pow(3.,-alpha)*pow(4.,alpha+beta)*
                         pow(1.0-velocity_normed,alpha)*pow(velocity_normed,beta)));
    return GSL_SUCCESS;
}

gsl_odeiv2_system accretion_column = {Normalized_Position_Derivative, nullptr, 1, &epsilon_s};
const gsl_odeiv2_step_type *ode_type = gsl_odeiv2_step_rkf45;
gsl_odeiv2_step *ode_stepper = gsl_odeiv2_step_alloc(ode_type, 1);
gsl_odeiv2_control *ode_control = gsl_odeiv2_control_y_new(absolute_err, relative_err);
gsl_odeiv2_evolve *ode_evolver = gsl_odeiv2_evolve_alloc(1);

void Runge_Kutta(double initial_veloicty, double initial_height){
    double vel = initial_veloicty;
    double height[1] = {initial_height};
    int n = 0;
    double initial_step = absolute_err;
    velocity[0] = vel*free_fall_velocity;
    altitude[0] = height[0]*shock_height;
    while (vel < shock_velocity){
        int status = gsl_odeiv2_evolve_apply(ode_evolver, ode_control, ode_stepper, &accretion_column,
                                             &vel, shock_velocity ,&initial_step, height);
        if(status != GSL_SUCCESS){
            break;
        }
        velocity[n] = vel*free_fall_velocity;
        altitude[n] = height[0]*shock_height;
        n += 1;
    }
    gsl_odeiv2_evolve_free(ode_evolver);
    gsl_odeiv2_control_free(ode_control);
    gsl_odeiv2_step_free(ode_stepper);

    for (int i = 0; i < n; i++){
        density[i] = accretion_rate/velocity[i];
        pressure[i] = accretion_rate*(free_fall_velocity - velocity[i]);
        temperature[i] = boltz_const_kev*shock_temperature*velocity[i]*(free_fall_velocity-velocity[i])*(16./(3.*pow(free_fall_velocity, 2)));
        electron_dens[i] = density[i]*electron_ion_ratio/(hydrg_mass*col_mol_mass);
    }
    num_grid_points =  n;
}

void Shock_Height_Shooting(double tolerance, int max_itter){
    double old_height, error;

    error = 1e100;
    int i = 0;

    Runge_Kutta(initial_veloicty, initial_height);
    shock_height = altitude[num_grid_points-1];
    old_height = altitude[num_grid_points-1];

    while(i < max_itter && error > tolerance){
        free_fall_velocity = sqrt(2.*grav_const*mass*(1/(wd_radius+shock_height) - 1/mag_radius));
        epsilon_s = Root_Finder(1e5, 100000, 1e-6);

        Runge_Kutta(initial_veloicty, initial_height);

        shock_height = altitude[num_grid_points-1];
        if (wd_radius + shock_height > mag_radius) {
            shock_height = old_height;
            cerr << "Accretion column height has exceded magnetospheric radius" << endl;
        }
        error = abs((old_height-shock_height)/old_height);
        old_height = shock_height;
        i++;
    }
}

void MCVspec_Spectrum(const RealArray& energy, const int spectrum_num, RealArray& flux, const string& init_string){

    int n = flux.size();

    double segment_height = altitude[0];

    valarray<double> flux_from_layer(n);
    valarray<double> flux_error(n);

    valarray<double> apec_parameters(3);
    valarray<double> refl_parameters(5);

    refl_parameters[0] = 1-sqrt(1.0-1.0/pow(1+altitude[num_grid_points-1]/wd_radius,2));
    refl_parameters[1] = 0.;
    refl_parameters[2] = wd_abund;
    refl_parameters[3] = wd_abund;
    refl_parameters[4] = cos_incl;

    for(int i=0; i<num_grid_points; i++){
        for(int j=0; j<n; j++){
            flux[j] += flux_from_layer[j];
            flux_from_layer[j] = 0;
        }

        apec_parameters[0] = temperature[i]*boltz_const_kev;
        apec_parameters[1] = col_abund;
        apec_parameters[2] = 0.0;

        if (apec_parameters[0] > 64.0){
            CXX_bremss(energy, apec_parameters, spectrum_num, flux_from_layer, flux_error, init_string);
        }
        else{
            CXX_apec(energy, apec_parameters, spectrum_num, flux_from_layer, flux_error, init_string);
        }

        // normalizes spectrum appropriately
        if(i!=0){
            segment_height = altitude[i]-altitude[i-1];
        }
        for(int j=0; j<n; j++){
            flux_from_layer[j] *= apec_norm*segment_height*pow(electron_dens[i]*1e-7,2)/electron_ion_ratio;
        }

        // apply reflect to each slice
        if (reflection_sel == 2){
            refl_parameters[0] = 1-sqrt(1.0-1.0/pow(1+altitude[i]/wd_radius,2));
            CXX_reflect(energy, refl_parameters, spectrum_num, flux_from_layer, flux_error, init_string);
        }
    }

    for(int j=0; j<n; j++){
        flux[j] += flux_from_layer[j];
        flux_from_layer[j] = 0;
    }

    if(reflection_sel==1){
        CXX_reflect(energy, refl_parameters, spectrum_num, flux, flux_error, init_string);
    }
}
