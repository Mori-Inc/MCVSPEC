#include "Cataclysmic_Variable.hh"
#include "constants.hh"
#include "integration.hh"
#include "mass_radius.hh"
#include "gaunt.hh"
#include <cmath>
#include <iostream>
#include <valarray>

using std::cout;
using std::endl;
using std::abs;

static double previous_shock_height = 0;

Cataclysmic_Variable::Cataclysmic_Variable(double m, double r, double b, double mdot, double inv_r_m, double corot_rat, double area, double theta, double n, double dist, int reflection):
    mass(m), radius(r), b_field(b),  inverse_mag_radius(inv_r_m), corotation_ratio(corot_rat), distance(dist), accretion_rate(mdot), accretion_area(area), pressure_ratio(.75), incl_angle(theta), area_exponent(n),  refl(reflection)
{
    if(inverse_mag_radius>0){
        b_field = sqrt(32*accretion_rate*sqrt(grav_const*mass/pow(inverse_mag_radius,7)))/(radius*radius*radius);
    }
}

void Cataclysmic_Variable::Set_Cooling_Constants(){ // "constant" insofar as these values depend only on the input properties not on any derived properties
    avg_ion_mass = (abundances*atomic_mass).sum()*amu_to_g;
    avg_atomic_charge = (abundances*atomic_charge).sum();
    double avg_charge_squared = (abundances*atomic_charge*atomic_charge).sum();
    double avg_charge_sqr_over_mass = (abundances*atomic_charge*atomic_charge/atomic_mass).sum()/amu_to_g;
    density_const = avg_atomic_charge/(1 + m_e*avg_atomic_charge/avg_ion_mass);

    force_const = grav_const*mass/(radius*radius);
    cooling_ratio_const = 8.07e-2*avg_atomic_charge*pow(b_field, 2.85)*pow(avg_ion_mass/density_const,3.85);
    cooling_ratio_const /= avg_charge_squared*k_b*k_b*pow(accretion_area,0.425);
    coulomb_log_const = 0.5*log(2*m_e/(pi*alpha*c)) + 1.5*log(avg_ion_mass/(hbar*density_const));
    exchange_const = 4*(alpha*hbar*c)*(alpha*hbar*c)*sqrt(2*pi*m_e*pow((density_const/avg_ion_mass),5))*avg_charge_sqr_over_mass;
    bremss_const = sqrt(512*pi/(27*m_e*m_e*m_e))*alpha*alpha*alpha*hbar*hbar;
    bremss_const *= (avg_charge_squared/avg_atomic_charge)*sqrt(pow(density_const/avg_ion_mass,3));
}

double Cataclysmic_Variable::Get_Accretion_Rate(double luminosity, double mass, double radius, double inverse_mag_radius){
    double accretion_rate = luminosity/(grav_const*mass*((1./radius) - inverse_mag_radius));
    return accretion_rate;
}

double Cataclysmic_Variable::Get_Radius(double mass){
    int left_ind = 0;
    int i = mass_radius_length/2;
    int right_ind = mass_radius_length-1;
    while(right_ind-left_ind > 1){
        i = left_ind + (right_ind-left_ind)/2;
        if(mass>white_dwarf_mass[i]){
            left_ind = i;
        }
        else{
            right_ind = i;
        }
    }
    double delta_r = white_dwarf_radius[right_ind]-white_dwarf_radius[left_ind];
    double delta_m = white_dwarf_mass[right_ind]-white_dwarf_mass[left_ind];
    double radius = white_dwarf_radius[left_ind] + (delta_r/delta_m)*(mass-white_dwarf_mass[left_ind]);
    return radius;
}

void Cataclysmic_Variable::Guess_Shock_Height(){
    if(previous_shock_height != 0){
        Update_Shock_Height(previous_shock_height);
    }
    else{
        double integral = (39.*sqrt(3.) - 20*pi)/96.; // value of integral from EQ 7a of Wu 1994 DOI: 10.1086/174103
        shock_speed = sqrt(2*grav_const*mass*((1./radius) - inverse_mag_radius));
        Update_Shock_Height(pow(shock_speed,3.)*integral*accretion_area/(2*bremss_const*1.2*accretion_rate));
        Update_Shock_Height(pow(shock_speed,3.)*integral*accretion_area/(2*bremss_const*1.2*accretion_rate));
    }
}

void Cataclysmic_Variable::Update_Shock_Height(double h_s){
    shock_height = h_s;
    double rotation_correction = 0.5*corotation_ratio*corotation_ratio*corotation_ratio*inverse_mag_radius;
    rotation_correction *= (radius+shock_height)*(radius+shock_height)*inverse_mag_radius*inverse_mag_radius - 1.;
    shock_speed = sqrt(2*grav_const*mass*((1./(radius+shock_height)) - inverse_mag_radius + rotation_correction));
    shock_mdot = accretion_rate/(accretion_area*pow(1+shock_height/radius, area_exponent));
    non_dim_radius = radius/shock_height;
    cooling_ratio = cooling_ratio_const*pow(shock_speed,5.85)/pow(shock_mdot, 1.85);
}

void Cataclysmic_Variable::Flow_Equation(double vel,const valarray<double>& pos_pres_epres, valarray<double>& derivs) const{
    const double& v = vel;
    const double& x = pos_pres_epres[0];
    const double& p = pos_pres_epres[1];
    const double& pe = pos_pres_epres[2];
    const double& r = non_dim_radius;
    const double& hs = shock_height;
    const double& n = area_exponent;

    const double mdot = pow((1+r)/(x+r), n);
    const double dens = mdot/v;
    const double vff2 = shock_speed*shock_speed;
    const double vff3 = shock_speed*shock_speed*shock_speed;
    const double kT = (avg_ion_mass/density_const)*vff2*pe*v/mdot;
    const double coulomb_log = coulomb_log_const + 2.5*log(shock_speed) - 0.5*log(shock_mdot) + 0.5*log(pe*pe/(dens*dens*dens));
    const double gff = gaunt::gaunt_factor(kT);

    const double grav = (hs/vff2)*force_const*dens/((1+x/r)*(1+x/r));
    const double chi = (1+avg_atomic_charge)/avg_atomic_charge;
    const double exch = (shock_mdot*hs/vff3)*exchange_const*coulomb_log*sqrt(dens*dens*dens*dens*dens/pe)*(p/pe - chi);
    const double cyc = (cooling_ratio/gff)*pe*pe*pow(dens,-3.85)*pow(1+x/r,-8.55-0.425*n);
    const double rad = (shock_mdot*hs/vff3)*bremss_const*gff*sqrt(dens*dens*dens*pe)*(1+cyc);

    double dx_dv = (5*p - 3*mdot*v)/(2*rad + 3*v*grav - 5*p*v*n/(r+x));
    double dp_dv = -grav*dx_dv - mdot;
    double dpe_dv = (2./(3*v))*(rad - exch/vff2)*dx_dv - (5*pe/3.)*(1/v + n*dx_dv/(r+x));

    derivs[0] = dx_dv;
    derivs[1] = dp_dv;
    derivs[2] = dpe_dv;
}

double Cataclysmic_Variable::Get_Landing_Altitude(){
    double t = 0.25;
    valarray<double> y = {1., 0.75, 0.75*(pressure_ratio/(pressure_ratio+1))};
    accretion_column.Integrate(t, 1e-4, y);
    valarray<double> slope(3);
    Flow_Equation(t,  y, slope);
    return y[0] - slope[0]*t;
}

void Cataclysmic_Variable::Shock_Height_Shooting(){
    Update_Shock_Height(shock_height);
    upper_bound = shock_height;
    lower_bound = shock_height;
    double xf_upper = Get_Landing_Altitude();
    double xf_lower = xf_upper;
    if(xf_upper<0){
        upper_bound *= 1.2;
        Update_Shock_Height(upper_bound);
        xf_upper = Get_Landing_Altitude();
        while(xf_upper<0){
            lower_bound = upper_bound;
            xf_lower = xf_upper;
            upper_bound *= 1.2;
            Update_Shock_Height(upper_bound);
            xf_upper = Get_Landing_Altitude();
        }
    }
    else if(xf_lower>0){
        lower_bound *= 0.8;
        Update_Shock_Height(lower_bound);
        xf_lower = Get_Landing_Altitude();
        while(xf_lower>0){
            upper_bound = lower_bound;
            xf_upper = xf_lower;
            lower_bound *= 0.8;
            Update_Shock_Height(lower_bound);
            xf_lower = Get_Landing_Altitude();
        }
    }
    else{
        return;
    }

    double k1 = 0.2/(upper_bound-lower_bound);
    double n0 = 1;
    double nmax = log2((upper_bound-lower_bound)/(2*h_s_tolerance)) + n0;
    int i=0;
    double new_bound, new_altitude, midpoint, regula_falsi, truncation, projection, dir;

    while(upper_bound-lower_bound > 2*h_s_tolerance){
        midpoint = (upper_bound+lower_bound)/2;
        regula_falsi = (xf_upper*lower_bound - xf_lower*upper_bound)/(xf_lower-xf_upper);
        dir = (0. < (midpoint-regula_falsi)) - ((midpoint-regula_falsi) < 0.);
        truncation = k1*(upper_bound-lower_bound)*(upper_bound-lower_bound); // k2 = 2

        if(truncation <= abs(midpoint-regula_falsi)){
            new_bound = regula_falsi + dir*truncation;
        }
        else{
            new_bound = midpoint;
        }

        projection = h_s_tolerance*(pow(2,nmax-i)) - (upper_bound-lower_bound)/2;
        if(abs(new_bound-midpoint) > projection){
            new_bound = midpoint - dir*projection;
        }

        Update_Shock_Height(new_bound);
        new_altitude = Get_Landing_Altitude();
        if(new_altitude>0){
            upper_bound = new_bound;
            xf_upper = new_altitude;
        }
        else if(new_altitude<0){
            lower_bound = new_bound;
            xf_lower = new_altitude;
        }
        else{
            upper_bound = new_altitude;
            lower_bound = new_altitude;
        }
        i++;
    }
    Update_Shock_Height((upper_bound+lower_bound)/2);
    previous_shock_height = shock_height;
}

void Cataclysmic_Variable::Build_Column_Profile(){
    // generate a velocity grid that is has a spacing of ~ kT_grid_spacing
    // interpolate between each point in the RK grid to find a set of velocities to evaluate our integral at

    // determine de-dimensionalized grid size
    const double kTe_const = erg_to_kev*avg_ion_mass*shock_speed*shock_speed/density_const;
    const double kTi_const = avg_atomic_charge*kTe_const;
    const double dkTe = kT_grid_spacing/kTe_const;
    const double dkTi = kT_grid_spacing/kTi_const;
    const double& dx = altitude_grid_spacing;
    const double& r = non_dim_radius;
    const double& n = area_exponent;

    // prep vars for integration
    double t = 0.25;
    valarray<double> y = {1., 0.75, 0.75*(pressure_ratio/(pressure_ratio+1))};
    const double& v = t;
    const double& x = y[0];
    const double& p = y[1];
    const double& pe = y[2];
    accretion_column.Initialize(t, 1e-4, y);

    double ascale = pow(1+x/r,n);
    double kTe_new = pe*v*ascale;
    double kTi_new = (p-pe)*v*ascale;

    bool found_grid=false;
    double kTe_grid=kTe_new, kTi_grid=kTi_new, x_grid=x;
    double target, v_high, v_low, v_mid;
    double v_x=0, v_e=0, v_i=0;
    valarray<double> y_mid(3);
    const double& x_mid = y_mid[0];
    const double& p_mid = y_mid[1];
    const double& pe_mid = y_mid[2];
    double kTe_mid, kTi_mid;
    vector<valarray<double>> grid;
    grid.push_back({v,x,p,pe});

    double kTe_old, kTi_old, x_old, v_old;

    while(t>1e-4 && kTe_new > 0.5*dkTe){

        x_old = x;
        kTe_old = kTe_new;
        kTi_old = kTi_new;
        v_old = v;

        accretion_column.Dense_Step(t, y);

        ascale = pow(1+x/r,n);
        kTe_new = pe*v*ascale;
        kTi_new = (p-pe)*v*ascale;

        if(x_grid-dx < x_old && x_grid-dx >= x){
            found_grid = true;
            target = x_grid-dx;
            if(target == x){
                v_x = v;
            }
            else{
                v_high = v_old;
                v_low = v;
                while(abs(v_high-v_low) > 1e-4){
                    v_mid = (v_high+v_low)/2;
                    accretion_column.Interpolate(v_mid, y_mid);
                    if(x_mid < target){
                        v_low = v_mid;
                    }
                    else{
                        v_high = v_mid;
                    }
                }
                v_x = (v_high+v_low)/2;
            }
        }
        if(((kTe_old-kTe_grid-dkTe)*(kTe_new-kTe_grid-dkTe)<=0)||
           ((kTe_old-kTe_grid+dkTe)*(kTe_new-kTe_grid+dkTe)<=0)){

            found_grid = true;
            target = kTe_grid+dkTe;
            if((kTe_old-kTe_grid-dkTe)*(kTe_new-kTe_grid-dkTe)<=0){
                target = kTe_grid-dkTe;
            }
            if(target == kTe_new){
                v_e = v;
            }
            else{
                v_high = v_old;
                v_low = v;
                if(kTe_old < kTe_new){
                    v_high = v;
                    v_low = v_old;
                }
                while(abs(v_high-v_low) > 1e-4){
                    v_mid = (v_high+v_low)/2;
                    accretion_column.Interpolate(v_mid, y_mid);
                    kTe_mid = pe_mid*v_mid*pow(1+x_mid/r,n);
                    if(kTe_mid < target){
                        v_low = v_mid;
                    }
                    else{
                        v_high = v_mid;
                    }
                }
                v_e = (v_high+v_low)/2;
            }
        }
        if(((kTi_old-kTi_grid-dkTi)*(kTi_new-kTi_grid-dkTi)<=0)||
           ((kTi_old-kTi_grid+dkTi)*(kTi_new-kTi_grid+dkTi)<=0)){

            found_grid = true;
            target = kTi_grid+dkTi;
            if((kTi_old-kTi_grid-dkTi)*(kTi_new-kTi_grid-dkTi)<=0){
                target = kTi_grid-dkTi;
            }
            if(target == kTi_new){
                v_i = v;
            }
            else{
                v_high = v_old;
                v_low = v;
                if(kTi_old < kTi_new){
                    v_high = v;
                    v_low = v_old;
                }
                while(abs(v_high-v_low) > 1e-4){
                    v_mid = (v_high+v_low)/2;
                    accretion_column.Interpolate(v_mid, y_mid);
                    kTi_mid = (p_mid-pe_mid)*v_mid*pow(1+x_mid/r,n);
                    if(kTi_mid < target){
                        v_low = v_mid;
                    }
                    else{
                        v_high = v_mid;
                    }
                }
                v_i = (v_high+v_low)/2;
            }
        }

        if(found_grid){
            found_grid=false;
            v_mid = max(v_x,max(v_e,v_i));
            accretion_column.Interpolate(v_mid, y_mid);
            grid.push_back({v_mid,x_mid,p_mid,pe_mid});
            v_x = 0;
            v_e = 0;
            v_i = 0;
        }
    }

    int n_points = grid.size();
    velocity.resize(n_points);
    altitude.resize(n_points);
    total_pressure.resize(n_points);
    electron_pressure.resize(n_points);
    electron_density.resize(n_points);
    ion_density.resize(n_points);
    electron_temperature.resize(n_points);
    ion_temperature.resize(n_points);
    volume.resize(n_points);

    double mdot;

    volume[0] = pow(1+grid[0][1]/r, n+1);
    for(uint i=0; i<n_points; i++){
        velocity[i] = grid[i][0]*shock_speed;
        altitude[i] = grid[i][1]*shock_height;
        total_pressure[i] = shock_mdot*shock_speed*grid[i][2];
        electron_pressure[i] = shock_mdot*shock_speed*grid[i][3];
        mdot = (accretion_rate/accretion_area)*pow(1+altitude[i]/radius,-n);
        electron_density[i] = (mdot/velocity[i])*density_const/avg_ion_mass;
        ion_density[i] = electron_density[i]/avg_atomic_charge;
        electron_temperature[i] = erg_to_kev*electron_pressure[i]/electron_density[i];
        ion_temperature[i] = erg_to_kev*avg_atomic_charge*(total_pressure[i]-electron_pressure[i])/electron_density[i];
        volume[i] -= pow(1+(grid[i][1]+grid[i+1][1])/(2*r), n+1);
        volume[i] *= accretion_area*radius/(n+1);
        if(i==n_points-1){
            break;
        }
        volume[i+1] = pow(1+(grid[i][1]+grid[i-1][1])/(2*r), n+1);
    }
    volume[n_points-1] = (accretion_area*radius/(n+1))*(pow(1+(grid[n_points-1][1]+grid[n_points-2][1])/(2*r), n+1)-1);

    double dens, gff, cyc, brems, bremss_weight=0;
    cyclotron_ratio = 0;
    for(uint i=0; i<n_points; i++){
        gff = gaunt::gaunt_factor(electron_temperature[i]/erg_to_kev);
        dens = pow((1+r)/(grid[i][1]+r), n)/grid[i][0];
        cyc = (cooling_ratio/gff)*grid[i][3]*grid[i][3]*pow(dens,-3.85)*pow(1+grid[i][1]/r,-8.55-0.425*n);
        brems = volume[i]*gff*sqrt(dens*dens*dens*grid[i][3]);
        cyclotron_ratio += brems*cyc;
        bremss_weight += brems*(1+cyc);
    }
    cyclotron_ratio /= bremss_weight;
}

void Cataclysmic_Variable::Print_Properties(){
    cout << "===================================================" << endl;
    cout << "                   mCV Properties                  " << endl;
    cout << "===================================================" << endl;
    cout << " mass:               " << mass/m_sol << " M_solar" << endl;
    cout << " radius:             " << radius/r_sol << " R_solar" << endl;
    cout << " B_field:            " << b_field/1e6 << " MG" << endl;
    if(inverse_mag_radius != 0){
        cout << " R_m/R:              " << (1./inverse_mag_radius)/radius << endl;
    }
    cout << " accretion rate:     " << accretion_rate << " g/s" << endl;
    cout << " accretion rate:     " << shock_mdot << "-->" << accretion_rate/accretion_area << " g/cm2/s" << endl;
    cout << " shock height:       " << shock_height/radius << " (h/R_wd)" << endl;
    cout << " shock temperature:  " << electron_temperature[0] << " keV" << endl;
    cout << " cooling ratio:      " << cooling_ratio << endl;
    cout << " cycl to brems flux: " << cyclotron_ratio << endl;
}
