#include <cmath>
#include <iostream>
#include <vector>

#include <xsTypes.h>
#include <XSFunctions/functionMap.h>
#include <XSFunctions/funcWrappers.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>
/*
  DEQSPEC
  ----

  Oct 2023: Accretion column temperature/density profile calculator for Intermediate Polars.

*/
const int VGRIDMAX=250;
const int RGRIDMAX=5;
int J,VGRID, RGRID;

std::valarray<Real> RHO(VGRIDMAX);
std::valarray<Real> P(VGRIDMAX);
std::valarray<Real> TK(VGRIDMAX);
std::valarray<Real> X(VGRIDMAX);
std::valarray<Real> NELEC(VGRIDMAX);
std::valarray<Real> SOLN(VGRIDMAX);

Real B, B_approx, B_derived;
Real DISTNORM;
Real R_ratio, M, METABUN;
Real inter1,inter2,inter3,inter4,inter5,inter6;//intermediate/placeholder variables for calculation
//declarations of constants
//Real Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma;
//declaration of white dwarf parameters
Real M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, MDOT, XS0, coeff;
Real R_m;
Real hstar;
Real cosAngle;
Real shockratio;
Real f;
Real L;
Real WDABUN;
Real reflectOn;

Real Msun = 1.989100e+33;//solar mass [gram]
Real Rsun = 6.9599e10;//solar radius [cm]
Real G = 6.672590e-8;//Newton's gravitational constant
Real mue = 1.173;//mean number of electrons
Real mmw = 2.0;//Mean molecular mass
Real mu = 0.615;//mean molecular weight (some papers use mu = 0.5 and Cropper 1999 did comparison)
Real A = 6.99e16;//used by Saxton (this coefficient gives the flux norm close to Suleimanov's model for the same WD mass input - 18% difference)
Real k = 1.380658e-16;//Boltzmann constant
Real mH = 1.672623e-24;//Mass of hydrogen
Real alpha = 2.0;//Thermodynamic Constant
Real beta = 3.85;//Thermodynamic Constant
Real gamma = 5.0/3.0;//Adidabatic Constant

//FUNCTION DECLARATIONs
Real epsilon_func(Real x, const Real coeff);
Real epsilon_func_prime(Real x, const Real coeff);
Real kernel(Real tau);

//Solves for Epsilon_s using Newton's Method
Real Newton_Solver(Real (*func)(Real, Real), Real (*func_prime)(Real, Real), 
                    Real x0, Real coeff, Real tol, int maxIter);

//Accretion Column ODE 
void DEQSPEC_Shooting(int VGRID, RealArray& RHO, RealArray& P, RealArray& TK,
                    RealArray& X, RealArray& NELEC, RealArray& SOLN);

//Shock Height Correction
void XS_CORR(Real h_init, int ITER, Real tol, int VGRID, RealArray& RHO, RealArray& P, RealArray& TK,
                    RealArray& X, RealArray& NELEC, RealArray& SOLN);

//Generator of Spectrum From Accretion Column Profile
void DEQSPEC_MEWE_SPECTRUM(int vgrid, const RealArray& X, const RealArray& TK, const RealArray& NELEC, Real METABUN, 
                            Real DISTNORM, const RealArray& energyArray, int spectrumNumber, RealArray& fluxArray, int NE, const string& initString);

/*---------------------------
        Primary routine
----------------------------*/
extern "C"
void DEQSPEC(const RealArray& energyArray, const RealArray& params, int spectrumNumber, RealArray& fluxArray, RealArray& fluxErrArray,
	  const string& initString)
{
    reflectOn = params[0];
    R_ratio = params[1]; // Magnetospheric radius/WD radius ratio
    f = params[2];//fractional accretion area
    L = params[3]*std::pow(10,33.0); // luminosity
    M = params[4]; // WD mass in solar mass
    METABUN = params[5]; // ACCRETION COLUMN Abundance
    WDABUN = params[6]; // WD SURFACE Abundance
    cosAngle = params[7]; // cos i


    VGRID = 50; //Number of vertical grids fixed to 50
    RGRID = 1;//Number of radial zones in the accretion column

    DISTNORM = 2.62511E-34;

    fluxArray.resize(energyArray.size()-1,0);
    fluxErrArray.resize(energyArray.size()-1,0);
    int NE = fluxArray.size();

    /*for(int i=0;i<NE;i++){
        fluxErrArray[i]=0.0;
    }*/

    //Setting VGRID to max value allowed if it is too large
    if (VGRID>VGRIDMAX){
        std::cerr << "Too many vertical elements, max= " << VGRIDMAX << std::endl;
        VGRID=VGRIDMAX;
    }

    //Computing white dwarf characteristics
    M_3 = (5.816*Msun)/(std::pow(mmw,2.0)); //Chandresekar mass [grams]
    M_wd=M*Msun; //WD mass [grams]
    R_wd=Rsun*(0.0225/mmw)*std::pow(1.0-std::pow(M_wd/M_3,4.0/3.0),0.5)/std::pow(M_wd/M_3,1.0/3.0); //\WD radius [cm]
    R_m=R_wd*R_ratio;

    MDOT = L/(G*M_wd*(1/R_wd-1/R_m));
    MDOT0 = MDOT/(4*3.141592654*std::pow(R_wd,2.0)*f);
    std::cout << "m dot = " << MDOT0 << std::endl;
    //intermediate steps
    inter1 = R_m/((2.75e8)*(std::pow(M,1.0/7.0))*(std::pow((R_wd/(1.e9)),-2.0/7.0)));
    inter2 = (1.e33)/M_wd;
    inter3 = inter2* R_wd / G *(std::pow(inter1,(-7./2.)));
    inter4 = inter3 / (8.6e17*(std::pow(MDOT0,7./5.))*(std::pow(M,1./7.))*(std::pow(R_wd/(1.e9),9./5.)));
    inter5 = (std::pow(inter4,(-1./2.)))/R_wd * (1.e30)/R_wd /R_wd;
    inter6 = inter5/1.e6;

    B = std::pow(inter6,5./7.);
    //std::cout << "B [MG] = " << B << std::endl;

    if(R_wd<R_m){
        vff=std::pow((2.*G*M_wd)*((1./R_wd)-(1./R_m)),0.5);//free fall velocity with R_m correction, Suleimanov (2016)
    }else{
        vff = std::pow((2.*G*M_wd)/(R_wd),0.5);
    }
    n_elec_shock = 4.0*7.013e23*MDOT0/vff;//electron number density
    Tshock = (3*G*M_wd*mu*mH)/(8*k) * (1/R_wd - 1/R_m);//Adjusted shock temperature

    XS0=std::pow(vff,3)*0.049/(2.*A*MDOT0);//shock height [cm] when B = 0
    coeff = std::pow(9.1e-3*(B/10.),2.85)*std::pow((Tshock/1.e8),2.0)*std::pow((n_elec_shock/1.e16),-1.85)*std::pow((XS0/1.e7),-0.85);

    ES0 = Newton_Solver(epsilon_func, epsilon_func_prime, 1.e5, coeff, 1.e-7, 200);

    shock_height = 7.59e6*(4.0/MDOT0)*std::pow(M/0.5,3./2.)*std::pow(R_wd/1.e9,-3./2.)/std::pow(1.+ES0,0.5);//shock height [cm] from equation 7c in\ Wu et al. 1994 paper
    
    int ITER=10;
    DEQSPEC_Shooting(VGRID,RHO,P,TK,X,NELEC,SOLN);

    if ((X[VGRID-1]/R_wd)>0.01 && (X[VGRID-1]/R_wd)<10){
        XS_CORR(X[VGRID-1], ITER, 0.005, VGRID, RHO, P, TK, X, NELEC, SOLN);
    }

    B_approx = 52.2*std::pow(ES0,0.35)*std::pow(Tshock/1.0E8,-0.7)*std::pow(NELEC[VGRID-1]/1.0E16,0.65)*std::pow(shock_height/1.0E7,0.3);
    B_derived = 52.2*std::pow(ES0,0.35)*std::pow(TK[VGRID-1]/1.0E8,-0.7)*std::pow(NELEC[VGRID-1]/1.0E16,0.65)*std::pow(X[VGRID-1]/1.0E7,0.3);

    //std::cout << "White dwarf radius [10^7 cm] = " << R_wd * 1.0e-7 << std::endl;
    //std::cout << "Shock height [10^7 cm] (numerical) = " << X[VGRID-1] * 1.0e-7 << std::endl;
    //std::cout << "Shock height [10^7 cm] (approximate) = " << shock_height * 1.0e-7 << std::endl;  // Same value as numerical
    std::cout << "h/R_wd = " << shock_height/ R_wd << std::endl;
    //std::cout << "Tshock [keV] = " << Tshock * 8.618e-8 << std::endl;
    //std::cout << "B [MG] (numerical) = " << B << std::endl;
    //std::cout << "B [MG] (approximate) = " << B_approx << std::endl;
    //std::cout << "B [MG] (derived) = " << B_derived << std::endl;
    //std::cout << "V_ff [10^8 cm/s] = " << vff*1.0e-8 << std::endl;

    DEQSPEC_MEWE_SPECTRUM(VGRID, X, TK, NELEC, METABUN, DISTNORM, energyArray, spectrumNumber, fluxArray, NE, initString);
}

//function inside the integral in Wu 1994
Real kernel(Real tau){
    Real constant;
    Real ker;
    constant = std::pow(vff,3)/(2.*A*MDOT0);
    ker=constant*( (std::pow(tau,1.5)*(5.-8.*tau))/(std::pow(1-tau,0.5))/
         ((1+(std::pow(3.,-alpha))*(std::pow(4,alpha+beta))*(ES0)*(std::pow(1.-tau,alpha)*(std::pow(tau,beta))))));
    return ker;
}

//function for f(x) = 0 (from Wu et al. 1994 and 1995 papers). Using equation (10) in Wu94 paper and xs/xs0 = (1+es0)^(-1/2) and we solve for es0 for a given B-field.
Real epsilon_func(const Real x, const Real coeff)
{
    return coeff*std::pow(1.+x,0.425) - x;
}
//1st derivative of the above function
Real epsilon_func_prime(const Real x, const Real coeff)
{
    return 0.425*coeff*std::pow(1.+x,-0.575)-1.0;
}

/*
Solving f(x) = 0 using the Newton's method

epsilon_func: f(x)
epsilon_func_prime: 1st derivative of f(x)
start: initial guess for x
epsilon: tolerance
max_iter: maximum number of iterations
root: solution for x 
*/
struct Params {
    Real coeff;
};
Real func_wrapper(Real x, void *params) {
    Params *p = static_cast<Params*>(params);
    return epsilon_func(x, p->coeff);
}

Real func_prime_wrapper(Real x, void *params) {
    Params *p = static_cast<Params*>(params);
    return epsilon_func_prime(x, p->coeff);
}

void func_and_prime_wrapper(Real x, void *params, Real *y, Real *dy) {
    Params *p = static_cast<Params*>(params);
    *y = epsilon_func(x, p->coeff);
    *dy = epsilon_func_prime(x, p->coeff);
}

Real Newton_Solver(Real (*func)(Real, Real), Real (*func_prime)(Real, Real), 
                   Real x0, Real coeff, Real tol, int maxIter) {
    Params params = {coeff};

    gsl_function_fdf F;
    F.f = &func_wrapper;
    F.df = &func_prime_wrapper;
    F.fdf = &func_and_prime_wrapper;
    F.params = &params;

    gsl_root_fdfsolver *solver = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_newton);
    gsl_root_fdfsolver_set(solver, &F, x0);

    Real x = x0;
    int iter = 0;
    int status;

    do {
        iter++;
        status = gsl_root_fdfsolver_iterate(solver);
        x = gsl_root_fdfsolver_root(solver);
        status = gsl_root_test_residual(func_wrapper(x, &params), tol);
    } while (status == GSL_CONTINUE && iter < maxIter);

    gsl_root_fdfsolver_free(solver);

    if (status != GSL_SUCCESS) {
        std::cerr << "Failed to converge after " << maxIter << " iterations." << std::endl;
    }

    return x;
}

//method to calculate temperature/density profile in accretion column
void DEQSPEC_Shooting(int VGRID, RealArray& RHO, RealArray& P, RealArray& TK,
                    RealArray& X, RealArray& NELEC, RealArray& SOLN){
    Real tau, y, tau_f, m1, m2, m3, m4, tau_i, y_i, n;
    int number_of_steps; //counter variables

    std::valarray<Real> vfinal(VGRID + 1), taufinal(VGRID);
    std::valarray<Real> TkeV(VGRID);

    tau_i = 0;//lower bound of integration
    y_i=0;//first test value in the RK4 method
    tau_f=0.25; //final integration value
    number_of_steps=VGRID; //number of integration steps
    n=(tau_f-tau_i)/static_cast<float>(VGRID);//integration step size

    /*--------------------------------------------------------
                            RK4 Method
    --------------------------------------------------------*/

    tau=tau_i;
    y=y_i;

    for (int steps = 1; steps <= VGRID; steps++) {
        m1 = kernel(tau);
        m2 = kernel(tau + 0.5 * n);
        m3 = kernel(tau + 0.5 * n);
        m4 = kernel(tau + n);

        tau = tau + n;
        y = y + n * (m1 + 2.0 * m2 + 2.0 * m3 + m4)/6.0;

        // These three statements fill the soln, taufinal, and vfinal grid/array.
        taufinal[steps - 1] = tau_i + (tau_f - tau_i) * static_cast<Real>(steps) / VGRID;//normalized velocity grid
        SOLN[steps - 1] = y;//non-normalized position grid
        vfinal[steps - 1] = taufinal[steps - 1]*vff;
    }

    for (int r = 0; r < VGRID; r++) {
        X[r] = SOLN[r] - SOLN[0]; // position grid
        
        RHO[r] = MDOT0 / vfinal[r]; // density grid
        // P(r) from Wu 1994 paper (KM)
        P[r] = MDOT0 * vff * (1.0 - taufinal[r]);// pressure grid

        // TK(r) is the actual temperature grid
        TK[r] = (16.0 * taufinal[r] * (1.0 - taufinal[r]) * (1.0 / 3.0)) * Tshock; // actual temperature grid

        TkeV[r] = TK[r] * 8.6173e-8; // the actual temperature grid in keV
        NELEC[r] = RHO[r] * 7.01e23; // electron density grid
    }
}


void DEQSPEC_MEWE_SPECTRUM(int VGRID, const RealArray& X, const RealArray& TK, const RealArray& NELEC, Real METABUN, 
                            Real DISTNORM, const RealArray& energyArray, int spectrumNumber, RealArray& fluxArray, int NE, const string& initString){
    Real PI, KK, NENH;
    std::valarray<Real> PARAM1(3);
    std::valarray<Real> flux(NE), fluxError(NE);
    PI=3.141592654;
    KK=11.6048e6;
    NENH=1.21;
    for(int j=0;j<NE;j++){
        fluxArray[j]=0.0;
    }
        
    std::valarray<Real> REFLECTPARAM1(5);
    shockratio = X[VGRID-1] / R_wd;
    if (reflectOn==1){
        REFLECTPARAM1[0]=1-std::pow(1.-1./std::pow(1+shockratio,2),0.5);
        //"reflect param = " << REFLECTPARAM1[0] << std::endl;
        REFLECTPARAM1[1]=0.0;
        REFLECTPARAM1[2]=WDABUN;
        REFLECTPARAM1[3]=WDABUN;
        REFLECTPARAM1[4]=cosAngle;
    }
    //Main loop. Loop over each vertical element to calculate and add the flux into the array
    for(int k=0;k<VGRID;k++){
        
        /* Calculates the Mewe spectrum for each vertical element on the energy
        grid passed into it, using
        a) PARAM(1) the temperature in keV of the vertical element
        b) PARAM(2) the density of electrons nH (cm^-3) of the vertical eleme
        c) PARAM(3) the heavy metal abundance
        d) PARAM(4) is the redshift, here set to zero

        MEKAL parameter set
                PARAM1(1) = TK(J)/KK
                PARAM1(2) = NELEC(J)
                PARAM1(3) = METABUN
                PARAM1(4) = 0.0
                PARAM1(5) = 0.0
                PARAM1(6) = 0.0

        APEC parameter set
        */
        PARAM1[0]=TK[k]/KK;
        PARAM1[1]=METABUN;
        PARAM1[2]=0.0;
        if(TK[k]/KK<86){
            apec(energyArray, PARAM1, spectrumNumber,flux, fluxError, initString);
        }else{
            //std::cout << "is this ever used? " << std::endl;
            //xsbrms_((float*)&energyArray[0], NE, (float*)&PARAM1[0], spectrumNumber, (float*)&flux[0], (float*)&fluxError[0]);
            CXX_bremss(energyArray, PARAM1, spectrumNumber, fluxArray, fluxError,initString);
        }

        /*
     ! Now multiplies the calculated spectrum by the volume and density^2
     ! of the particular element for each energy bin. This is required becaus
     ! the XSMEKL code has the density stripped out of the flux it returns,
     ! expecting this to be applied in the XSPEC normalisation. The 10^14 ari
     ! because of the units of 10^50 cm^3 at a distance of 1 pc = 3.086x10^18
     ! (ie 10^50/(10^(18^2))).
     !
     ! Explicitly we have:
     !  lumin(fmekal) = flux(fmekal)*(4pi*(3.086^18)^2)
     ! so for 1 cm^3 of gas
     !  lumin(fmekal) = flux(fmekal)*(4pi*(3.086^18)^2)/10^50
     ! as this passes to XSMEKL through XSVMKL, flux(fmekal) gets multiplied
     ! by 4pi*3.086^2/Ne^2 for XSPEC normalisations, so this needs to be reve
     ! so that to get the correct luminosity from XSMEKL we must have
     !  lumin(xsmekl) = flux(fmekal)*[(4pi*(3.086^18)^2)/10^50]*[Ne^2/(4pi*3.
     !                = flux(fmekal)*Ne^2*10^-14
     ! and then to get the fluxes from a volume element dV we then have
     !  flux(xsmekl)  = lumin(xsmekl)*dV/(4pi*dist^2)   (dist in cm)
     !                = flux(fmekal)*Ne^2*(10^-14)*dV/(4pi*dist^2)
     !
     ! On 1 June 1998 this routine was changed to reflect the fact that the
     ! XSVMKL normalisation is acutally 4pi*3.086^2/(NeNh)
     !                                = 4pi*3.086^2/(Ne^2/1.21)
     ! rather than 4pi*3.086^2/(Ne^2)
     !
     !DIST = (DISTNORM*3.085678E18)**2
*/
        for(int i=0;i<NE;i++){
            if(k==0){
                flux[i]=DISTNORM*X[k]*(std::pow(NELEC[k]*1e-7,2)/NENH)*flux[i];
            }else{
                flux[i]=DISTNORM*(X[k]-X[k-1])*(std::pow(NELEC[k]*1e-7,2)/NENH)*flux[i];
            }
        }
        /* adds the flux for this vertical element to the existing flux in each
     ! energy bin
     !  */
        for(int l=0;l<NE;l++){
            fluxArray[l]=fluxArray[l]+flux[l];
        }
    }
    if (reflectOn==1){
    CXX_reflect(energyArray, REFLECTPARAM1, spectrumNumber, fluxArray, fluxError, initString);
    }
}

void XS_CORR(Real h_init, int ITER, Real tol, int VGRID, RealArray& RHO, RealArray& P, RealArray& TK,
                    RealArray& X, RealArray& NELEC, RealArray& SOLN){
    Real hstar;
    Real vff_old,h_old;

    hstar=h_init;
    h_old=h_init;
    vff_old = vff;
    
    for (int i=0;i<ITER;++i){
        vff = std::pow((2.*G*M_wd)*((1./(R_wd+hstar))-(1./R_m)),(1./2.));
        n_elec_shock=4.0*7.013e23*MDOT0/vff;
        Tshock = (3*G*M_wd*mu*mH)/(8*k) * (1/(R_wd+hstar) - 1/R_m);
        XS0 = std::pow(vff,3)*0.049/(2.*A*MDOT0); //shock height [cm] when B = 0
        coeff = 9.1e-3*std::pow(B/10.,2.85)*std::pow(Tshock/1.e8,2.0)*std::pow(n_elec_shock/1.e16,-1.85)*std::pow(XS0/1.e7,-0.85);
        ES0 = Newton_Solver(epsilon_func,epsilon_func_prime,1.e5,coeff,1.e-6,200);

        DEQSPEC_Shooting(VGRID,RHO,P,TK,X,NELEC,SOLN);
        hstar = X[VGRID-1];
        //std::cout << "shock_height " << X[VGRID - 1] << std::endl;

        if (R_wd + hstar > R_m) {
            shock_height = h_old;
            std::cout << "ACCRETION COLUMN IS TOO TALL" << std::endl;
            return;
        }

        if (std::abs((h_old - hstar) / h_old) <= tol) {
            return;
        }

        h_old=hstar;
        vff_old=vff;
    }
}