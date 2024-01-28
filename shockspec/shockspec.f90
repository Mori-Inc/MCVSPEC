!                                                                       
!  SHOCKSPEC                                                                  
!  ----                                                                 
!                                                                       
!  Feb 2023: Accretion column temperature/density profile calculator for magnetic CVs with x_s correction.        
!                                                                       
!  INPUT PARAMETERS: 
!  PARAM(1) -  
!  PARAM(2) - 
!  PARAM(3) -
!  PARAM(4) - 
!  PARAM(5) - 
!  PARAM(6) - 
!  PARAM(7) - 


SUBROUTINE SHOCKSPEC(EAR,NE,PARAM,IFL,PHOTAR,PHOTER) 

  IMPLICIT NONE 
                                                                        
  INTEGER    NE,IFL,ITER
  
  INTEGER    VGRIDMAX,RGRIDMAX 
  
  PARAMETER  (VGRIDMAX=250) 
  PARAMETER  (RGRIDMAX=5) 
  
  REAL     EAR(0:NE),PARAM(*),PHOTAR(NE),PHOTER(NE) 
  REAL     METABUN,ANG 
  
!  REAL     ES0,MDOT0
!  REAL     RC     
!  REAL     ES(0:RGRIDMAX),MDOT(0:RGRIDMAX),AREA(0:RGRIDMAX) 
  REAL     M, B, B_approx,B_derived!,AREA0 
  REAL     DISTNORM 
  REAL     RHO(VGRIDMAX),P(VGRIDMAX),TK(VGRIDMAX),                &
       &               X(VGRIDMAX),NELEC(VGRIDMAX), SOLN(VGRIDMAX)                        
  !REAL     EMIN,EMAX,CONSTANT 
  INTEGER    J,VGRID,RGRID 
  REAL, EXTERNAL :: func, func_prime
  
  REAL h_star

  REAL      Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma
  REAL M_3 , M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, XS0, coeff                                         
  common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, XS0, coeff, B       
  common /constants/ Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma
   
  Msun     = 1.989100e+33 ! solar mass [gram]                                                                                   
  Rsun     = 6.9599e10 ! solar radius [cm]                                                                                      
  G        = 6.672590e-8 ! Gravitational constant                                                                               
  mue      = 1.173 ! - mean number of electrons                                                                                 
  mmw      = 2.0    !- Mean molecular mass-                                                                                                                    
  mu       = 0.615  ! - mean molecular weight (some papers use mu = 0.5 and Cropper 1999 did comparison)
  !A        = 3.9e16 ! Wu et al. 1994 paper used this value (for bremmstraulung coefficient)    
  A        = 4.91e16 ! used by Cropper's code instead of Wu's 3.9e16 (this coefficient gives the flux norm close to Suleimanov's model for the same WD mass input - 18% difference). 
  k        = 1.380658e-16 ! -Boltzmann constant-                                                                                
  mH       = 1.672623e-24 ! -Mass of hydrogen-                                                                                  
  alpha    = 2.0 ! -Thermodynamic Constant-                                                                                     
  beta     = 3.85 ! -Thermodynamic Constant-                                                                                    
  gamma    = 5.0/3.0  ! -Adidabatic Constant- 


!                                                                       
!                                                                       
! RGRID is the number of radial zones in the accretion column.          
!                                                                       
  RGRID = 1 
!                                                                       
! These are the input parameters:                                       
!  Mdot(0) is the specific accretion rate on the central axis of the col
!  Es(o) is the ratio of cyclotron to brems cooling at the shock front  
!  Rc is the radius of the column in 10^9 cm                            
!  M is the mass of the primary in solar masses                         
!  METABUN is the metal abundance: solar = 1.0                          
!  VGRID is the number of vertical layers in the accretion column       
!  ANG is the viewing angle used to modify the albedo in degrees        
!                                                                       
  B     = PARAM(1) ! Surface B-field [10^6 G] 
  MDOT0  = PARAM(2) ! mass accretion rate [g/cm2/s] 
!  RC      = PARAM(3) 
  M       = PARAM(3) ! WD mass in solar mass 
  METABUN = PARAM(4) ! Abundance 
  VGRID   = 50 ! Number of vertical grids fixed to 50 
  ANG     = 1 ! PARAM(7) unused parameter 
  !                                                                       
! Flux normalization factor for APEC model as we define our flux norm as (R [km] / d [kpc])^2
! Note: PI*R_C^2/(4*PI*Distance^2) where R_C = accretion column radius [cm] and distance [cm]                                
!                                                                       
  DISTNORM = 2.62511E-34 

!  DISTNORM = 1.0 
!                                                                       
!        IF(J.EQ.RGRID) THEN                                            
!           DO L=1,NE                                                   
!              PRINT '(1X,I5,3(1PE12.4))',L,EAR(L-1),EAR(L),PHOTAR(L)   
!           ENDDO                                                       
!        ENDIF                                                          
                                                                        
! Change units of RC from 10^9 cm to cm                                 
!                                                                       
!  RC = RC*1.0E9 
!  RC = 1.0e9 ! Note: RC used to be one of the input parameters (PARAM(3)) but we should fold in RC into the flux normalization like PI*RC^2/(4*PI*Distance^2) just like Suleimanov. 
!                                                                       
! Zero error array                                                      
!                                                                       
  DO J=1,NE 
     PHOTER(J)=0.0 
  ENDDO
  !                                                                       
  
! **This checks to make sure that there aren't too many vertical layers 
! Checks that VGRID is not too large                                    
!                                                                       
  IF(VGRID.GT.VGRIDMAX)THEN 
     PRINT *,'Too many vertical elements, max = ',VGRIDMAX 
     VGRID = VGRIDMAX 
  ENDIF
  

!           Calculating Neccesary WD Characteristics                                             
  M_3  = (5.816*Msun)/(mmw**2.0) ! -Chandresekar Mass [gram]  
  M_wd = M*Msun ! WD mass [grams] 
  R_wd = Rsun*(0.0225/mmw)*SQRT(1.0-(M_wd/M_3)**(4.0/3.0))/((M_wd/M_3)**(1.0/3.0))  ! WD radius [cm]                                                                                
  vff = ((2.*G*M_wd)/(R_wd+shock_height*0.0))**(1./2.) ! -Free Fall Velocity- at the shock height (we should add shock_height to R_wd)              
  n_elec_shock = 4.0*7.013e23*MDOT0/vff !4.0*MDOT0*mue/vff/mH ! electron number density at shock                  
  !      constant = (1.55e8)*((c/4.)**(-1))*((M_wd/(0.5*Msun))**1.5)*((Rwd/1.e9)**(-1.5)) ! Mult. Constant in Wu 1994-                                                                       
  TShock = (3*G*M_wd*mu*mH)/(8*k*(R_wd+shock_height*0.0)) ! shock temperature [K] from Yu et al. 2018 paper formula (1)                                                                                  
!  AREA0 = 3.14159*RC**2   ! the area of the accretion column [cm^2]  

  XS0 = vff**3*0.049/(2.*A*MDOT0) ! shock height [cm] when B = 0 

  coeff = 9.1e-3*(B/10.)**2.85*(TShock/1.e8)**2.0*(n_elec_shock/1.e16)**(-1.85)*(XS0/1.e7)**(-0.85) 

  CALL Newton_solver(func, func_prime, 1.e9, 1.e-7, 100, ES0) ! solving for ES0 

  shock_height = 7.59e6*(MDOT0/4.0)**(-1.)*(M/0.5)**(3./2.)*(R_wd/1.e9)**(-3./2.)/(1.+ES0)**0.5 ! shock height [cm] from equation 7c in\ Wu et al. 1994 paper 

!           Checking WD characteristics calculation                                              
!  write(*,*) M_3, M_wd, R_wd, shock_height, vff, MDOT0 


!                                                                       
! Now generates the radial mass transfer rate and Es. If RGRID=1 then   
! there is no radial dependence.                                        
!                                                                       
      !CALL SHOCKSPEC_HORIZ_GENERATE(RGRID,ES0,MDOT0,RC,                    &
       !    &                            ES,MDOT,AREA)                         
  
!  DO J=1,RGRID 
     !        write(*,*) VGRID,ES0, M,MDOT(J), RHO, P, TK, X, NELEC, SOLN 
  CALL SHOCKSPEC_SHOOTING(VGRID,RHO,P,TK,X,NELEC,SOLN)          
  ITER = 10

  IF(((X(VGRID)/R_wd).GT.(0.03)).AND.((X(VGRID)/R_wd).LT.(2))) THEN
    !CALL Shock_Shooting(shock_height, h_star, ITER, .005, func, func_prime)
    CALL SHOCK_SHOOT(X(VGRID), h_star, ITER, .01, func, func_prime, VGRID, RHO, P, TK, X, NELEC,SOLN)
  ENDIF  
!                                                                       
! Calculate the magnetic field. This is done using eqn(10) of Wu, Chanmu
! and Shaviv (1994) ApJ 426, 664, inverted to determine B from the other
! parameters. KM note: I think B is in unit of MegaGauss.                                                             
!                                                                       
  B_approx = 52.2*(ES0**0.35)*(TShock/1.0E8)**(-0.7)*(NELEC(VGRID)/1.0E16)**0.65*(shock_height/1.E7)**0.3 

  B_derived = 52.2*ES0**0.35*                                          &
       &             (1.0E8/TK(VGRID))**0.7*                              &
       &             (NELEC(VGRID)/1.0E16)**0.65*                         &
       &             (X(VGRID)/1.0E7)**0.3 
  PRINT '(A, F10.2)', 'White dwarf radius [10^7 cm] = ', R_wd*1e-7
  !PRINT '(A, F10.2)', 'Shock height [10^7 cm] (numerical) = ', X(VGRID)*1e-7
  PRINT '(A, F10.2)', 'Shock height [10^7 cm] (approximate) = ', shock_height*1e-7 
  !PRINT '(A, F10.2)', 'Tshock [keV] = ', TShock*8.618e-8 
  !PRINT '(A, F10.2)', 'B [MG] (numerical) = ', B_derived 
  !PRINT '(A, F10.2)', 'B [MG] (approximate) = ', B_approx 
  !PRINT '(A, F10.2)', 'V_ff [10^8 cm/s] = ', vff*1e-8  
  !PRINT '(A, F10.2)', 'N_e [10^14 cm^-3] = ', NELEC(VGRID)*1e-14
  !PRINT '(A, F10.2)', 'N_e [10^14 cm^-3] = ', n_elec_shock*1e-14
  write(*,*) ES0, TK(VGRID), TSHOCK, NELEC(VGRID), X(VGRID), shock_height  
  !                                                                       
  CALL SHOCKSPEC_MEWE_SPECTRUM(VGRID,X,TK,NELEC,METABUN,    &
       &              DISTNORM,EAR,NE,IFL,PHOTAR)                     
     
     
  !ENDDO
  
  RETURN 
END SUBROUTINE SHOCKSPEC

! Shooting method to calculate temperature/density profile in accretion column 

subroutine SHOCKSPEC_SHOOTING(VGRID,RHO,P,TK,X,NELEC,SOLN) 

!           Declaration of RK4 Variables
  REAL   tau, y, tau_f, m1, m2, m3, m4, f, tau_i, y_i, n
  
  !           Declaring Common Block Variables - this passes them to the function f 

  common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, XS0, coeff, B    
  common /constants/ Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma
  
  !           Declaration of counter variables
  INTEGER steps, r, number_of_steps, VGRID  
  
  !           Declaration of arrays
  REAL    vfinal(VGRID+1), soln(VGRID), taufinal(VGRID), MDOT0 
  REAL    RHO(VGRID),P(VGRID),TK(VGRID),X(VGRID),NELEC(VGRID),TkeV(VGRID) 
  
  !           Calculating Integration Paramters

  tau_i    = 0 ! -This is the lower bound of integration-
  y_i      = 0 ! -This is the first test value in the RK4 test method-
  tau_f    = 0.25 ! -This is the final integration value-
  number_of_steps = VGRID ! -This is the number of integration steps-
  n = ((tau_f - tau_i)/float(VGRID)) ! -This is the integration step size-
  
  !------------------------------------------------------------------------------
  !                                  RK4 Method
  !------------------------------------------------------------------------------
  
  
  !           Assiging Initial Variable Values            
  tau = tau_i
  y = y_i
 
!  write(*,*) MDOT0, ES0 
 
  !           Performing RK4 Method
  do steps = 1,VGRID 
     
     m1  = n*f(tau)
     m2  = n*f(tau+0.5*n)
     m3  = n*f(tau+0.5*n)
     m4  = n*f(tau)
     
     tau   = tau+n
     y   = y+(m1+2.0*m2+2.0*m3+m4)*(1/6.0)
 
    
     ! These three statements fill the soln, tau final, and vfinal grid/array so that we can
     ! calculate the rest of the properties based on the initial conditions.
     
     taufinal(steps) = (tau_i) + ((tau_f)-(tau_i)) * float(steps)/(VGRID) !* float(steps)/float(h) ! -This is the normalized velocity grid- CHECK THIS LINE! 
     soln(steps)     = y ! -This is the non-normalized position grid-
     vfinal(steps)   = taufinal(steps)*vff 
!     vfinal(steps)   = (tau_i*vff) + ((tau_f*vff)-(tau_i*vff)) * float(steps)/float(number_of_steps) ! -This is the non-normalized velocity grid-
     
!     write(*,*) soln(steps), vfinal(steps) 

    
     !           Checking the values for the calulated constants and grids
     !write(*,*),'DEBUG Part3 ',steps, vff, Rwd, M_wd, G, constant, ES_IN, n
     !PRINT '(1X,I5,4(1PE12.3))', steps, taufinal(steps), soln(steps), vfinal(steps) 
     
  end do
  
  ! These statement fill their respective arrays with the formulas prescribed in Wu 1994.
  
  do r = 1, VGRID
     
     X(r)         = soln(r)-soln(1) ! -This is the position grid-
     RHO(r)       = MDOT0/vfinal(r) ! -This is the density grid-
     !P(r)         = MDOT0*(vff**2.)*(taufinal(r)-soln(r)) ! -This is the pressure grid-
     P(r)         = MDOT0*vff*(1.-taufinal(r)) ! from Wu 1994 paper (KM) 
!     P(r)         = MDOT0*(soln(r) - taufinal(r)) ! Gavin's original version. I don't understand this formula (KM)  
     TK(r)        = (16.*(taufinal(r))*(1. - taufinal(r))*(1.0/3.0)) * TShock ! -This is the actual temperature grid-
     !TK(r)        = (mu*mH/k)*P(r)/RHO(r) ! Gavin's original code uses pressure and density 
     
     TkeV(r)      = TK(r)*(8.6173e-8) ! -This is the actual temperature grid in keV-
     NELEC(r)     = (RHO(r)*7.01e23)! 5.09687e23)!mue/mH)  ! electron density grid. I multiplied the actual number otherwise dividing by mH*mue outputs weird numbers (due to not good enough precision?). I think Wu 1994 assumed mue = 1.0 (all ionized hydrogen gas) I think Gavin's code got this wrong by multiplying 1/(mH*mue).     
     
     !              Checking the respective arrays
!     write(*,*) r, taufinal(r),  vfinal(r), MDOT0, X(r), RHO(r), P(r), TkeV(r), 6.42e-16*P(r)/RHO(r), NELEC(r) 
     
  end do
  
END subroutine SHOCKSPEC_SHOOTING

! This is the function inside the integral in Wu 1994
real function f(tau)
  
  !    common  /xsp/ ES_IN,c,Rwd! ,alpha,beta,gamma,G,A ! common block variables 
  common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, XS0, coeff, B   
  common /constants/ Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma 
  
  REAL constant, MDOT0  ! we need to declare MDOT0 again here otherwise MDOT0 returns a weird number due to type mismatch? 
  !constant = (1.55e8)*((MDOT0/4.)**(-1))*((M_wd/(0.5*Msun))**1.5)*((R_wd/1.e9)**(-1.5))
  
  constant = vff**3/(2.*A*MDOT0) ! replaced by KM (somehow the above constant using MDOT, M_wd, R_wd etc. didn't match with Wu results well

!  write(*,*) MDOT0, M_wd/(0.5*Msun)  

  f = constant*((((tau**1.5)*(5.-8.*tau))/(((1-tau)**0.5))/ &
         ((1+(3.**(-alpha))*(4**(alpha+beta))*(ES0)*((1.-tau)**alpha)*(tau**beta)))))
 
!  write(*,*) f 
 
end function f
                         
                                                                        
! --------------------------------------------------------------------  

! We are not using this subroutine as we don't consider any radial dependence. 
                                                                        
SUBROUTINE SHOCKSPEC_HORIZ_GENERATE(RGRID,ES0,MDOT0,RC,              &
     &                          ES,MDOT,AREA)                           
  
  IMPLICIT NONE 
  
  INTEGER    J,RGRID 
  REAL     ES0,MDOT0,RC 
  REAL     ES(0:RGRID),MDOT(0:RGRID),AREA(0:RGRID) 
  REAL     DR,PI 
  
  PI = 3.141592654 
  DR = RC/FLOAT(RGRID) 
  !                                                                       
  ! The prescription used for the radial dependence is that in Wu (1994)  
! Proc ASA 11(1) 61-64. The dependence is azimuthally symmetric and the 
! mass transfer rate at the edge approaches 1/4 of that at the centre as
! RGRID becomes large. If there is only one radial grid point (RGRID=1) 
! then it is naturally all at the central mass transfer rate.           
!                                                                       
  DO 100 J=1,RGRID 
     AREA(J) = PI*((DR*J)**2-(DR*(J-1))**2) 
     MDOT(J) = MDOT0*(RC/(RC+DR*(J-1)))**2 
     ES(J)   = ES0*((RC+DR*(J-1))/RC)**(3.7) 
     !PRINT '(1X,I5,4(1PE12.3))',J,DR*J,AREA(J),MDOT(J),ES(J)        
100 END DO
  
  RETURN 
END SUBROUTINE SHOCKSPEC_HORIZ_GENERATE
!                                                                             
!                                                                                                                          
! ====================================================================  
!  end of shooting subroutines...                                       
! --------------------------------------------------------------------  
                                                                        
SUBROUTINE SHOCKSPEC_MEWE_SPECTRUM(VGRID,X,TK,NELEC,METABUN,    &
     &                        DISTNORM,EAR,NE,IFL,PHOTAR)           
  
  IMPLICIT NONE 
  
  INTEGER    IFL,NE 
  INTEGER    VGRID,J,L 
  REAL     EAR(0:NE),PHOTAR(NE)!,ANG!,FACT 
  REAL     PARAM1(3) 
  REAL     METABUN,FLX(NE),FLXERR(NE) 
  REAL     TK(VGRID),X(VGRID),NELEC(VGRID) 
  REAL     KK,PI,NENH 
  REAL     DISTNORM!,DIST,NORM_FACTOR 
  
!  REAL     GAMMA,THETA,ALPHA,ALB 
  
  PI = 3.141592654 
  !                                                                       
  ! Constant to convert T(Kelvin) to kT                                   
  !                                                                       
  KK  = 11.6048E6 
  !                                                                       
  ! Constant to convert Ne to Nh                                          
  !                                                                       
  NENH = 1.21 
  !                                                                       
  ! First zero the flux array                                             
  !                                                                       
  DO L=1,NE 
     PHOTAR(L) = 0.0 
  ENDDO
  !                                                                       
  ! Main loop to label 100.                                               
  ! Loop over each vertical element to calculate and add the flux into the
  ! array                                                                 
  !                                                                       
  DO 100 J=1,VGRID 
     !                                                                       
     ! Calculates the Mewe spectrum for each vertical element on the energy  
     ! grid passed into it, using                                            
     !  a) PARAM(1) the temperature in keV of the vertical element           
     !  b) PARAM(2) the density of electrons nH (cm^-3) of the vertical eleme
     !  c) PARAM(3) the heavy metal abundance                                
     !  d) PARAM(4) is the redshift, here set to zero                        
     !                    
     
     ! MEKAL parameter set                                                    
     !       PARAM1(1) = TK(J)/KK 
     !       PARAM1(2) = NELEC(J) 
     !       PARAM1(3) = METABUN 
     !       PARAM1(4) = 0.0 
     !       PARAM1(5) = 0.0 
     !       PARAM1(6) = 0.0 
     
     ! APEC parameter set 
     PARAM1(1) = TK(J)/KK
     PARAM1(2) = METABUN 
     PARAM1(3) = 0.0 
     
     !        PRINT '(1x,i4,5(1pe12.4))',                                    
     !     :              J,PARAM1(1),PARAM1(2),PARAM1(3),PARAM1(4),X(J)     
     !CALL XSMEKL(EAR,NE,PARAM1,IFL,FLX,FLXERR) ! for older versions of HEASOFT 
     !       CALL MEKAL(EAR,NE,PARAM1,IFL,FLX,FLXERR) ! compatible with Mekal model in more recent HEASOFT 
     
     IF(TK(J)/KK.LT.86.)THEN
        CALL APEC(EAR,NE,PARAM1,IFL,FLX,FLXERR)
!        CALL BREMSS(EAR,NE,PARAM1,IFL,FLX,FLXERR)
     !      PRINT '(1x,I5,5(1pe12.3))',J,(FLX(L),L=1,5)                     
     ELSE
        CALL BREMSS(EAR,NE,PARAM1,IFL,FLX,FLXERR) ! added by KM if kT > 86 keV as it exceeds APEC's temperature upper limit. 
     ENDIF
     
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
     
!     DIST=1. 
!     AREA=1.
!     NORM_FACTOR = 2.62511E-34 ! Normalization factor adjustment as we use (R [km] / d [kpc])^2 as our XSPEC flux norm 

     DO L=1,NE 
        IF(J.EQ.1)THEN 
            FLX(L) = DISTNORM*X(J)*((NELEC(J)**2)/NENH)*1.0E-14*FLX(L)   
!           FLX(L) = (AREA*X(J)         /(4.0*PI*DIST))*              &
!                &                              ((NELEC(J)**2)/NENH)*1.0E-14*FLX(L) 
        ELSE 
           FLX(L) = DISTNORM*(X(J)-X(J-1))*((NELEC(J)**2)/NENH)*1.0E-14*FLX(L)
!           FLX(L) = (AREA*(X(J)-X(J-1))/(4.0*PI*DIST))*              &
!                &                              ((NELEC(J)**2)/NENH)*1.0E-14*FLX(L) 
        ENDIF

!        write(*,*) AREA, X(J), NELEC(J), NENH, FLX(L)
     ENDDO
     !                                                                       
     ! adds the flux for this vertical element to the existing flux in each  
     ! energy bin                                                            
     !                                                                       
     DO L=1,NE 
        PHOTAR(L) = PHOTAR(L) + FLX(L) 
     ENDDO
     
100 END DO

  RETURN 
END SUBROUTINE SHOCKSPEC_MEWE_SPECTRUM


! function for f(x) = 0 (from Wu et al. 1994 and 1995 papers). Using equation (10) in Wu94 paper and xs/xs0 = (1+es0)^(-1/2) and we solve for es0 for a given B-field. 

 REAL FUNCTION func(x)
!   IMPLICIT NONE 
   common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, XS0, coeff, B
  
   REAL, INTENT(IN) :: x!, coeff  
   func = coeff*(1.+x)**0.425-x
 END FUNCTION func

! 1st derivative of the above function 

 REAL FUNCTION func_prime(x)
 !  IMPLICIT NONE 
   common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, XS0, coeff,B 
   REAL, INTENT(IN) :: x!, coeff  
   func_prime = 0.425*coeff*(1.+x)**(-0.575)-1.0
 END FUNCTION func_prime


! Solving f(x) = 0 using the Newton's method                                                                                    
!                                                                                                                         
! func: f(x) 
! func_prime: 1st derivative of f(x) 
! start: initial guess for x 
! epsilon: tolerance 
! max_iter: maximum number of iterations 
! root: solution for x 
  
SUBROUTINE Newton_solver(func, func_prime, start, epsilon, max_iter, root)

   IMPLICIT NONE

   REAL, EXTERNAL :: func, func_prime
   REAL, INTENT(IN) :: start, epsilon
   INTEGER, INTENT(IN) :: max_iter

   REAL, INTENT(INOUT) :: root

   INTEGER :: i
   REAL :: f_val, f_der

   ! Begin the iteration up to the maximum number specified                                                                      
   root = start

   DO i = 1, max_iter
      f_val = func(root)

      IF(ABS(f_val ) <= epsilon) THEN
          RETURN
      END IF

      f_der = func_prime(root)

      root = root - f_val/f_der

   END DO

END SUBROUTINE Newton_solver

SUBROUTINE SHOCK_SHOOT(h, hstar, max_iter, tol, func, func_prime, VGRID,RHO,P,TK,X,NELEC,SOLN)
  REAL h, hstar, tol
  INTEGER  max_iter, VGRID
  INTEGER i
  REAL, EXTERNAL :: func, func_prime


  common  /wdparams/ M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, XS0, coeff, B  
  common /constants/ Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma 

  REAL M_3, M_wd, R_wd, shock_height, vff, n_elec_shock, Tshock, ES0, MDOT0, XS0, coeff, B
  REAL Msun, Rsun, G, mue, mmw, A, mu, k, mH, alpha, beta, gamma

  REAL vff_old, h_old

  REAL RHO(VGRID),P(VGRID),TK(VGRID),X(VGRID),NELEC(VGRID),SOLN(VGRID)


  hstar = h
  h_old = h
  vff_old = vff
  DO i = 1, max_iter
    vff = ((2.*G*M_wd)*((1./(R_wd+hstar))))**(1./2.) ! Chuck's calculation for vff
    n_elec_shock = 4.0*7.013e23*MDOT0/vff
    Tshock = (3*G*M_wd*mu*mH)/(8*k*(R_wd+hstar))
    XS0 = vff**3*0.049/(2.*A*MDOT0) ! shock height [cm] when B = 0
    coeff = 9.1e-3*(B/10.)**2.85*(Tshock/1.e8)**2.0*(n_elec_shock/1.e16)**(-1.85)*(XS0/1.e7)**(-0.85)
    CALL Newton_solver(func, func_prime, 1.e5, 1.e-7, 100, ES0)

    CALL SHOCKSPEC_SHOOTING(VGRID,RHO,P,TK,X,NELEC,SOLN)
    hstar = X(VGRID)
    !print *,"shock_height ", X(VGRID)

    !IF ((R_wd+hstar).GT.(R_m)) THEN
    !  shock_height = h_old
    !  PRINT*,"ACCRETION COLUMN IS TOO TALL"
    !  RETURN
    !END IF

    IF (ABS((h_old-hstar)/h_old) <= tol) THEN
          RETURN
    END IF

      h_old = hstar
      vff_old = vff
  END DO
END SUBROUTINE SHOCK_SHOOT
