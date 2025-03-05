#ifndef TABLEAU_H
#define TABLEAU_H
const int order = 5;
const int n_stages = 6;
inline double a[6][order] = {{0,0,0,0,0},
                  {1./5.,0,0,0,0},
                  {3./40., 9./40., 0,0,0},
                  {44./45., -56./15., 32./9., 0,0},
                  {19372./6561., -25360./2187., 64448./6561., -212./729.,0},
                  {9017./3168., -355./33., 46732./5247., 49./176., -5103./18656.}};

inline double b[6] = {35./384., 0., 500./1113., 125./192., -2187./6784., 11./84.};

inline double c[7] = {0.,1./5.,3./10.,4./5.,8./9.,1.,1.};

inline double e[7] = {-71./57600., 0., 71./16695.,-71./1920., 17253./339200.,-22./525., 1./40.};

inline double p[7][4] = {{1.,-8048581381./2820520608., 8663915743./2820520608.,-12715105075./11282082432.},
                  {0., 0., 0., 0.},
                  {0., 131558114200./32700410799., -68118460800./10900136933., 87487479700./32700410799.},
                  {0.,-1754552775./470086768., 14199869525./1410260304.,-10690763975./1880347072.},
                  {0., 127303824393./49829197408.,-318862633887./49829197408., 701980252875./199316789632.},
                  {0.,-282668133./205662961., 2019193451./616988883.,-1453857185./822651844.},
                  {0., 40617522./29380423.,-110615467./29380423., 69997945./29380423.}};
#endif
