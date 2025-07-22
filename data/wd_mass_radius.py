import argparse
import numpy as np
from scipy.integrate import solve_ivp

parser = argparse.ArgumentParser(prog="White Dwarf mass-radius Solver",
                    description="Produces a list of masses and radii for a range of central densities")
parser.add_argument("-Z", type=float, required=True, help="charge of nucleus of WD material")
parser.add_argument("-A", type=float, help="mass of nucleus of WD material (in amu), superceded by mu if input")
parser.add_argument("-mu", type=float, help="A/Z ratio for WD material, if not specified A/Z is computed from other inputs")
args = parser.parse_args()

if args.mu is None:
    args.mu = args.A/args.Z

class ode_event():
    def __init__(self, func):
        self.func = func
    def __call__(self, t, y, args):
        return self.func(t,y,args)
    terminal = False

# physical constants in cgs units
c = 29979245800.0
hbar = 1.0545718176461565e-27
m_e = 9.1093837015e-28
amu = 1.66053906892e-24
alpha = 0.0072973525693
G = 6.674299999999999e-08
M_sun = 1.988409870698051e+33
R_sun = 69570000000.0

def dpressure(r, state, Z):
    x,m = state
    y = np.sqrt(1+x**2)
    dy = x/y
    b = x+y
    db = 1 + dy

    df = .125*(b**3 - b**-5) - .25*(b - b**-3) - 1.5*np.log(b)*(b - 2/b + b**-3)
    df2 = 3*b**-2 + .125*(3*b**2 + 5*b**-6) - .25*(7 + 9*b**-4) - 1.5*np.log(b)*(1 + 2*b**-2 - 3*b**-4)
    Q = alpha*np.cbrt(4*(Z**2)/(9*np.pi))

    dP_0 = dy*x**3
    dP_tf = -(54/35)*Q*Q*(x**3)*((7/(9*Q)) + dy - (dy**3)/5)
    dP_ex = -(alpha/(4*np.pi))*db*((dy**2 - dy + 2)*df - b*dy*df2)
    dP_cor = -0.0311*alpha*alpha*(x**2)

    dP = dP_0 + dP_tf + dP_ex + dP_cor

    return dP

def pressure(r, state, Z):
    x,m = state
    y = np.sqrt(1+x**2)
    dy = x/y
    b = x+y

    Q = alpha*np.cbrt(4*(Z**2)/(9*np.pi))

    f = (b**4 + b**-4)/32 + (b**2 + b**-2)/4 - 9/16 - .75*np.log(b)*(b**2 - b**-2) + 1.5*np.log(b)**2
    df = .125*(b**3 - b**-5) - .25*(b - b**-3) - 1.5*np.log(b)*(b - 2/b + b**-3)

    P_0 = .125*(y*(2*x**3 - 3*x) + 3*np.arcsinh(x))
    P_tf = -1.2*Q*((x**4)/4 + (9*Q/35)*dy*(x**4))
    P_ex = -(alpha/(4*np.pi))*(3*f - b*dy*df)
    P_cor = -(0.0311/3)*alpha*alpha*(x**3)

    return P_0 + P_tf + P_ex + P_cor

def momentum(r, state, Z):
    x,m = state
    return x

zero_pressure = ode_event(pressure)
zero_pressure.terminal = True
zero_momentum = ode_event(momentum)
zero_momentum.terminal = True

def equation_of_state(r, state, Z):
    x,m = state

    dP = dpressure(r, state, Z)

    dx_dr = -m*(x**3)/(dP*(r**2))
    dm_dr = 4*np.pi*(r**2)*(x**3)

    return dx_dr, dm_dr


pressre_coeff = ((m_e*c**2)/(3*np.pi**2))*(m_e*c/hbar)**3
density_coeff = (args.mu*amu/(3*np.pi**2))*(m_e*c/hbar)**3

n_points = 200
initial_x = np.logspace(-1,3,n_points//2)
mass = np.zeros(n_points)
radius = np.zeros(n_points)

r0 = 1e-32
r_max = density_coeff*np.sqrt(G/pressre_coeff)*R_sun
for i, x0 in enumerate(initial_x):
    m0 = (4/3)*np.pi*(x0*r0)**3
    res = solve_ivp(equation_of_state, (r0, r_max), [x0,m0], args=(args.Z,), rtol=1e-10, atol=1e-16, events=[zero_pressure, zero_momentum])

    x_f = res.y[0][-1]
    m_f = res.y[1][-1]
    r_f = res.t[-1]
    p_f = pressure(r_f, (x_f,m_f), args.Z)
    dp_dx = dpressure(r_f, (x_f,m_f), args.Z)
    dx_dr, dm_dr = equation_of_state(r_f, (x_f,m_f), args.Z)

    dx = min(-p_f/dp_dx, -x_f, key=abs)
    dr = dx/dx_dr
    dm = dm_dr*dr

    mass[i] = (m_f+dm)*np.sqrt((pressre_coeff/G)**3)/density_coeff**2
    radius[i] = (r_f+dr)*np.sqrt(pressre_coeff/G)/density_coeff

mass = mass/M_sun
radius = radius/R_sun

initial_x = np.concatenate((initial_x, np.zeros(n_points//2)))
for i in range(n_points//2, n_points):
    x0 = 10**((np.log10(initial_x[np.argmax(np.diff(mass[:i]))]) + np.log10(initial_x[np.argmax(np.diff(mass[:i]))+1]))/2)
    initial_x[i] = x0
    m0 = (4/3)*np.pi*(x0*r0)**3
    res = solve_ivp(equation_of_state, (r0, r_max), [x0,m0], args=(args.Z,), rtol=1e-10, atol=1e-16, events=[zero_pressure, zero_momentum])

    x_f = res.y[0][-1]
    m_f = res.y[1][-1]
    r_f = res.t[-1]
    p_f = pressure(r_f, (x_f,m_f), args.Z)
    dp_dx = dpressure(r_f, (x_f,m_f), args.Z)
    dx_dr, dm_dr = equation_of_state(r_f, (x_f,m_f), args.Z)

    dx = min(-p_f/dp_dx, -x_f, key=abs)
    dr = dx/dx_dr
    dm = dm_dr*dr

    mass[i] = ((m_f+dm)*np.sqrt((pressre_coeff/G)**3)/density_coeff**2)/M_sun
    radius[i] = ((r_f+dr)*np.sqrt(pressre_coeff/G)/density_coeff)/R_sun

    sorted_inds = np.argsort(mass[:i+1])
    initial_x[:i+1] = initial_x[:i+1][sorted_inds]
    radius[:i+1] = radius[:i+1][sorted_inds]
    mass[:i+1] = np.sort(mass[:i+1])

with open("mass_radius.txt", "w+") as out_file:
    out_file.write("Mass [Solar Mass], Radius [Solar Radii]\n")
    for m,r in zip(mass, radius):
        out_file.write(f"{m}, {r}\n")
