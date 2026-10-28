import numpy as np
import sympy as sym
import math
import constants as cn
import CoolProp.CoolProp as CP
from scipy.optimize import fsolve


# ----- Inputs ------
### Our inputs are limited to what we can set & record.
### Particularly, we can know:
### (1) the m which we set during fill (directly during testing & indirectly via calibration curve during launching),
### (2) the V which we set during fill (via F fill to hydrostatic jump setting piston position and therefore volumes), &
### (3) the P & T which we record during fill (via transducers on top bulkhead).

### After testing, these inputs can be refined, particularly those recordable (P & T).
### We may be able to generate models using this data, providing us P & T tables for given conditions.

## Initial Masses
### 2/3 <- (1) m_OX_o, (2) m_F_o, (3) m_OF_o
### where (1): initial OX mass, (2): initial F mass, (3): initial OX/F (O/F) mass ratio

## Initial Volumes
### 2/3 <- (1) V_OX_o, (2) V_F_o, (3) V_OF_o
### where (1): initial OX volume, (2): initial F volume, (3): initial OX/F (O/F) volume ratio

## Initial Thermodynamix
### 2/2 <- (1) P_OX_o, (2) T_OX_o
### where (1): initial OX pressure, (2): initial OX temp

# ----- Outputs -----
### Workin on it.

# ----- Standards & Constants -----
### We will stick to SI; any ambiguous units are assumed to be according to standard.
### Any non-standard units - US Customary or otherwise - must be converted to SI standard prior to use.
### SI Standard: s (seconds) [T], m (meter) [L], kg (kilogram) [M], A (ampere) [I], K (kelvin) [θ], mol (mole) [N], cd (candela) [J]

# ----- Fluids Function Definitions -----

# Specific Gravity (Density, Standard Density) [-]
def SG(rho: float) -> float:
    return rho / cn.CONST["rho_STD"]

# Dynamic Pressure (Density, Velocity) [Pa]
def DP(rho: float, v: float) -> float:
    return 0.5 * rho * v**2

# Hydrostatic Pressure (Density, Delta Z, Gravitational Constant) [Pa]
def HP(rho: float, dZ: float, g: float) -> float:
    return rho * dZ * g

# Velocity (Mass Flow, Density, Area) [m/s]
def v(mDot: float, rho: float, A: float) -> float:
    return mDot / (rho * A)

# Volumetric Flow (Mass Flow, Density) [m^3/s]
def Q(mDot: float, rho: float):
    return mDot / rho

# Reynold's Number (Density, Velocity, Inner Diameter, Dynamic Viscosity) [-]
def Re(rho: float, v: float, ID: float, mu: float) -> float:
    return rho * v * ID / mu

### viscosity calcs
def viscfuel(T,p):
    return CP.PropsSI('V','T', T, 'P', p, 'Ethanol')

def viscox(T,p):
    return CP.PropsSI('V','T',T,'P',p,'N2O')

#### velocity calcs
def velofuel(m_dot_fuel,T_fuel,P_fuel,pipe_area):
    dens = CP.PropsSI('D','T', T_fuel, 'P', P_fuel, 'Ethanol')
    vdot = m_dot_fuel/dens
    return vdot / pipe_area

def veloox(m_dot_ox,T_ox,P_ox,pipe_area):
    dens = CP.PropsSI('D','T', T_fuel, 'P', P_fuel, 'Ethanol')
    vdot = m_dot_fuel/dens
    return vdot / pipe_area

# Darcy Friction Factor (Reynold's Number, Roughness, Inner Diameter) [-]
def DFF(Re: float, rgh: float, ID: float) -> float:
    if Re < 2300:
        # Laminar Flow
        return 64 / Re
    else:
        # Turbulent Flow Approximation - Colebrook, Swamee-Jain Form
        return 0.25 / (math.log10(rgh / (3.7 * ID) + 5.74 / Re**0.9)**2)

# Mass Flow O/F Ratio (Mass Flow OX, Mass Flow F) [-]
def mDot_OF(mDot_OX: float, mDot_F: float) -> float:
    return mDot_OX / mDot_F

# Mass Flow Nozzle (Chamber Pressure, Throat Area, Characteristic Velocity) [kg/s]
def mDot_Nozzle(P_c: float, A_t: float, cstar: float) -> float:
    return P_c * A_t / cstar

# ----- Pressure Losses -----

def P_inj(P_tk: float, deltaP_line: float) -> float:
    P_inj = P_tk + deltaP_line
    return P_inj

def deltaP_inj(mdot: float, C_d: float, A: float, rho: float) -> float:
    deltaP_inj = - 0.5 * (mdot / (C_d * A))^2 / rho
    return deltaP_inj

# ----- Simulation -----

def cstar(P_c: float, OF: float) -> float:
    # Temporary; need to implement RPA interpolation.
    return 2639.4

## Residuals

def R_P(P, deltaP_line: float, P_c: float) -> float:
    return P - deltaP_line - P_c

def R_m(mdot_OX: float, mdot_F: float, P_c, A_t, cstar) -> float:
    R_m = mdot_OX + mdot_F - P_c*A_t/cstar
    return R_m

#### Line Pressure Drop 
### Here tons of constants needed -- listed below
#### NOTE: this needs to be a lot more fleshed out and specific -- rn I am trying to make a place holder to keep making engine solver

def line_loss(m_dot, rho, D, L, k, Cv, A_inj):
    dz = 0.9144 #### elevation change in m
    f = 0.04 ### darcy friction -- should be an input function but for simplicity leave it for rn
    # velocity and dynamic pressure
    A_pipe = (np.pi / 4) * D**2
    v = m_dot / (rho * A_pipe)
    q = 0.5 * rho * v**2
    G0 = 9.81

    dp_pipe = f * (L / D) * q   # pipe friction loss
    dp_k = k * q                # minor losses

    # valve loss: convert m_dot to GPM for the Cv equation
    SG = rho / 1000.0
    Q_gpm = (m_dot / rho) * 15850.32   # m^3/s -> GPM
    dp_v = SG * (Q_gpm / Cv) ** 2 * 6894.76  # psi -> Pa

    dp_inj = 0.5 * (m_dot / A_inj) ** 2 / rho  # injector loss
    dp_g = rho * G0 * dz                       # hydrostatic loss

    return dp_pipe + dp_k + dp_v + dp_inj + dp_g

def enginesolve(m_ox, m_fuel, p_fuel, p_ox, rho_ox, rho_f, T_ox, T_fuel, guess=None):
    # constants
    c_star = 2639.4          # m/s
    A_t = 7.85398e-5         # m^2
    A_e = 1.8816e-3          # m^2
    C_f = 1.5          # thrust coefficient, no idea if accurate
    g0 = 9.81                 # m/s^2

    #solve residuals
    ####NOTE: here I will input friction factors / other constants into the line_loss formula to get place holder values -- NEEDS more detail
    def deltap_line_ox(m_dot_ox, rho_ox): 
        return line_loss(m_dot_ox, rho_ox, 0.0535088, 0.762, 10, 9, 2.20893e-5)
    def deltap_line_fuel(m_dot_fuel, rho_f):
        return line_loss(m_dot_fuel, rho_f, 0.0535088, 0.762, 10, 5.228, 1.8816e-3)
    #### NOTE: DELTAPLINES ARE NOT FLESHED OUT -- COME BACK TO THIS
    ### fsolve mech here because sympy solver took too long
    def residuals(x):
        m_dot_ox, m_dot_fuel, P_c = x
        return [
            R_P(p_ox,   deltap_line_ox(m_dot_ox, rho_ox), P_c),
            R_P(p_fuel, deltap_line_fuel(m_dot_fuel, rho_f), P_c),
            R_m(m_dot_ox, m_dot_fuel, P_c, A_t, c_star),
        ]
    if guess is None:
        guess = [0.5, 0.2, 0.5 * p_ox]
    sol, info, ier, msg = fsolve(residuals, guess, full_output=True)
    if ier !=1:
        raise RuntimeError
    
    m_dot_ox = sol[0]
    m_dot_fuel = sol[1]
    Pc = sol[2]
    thrust = C_f * Pc * A_t
    ox_fuel_ratio = m_dot_ox / m_dot_fuel
    spec_impulse = thrust / ((m_dot_fuel+m_dot_ox)*g0)

    return (m_dot_ox, m_dot_fuel, Pc, thrust, ox_fuel_ratio, spec_impulse)