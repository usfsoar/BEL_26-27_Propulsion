import numpy as np
import math
import constants as cn

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

def R_P(P_inj: float, deltaP_inj: float, P_c: float) -> float:
    return P_inj + deltaP_inj - P_c

def R_m(mdot_OX: float, mdot_F: float, mdot_nozzle: float) -> float:
    R_m = mdot_OX + mdot_F - mdot_nozzle
    return R_m