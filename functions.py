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
    deltaP_inj = - 0.5 * (mdot / (C_d * A)) ** 2 / rho
    return deltaP_inj

# ----- Simulation -----

# ----- c* data [m/s] from NASA CEA (N2O / 100% ethanol), generated via rocketcea -----
# always list in ascending order pls
_CS_OF = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0]
_CS_PC = [2000000.0, 3000000.0, 4000000.0, 5000000.0]
_CS_DATA = [[1355.1, 1355.1, 1355.1, 1355.1], [1452.1, 1452.1, 1452.1, 1452.0], [1513.8, 1513.9, 1513.9, 1514.0], [1553.8, 1554.5, 1554.9, 1555.2], [1576.3, 1578.6, 1580.0, 1581.0], [1582.1, 1586.3, 1589.1, 1591.1], [1576.3, 1581.5, 1585.1, 1587.8], [1565.1, 1570.6, 1574.4, 1577.3], [1552.1, 1557.5, 1561.2, 1564.0], [1538.7, 1543.8, 1547.2, 1549.9], [1525.5, 1530.1, 1533.3, 1535.8], [1512.6, 1516.9, 1519.8, 1522.0], [1500.3, 1504.2, 1506.8, 1508.7]]

# Characteristic Velocity (Chamber Pressure [Pa], O/F [-]) [m/s]
def cstar(P_c: float, OF: float) -> float:
    """Interpolated c* [m/s]. Inputs outside the table are clamped to its edges."""
    OF = min(max(OF, _CS_OF[0]), _CS_OF[-1])
    P_c = min(max(P_c, _CS_PC[0]), _CS_PC[-1])
    col = [np.interp(OF, _CS_OF, [row[j] for row in _CS_DATA]) for j in range(len(_CS_PC))]
    return float(np.interp(P_c, _CS_PC, col))

## Residuals

def R_P(P_inj: float, deltaP_inj: float, P_c: float) -> float:
    return P_inj + deltaP_inj - P_c

def R_m(mdot_OX: float, mdot_F: float, mdot_nozzle: float) -> float:
    R_m = mdot_OX + mdot_F - mdot_nozzle
    return R_m

# ----- Guesses -----

def guesses(P_tk: float, deltaP_line: float, C_d: float, A: float, rho: float, P_c: float, OF: float) -> tuple:
    P_inj_guess = P_inj(P_tk, deltaP_line)
    deltaP_inj_guess = deltaP_inj(1.0, C_d, A, rho)
    mdot_nozzle_guess = mDot_Nozzle(P_c, A, cstar(P_c, OF))
    return P_inj_guess, deltaP_inj_guess, mdot_nozzle_guess

def residuals(P_inj: float, deltaP_inj: float, P_c: float, mdot_OX: float, mdot_F: float, mdot_nozzle: float) -> tuple: 
    R_P_guess= R_P(P_inj, deltaP_inj, P_c)
    R_m_guess = R_m(mdot_OX, mdot_F, mdot_nozzle)
    return R_P_guess, R_m_guess

def solve(P_tk: float, deltaP_line: float, C_d: float, A: float, rho: float, P_c: float, OF: float, mdot_OX: float, mdot_F: float) -> tuple:
    P_inj_guess, deltaP_inj_guess, mdot_nozzle_guess = guesses(P_tk, deltaP_line, C_d, A, rho, P_c, OF)
    R_P_guess, R_m_guess = residuals(P_inj_guess, deltaP_inj_guess, P_c, mdot_OX, mdot_F, mdot_nozzle_guess)
    return P_inj_guess, deltaP_inj_guess, mdot_nozzle_guess, R_P_guess, R_m_guess

def main(P_tk: float, deltaP_line: float, C_d: float, A: float, rho: float, P_c: float, OF: float, mdot_OX: float, mdot_F: float) -> tuple:
    P_inj_guess, deltaP_inj_guess, mdot_nozzle_guess, R_P_guess, R_m_guess = solve(P_tk, deltaP_line, C_d, A, rho, P_c, OF, mdot_OX, mdot_F)
    return P_inj_guess, deltaP_inj_guess, mdot_nozzle_guess, R_P_guess, R_m_guess

def run_simulation(P_tk: float, deltaP_line: float, C_d: float, A: float, rho: float, P_c: float, OF: float, mdot_OX: float, mdot_F: float) -> tuple:
    P_inj_guess, deltaP_inj_guess, mdot_nozzle_guess, R_P_guess, R_m_guess = main(P_tk, deltaP_line, C_d, A, rho, P_c, OF, mdot_OX, mdot_F)
    return P_inj_guess, deltaP_inj_guess, mdot_nozzle_guess, R_P_guess, R_m_guess   

def intialstates(P_tk: float, deltaP_line: float, C_d: float, A: float, rho: float, P_c: float, OF: float, mdot_OX: float, mdot_F: float) -> tuple:
    P_inj_guess, deltaP_inj_guess, mdot_nozzle_guess, R_P_guess, R_m_guess = run_simulation(P_tk, deltaP_line, C_d, A, rho, P_c, OF, mdot_OX, mdot_F)
    return P_inj_guess, deltaP_inj_guess, mdot_nozzle_guess, R_P_guess, R_m_guess

def run_intialstates(P_tk: float, deltaP_line: float, C_d: float, A: float, rho: float, P_c: float, OF: float, mdot_OX: float, mdot_F: float) -> tuple:
    P_inj_guess, deltaP_inj_guess, mdot_nozzle_guess, R_P_guess, R_m_guess = intialstates(P_tk, deltaP_line, C_d, A, rho, P_c, OF, mdot_OX, mdot_F)
    return P_inj_guess, deltaP_inj_guess, mdot_nozzle_guess, R_P_guess, R_m_guess