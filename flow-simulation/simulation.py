import numpy as np
from scipy.optimize import root

import functions as fn
import constants as cn
from classes import *

# Simulation Settings
dt = 0.01
t_max = 30
RAE = 0.001
i_max = 100

# Initialize Fluids & States
OX = Fluid(name = "OX", REF = "NitrousOxide")
F = Fluid(name = "F", REF = "Ethanol")

fluid_state_OX_o = OX.state_PT(P = P_OX_o, T = T_OX_o)
fluid_state_F_o = F.state_PT(P = P_F_o, T = T_F_o)

# Initialize Tanks & States
tank_state_OX_o = TankState(
    m = m_OX_o,
    V = V_OX_o,
    U = m_OX_o * fluid_state_OX_o.u,
    fluid_state = fluid_state_OX_o
)

tank_state_F_o = TankState(
    m = m_F_o,
    V = V_F_o,
    U = m_F_o * fluid_state_F_o.u,
    fluid_state = fluid_state_F_o
)

Tank_OX = Tank(
    name = "OX Tank",
    fluid = OX,
    tank_state = tank_state_OX_o
)

Tank_F = Tank(
    name = "F Tank",
    fluid = F,
    tank_state = tank_state_F_o
)

# Initialize Lines
Line_OXL = Line(
    name = "OX Line",
    source = Tank_OX
)

Line_F = Line(
    name = "F Line",
    source = Tank_F
)

# Initialize Line Components

# -----

# Injectors

Injector_OX = Injector(
    C_d = C_d_OX,
    A = A_OX
)

Injector_F = Injector(
    C_d = C_d_F,
    A = A_F
)