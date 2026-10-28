import numpy as np
import sympy as sym
from scipy.optimize import root
import matplotlib.pyplot as plt

import functions as fn
import constants as cn
from classes import *
from CoolProp.CoolProp import PropsSI

# Simulation Settings
dt = 0.01
t_max = 300
RAE = 0.001
i_max = 100

#### Initial Conditions
t_amb = 300
v_ox_initial = 0.1 #### m^3 so 100 L
v_fuel_initial = 0.06 ##### m^3 so 60 L
piston_loss_energy = 1 #### J
piston_loss_pressure = 1 ### Pa
p_ox_initial = PropsSI('P','T',t_amb,'Q',0,'N2O')  
p_fuel_initial = p_ox_initial - piston_loss_pressure
m_ox_initial = v_ox_initial * PropsSI('D','T',t_amb,'Q',0,'N2O') 
m_fuel_initial = v_fuel_initial * PropsSI('D','T',t_amb,'P',p_fuel_initial,'Ethanol')
dens_fuel = PropsSI('D','T',t_amb,'P',p_fuel_initial,'Ethanol') ### assumed constant bc liquid density change minimal
u_ox_initial = m_ox_initial*PropsSI('U','T',t_amb,'Q',0,'N2O')  #### total U NOT specific
u_fuel_initial = m_fuel_initial*PropsSI('U','T',t_amb,'P',p_fuel_initial,'Ethanol')

# Initialize Fluids & States
OX = Fluid(name = "OX", heos = "NitrousOxide")
F = Fluid(name = "F", heos = "Ethanol")

fluid_state_OX_o = OX.state_QT(T = t_amb, Q = 0)
fluid_state_F_o = F.state_PT(P = p_fuel_initial, T = t_amb)

# Initialize Tanks & States
tank_state_OX_o = TankState(
    m = m_ox_initial,
    V = v_ox_initial,
    U = m_ox_initial * fluid_state_OX_o.u,
    fluid_state = fluid_state_OX_o
)

tank_state_F_o = TankState(
    m = m_fuel_initial,
    V = v_fuel_initial,
    U = m_fuel_initial * fluid_state_F_o.u,
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
## constant values here wrong but just to let it work
C_d_OX = 999999999
C_d_F = 999999999
A_OX = 999999999
A_F = 999999999
Injector_OX = Injector(
    C_d = C_d_OX,
    A = A_OX
)

Injector_F = Injector(
    C_d = C_d_F,
    A = A_F
)

##### property arrays
time_array = np.arange(0,t_max + dt, dt)
m_dot_ox_array = np.zeros(len(time_array))
m_dot_fuel_array = np.zeros(len(time_array))
chamber_pressure_array = np.zeros(len(time_array))
thrust_array = np.zeros(len(time_array))
ox_fuel_ratio_array = np.zeros(len(time_array))
specific_impulse_array = np.zeros(len(time_array))
m_ox_array = np.zeros(len(time_array))
m_fuel_array = np.zeros(len(time_array))
p_ox_array = np.zeros(len(time_array))
p_fuel_array = np.zeros(len(time_array))
u_ox_array = np.zeros(len(time_array))
u_fuel_array = np.zeros(len(time_array))
T_fuel_array = np.zeros(len(time_array))
T_ox_array = np.zeros(len(time_array))
v_ox_array = np.zeros(len(time_array))


### adding initials to arrays 
T_ox_array[0] = t_amb
T_fuel_array[0] = t_amb
p_ox_array[0] = p_ox_initial
p_fuel_array[0] = p_fuel_initial
m_ox_array[0] = m_ox_initial
m_fuel_array[0] = m_fuel_initial
u_ox_array[0] = u_ox_initial
u_fuel_array[0] = u_fuel_initial
v_ox_array[0] = v_ox_initial

#### fuel density assumed constant
fueldens = m_fuel_initial / v_fuel_initial

###### sim here
for i in range(1,(len(time_array))):
    engine_solution = fn.enginesolve(m_ox_array[i-1], m_fuel_array[i-1], p_fuel_array[i-1], p_ox_array[i-1], PropsSI('D','T', T_ox_array[i-1], 'Q', 0, 'N2O'), fueldens, T_ox_array[i-1], T_fuel_array[i-1], guess=None)
    m_dot_ox_array[i] = engine_solution[0]
    m_dot_fuel_array[i] = engine_solution[1]
    chamber_pressure_array[i] = engine_solution[2]
    thrust_array[i] = engine_solution[3]
    ox_fuel_ratio_array[i] = engine_solution[4]
    specific_impulse_array[i] = engine_solution[5]
    m_ox_array[i] = m_ox_array[i-1] - m_dot_ox_array[i-1]*dt
    m_fuel_array[i] = m_fuel_array[i-1] - m_dot_fuel_array[i-1]*dt
    deltav_fuel = (-m_dot_fuel_array[i-1]*dt) / fueldens
    v_ox_array[i] = v_ox_array[i-1] - deltav_fuel
    work_into_fuel = p_fuel_array[i-1] * -deltav_fuel   ##### yields positive value
    work_out_of_ox = work_into_fuel + piston_loss_energy  ### also positive
    delta_u_ox = -work_out_of_ox - (dt*m_dot_ox_array[i-1]*(u_ox_array[i-1] / m_ox_array[i-1]))
    u_ox_array[i] = u_ox_array[i-1] + delta_u_ox
    delta_u_f = work_into_fuel - (dt*m_dot_fuel_array[i-1]*(u_fuel_array[i-1] / m_fuel_array[i-1]))
    u_fuel_array[i] = u_fuel_array[i-1] + delta_u_f
    p_ox_array[i] = PropsSI('P','U',(u_ox_array[i]/m_ox_array[i]),'D',(m_ox_array[i]/v_ox_array[i]),'N2O')
    p_fuel_array[i] = p_ox_array[i] - piston_loss_pressure
    T_ox_array[i] = PropsSI('T','U',(u_ox_array[i]/m_ox_array[i]),'P',p_ox_array[i],'N2O')
    T_fuel_array[i] = PropsSI('T','U',(u_fuel_array[i]/m_fuel_array[i]),'P',p_fuel_array[i],'Ethanol')

### plots
fig,axes = plt.subplots(nrows=2, ncols=3, figsize=(10,8))
axes[0,0].plot(time_array,chamber_pressure_array)
axes[0,0].set_title('Chamber Pressure')

axes[0,1].plot(time_array,thrust_array)
axes[0,1].set_title('Thrust')

axes[0,2].plot(time_array,ox_fuel_ratio_array)
axes[0,2].set_title('Oxidizer-Fuel Ratio')

axes[1,0].plot(time_array,specific_impulse_array)
axes[1,0].set_title('Specific Impulse')

axes[1,1].plot(time_array,p_ox_array)
axes[1,1].set_title('Oxidizer Pressure')

axes[1,2].plot(time_array,p_fuel_array)
axes[1,2].set_title('Fuel Pressure')

plt.tight_layout()

plt.show()