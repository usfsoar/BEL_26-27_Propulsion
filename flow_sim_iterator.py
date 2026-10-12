import numpy as np
import matplotlib.pyplot as plt
from enginesolver import engine_solver
from CoolProp.CoolProp import PropsSI

##### sim steps
dt = 0.01
t_max = 15
num_steps = int(t_max / dt)

#####initial values here --- NOTE: we assume initial full tanks of ox and fuel, using vox and vfuel
v_ox_initial= 0.1 ##### m^3 which is L/1000
v_fuel_initial= 0.06   ##### also m^3
piston_loss_energy= 1 ###### J
piston_loss_pressure= 1  ###### Pa ####NOTE:these constants are related but for now I leave them both as inputs
T_amb=300

time_array = np.linspace(0,t_max,num_steps)
values_matrix = np.zeros((15, num_steps))  ##### ORDERING: 0) m_dot_ox 1) m_dot_fuel 2) p_chamber 3) thrust 4) ox_fuel_ratio
#### 5) specific_impulse 6) m_ox 7) m_fuel 8) p_ox 9) p_fuel 10) u_ox 11) u_fuel 12) t_ox 13) t_fuel 14) v_ox_tank

##### Calculated / Other Initial Conditions  NOTE: calculated m_dots would be more accurate but this is good with small time steps
p_ox_initial = PropsSI('P','T',T_amb,'Q',0,'N2O')  ####here we need refprop call
p_fuel_initial = p_ox_initial - piston_loss_pressure
m_ox_initial = v_ox_initial * PropsSI('D','T',T_amb,'Q',0,'N2O') ###assume sat liquid ox initial state  units =kg
m_fuel_initial = v_fuel_initial * PropsSI('D','T',T_amb,'P',p_fuel_initial,'Ethanol')
values_matrix[8,0] = p_ox_initial
values_matrix[9,0] = p_fuel_initial
values_matrix[6,0] = m_ox_initial
values_matrix[7,0] = m_fuel_initial
dens_fuel = PropsSI('D','T',T_amb,'P',p_fuel_initial,'Ethanol') ### assumed constant bc liquid density change minimal
u_ox_initial = m_ox_initial*PropsSI('U','T',T_amb,'Q',0,'N2O')
u_fuel_initial = m_fuel_initial*PropsSI('U','T',T_amb,'P',p_fuel_initial,'Ethanol')
values_matrix[10,0] = u_ox_initial
values_matrix[11,0] = u_fuel_initial
values_matrix[12,0] = T_amb
values_matrix[13,0] = T_amb
values_matrix[14,0] = v_ox_initial
#######iterations here now
for i in range(1, num_steps):
    prev_values = values_matrix[:,i-1]
    solv_array = engine_solver(prev_values[6],prev_values[7],prev_values[8],prev_values[9],(PropsSI('D','U',(prev_values[10]/prev_values[6]),'P',prev_values[8],'N2O')),dens_fuel)
    ######solv_array returns (m_dot_ox, m_dot_f, chamber pressure, thrust, ox-fuel ratio, specific impulse)
    values_matrix[0,i] = solv_array[0]  #### m_dot_ox
    values_matrix[1,i] = solv_array[1]  #### m_dot_fuel
    values_matrix[2,i] = solv_array[2]  #### chamber_pres
    values_matrix[3,i] = solv_array[3]  #### thrust
    values_matrix[4,i] = solv_array[4]  #### ox-fuel ratio
    values_matrix[5,i] = solv_array[5]  #### specific impulse
    values_matrix[6,i] = -prev_values[0]*dt + prev_values[6]  #### new m_ox = old m_ox - old m_dot_ox * time step
    values_matrix[7,i] = -prev_values[1]*dt + prev_values[7]  #### new m_f = old m_f - old m_dot_f * time step
    delta_v_fuel = (-prev_values[1]*dt)/dens_fuel
    values_matrix[14,i] = prev_values[14] - delta_v_fuel
    work_into_fuel = prev_values[8]*-delta_v_fuel    #### work received by fuel = p_fuel*delta_v   NOTE: positive
    work_out_of_ox = work_into_fuel + piston_loss_energy
    delta_u_ox = -work_out_of_ox-(prev_values[0]*dt)*(prev_values[10]/prev_values[6]) #### dU (ox) = Wpiston out + m_dot_ox*(old U_ox/old m_ox)
    values_matrix[10,i] = prev_values[10] + delta_u_ox
    delta_u_fuel = (-prev_values[1]*dt)*(prev_values[11]/prev_values[7]) + work_into_fuel
    values_matrix[11,i] = prev_values[11] + delta_u_fuel
    values_matrix[8,i] = PropsSI('P','U',(values_matrix[10,i]/values_matrix[6,i]),'D',(values_matrix[6,i]/values_matrix[14,i]),'N2O') ##### Q is not accurate but we just need Psat
    values_matrix[9,i] = values_matrix[8,i] - piston_loss_pressure
    values_matrix[12,i] = PropsSI('T','U',(values_matrix[10,i]/values_matrix[6,i]),'P',values_matrix[8,i],'N2O')
    values_matrix[13,i] = PropsSI('T','U',(values_matrix[11,i]/values_matrix[7,i]),'P',values_matrix[9,i],'Ethanol')

######plotting info here
chamber_pressure_plot = values_matrix[3,:]
thrust_plot = values_matrix[4,:]
oxidizer_fuel_ratio_plot = values_matrix[5,:]
specific_impulse_plot = values_matrix[6,:]
oxidizer_pressure_plot = values_matrix[9,:]
fuel_pressure_plot = values_matrix[10,:]

fig,axes = plt.subplots(nrows=2, ncols=3, figsize=(10,8))
axes[0,0].plot(time_array,chamber_pressure_plot)
axes[0,0].set_title('Chamber Pressure')

axes[0,1].plot(time_array,thrust_plot)
axes[0,1].set_title('Thrust')

axes[0,2].plot(time_array,oxidizer_fuel_ratio_plot)
axes[0,2].set_title('Oxidizer-Fuel Ratio')

axes[1,0].plot(time_array,specific_impulse_plot)
axes[1,0].set_title('Specific Impulse')

axes[1,1].plot(time_array,oxidizer_pressure_plot)
axes[1,1].set_title('Oxidizer Pressure')

axes[1,2].plot(time_array,fuel_pressure_plot)
axes[1,2].set_title('Fuel Pressure')

plt.tight_layout()

plt.show()
