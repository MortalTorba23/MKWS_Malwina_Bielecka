import cantera as ct
import numpy as np
import matplotlib.pyplot as plt

def combust(T: float, P: float, spcies_dict: dict, reactor: str, dt: float= 1e-6, till_ignition: bool=True):
    gas = ct.Solution('gri30.yaml')
    gas.TPX = T, P, spcies_dict
    
    if reactor == 'constant temperature':
        r = ct.IdealGasReactor(gas)
    elif reactor == 'constant pressure':
        r = ct.IdealGasConstPressureReactor(gas)
    else:
        raise TypeError(f"Improper reactor type '{reactor}'")
    
    sim = ct.ReactorNet([r])
    states = ct.SolutionArray(gas, extra=['time_in_ms', 'time_in_sec'])

    simulation_time = 0.01  
    n_steps = int(simulation_time/dt)
    ignition_delay = 0.0
    ignited = False
    
    time = 0.0
    for n in range(n_steps):
        time += dt
        sim.advance(time)
        states.append(r.thermo.state, time_in_sec= time, time_in_ms= time*1e3)
        
        if not ignited:
            if r.T >= (T + 400):
                ignition_delay = time
                ignited = True
                if till_ignition:
                    break

    return gas, states, ignition_delay

# PART I-const temperature
fig1, axs = plt.subplots(2,1, figsize=(10, 10), layout='tight')

P = np.linspace(1,5,9)
p_plot = []
igd_plot_p = []

for t in P:
    _, states, igd = combust(1250, t*ct.one_atm, {'H2':2, 'O2':1, 'N2': 3.76}, 'constant temperature', dt= 1e-6)
    
   
    if igd > 0:
        lbl = f'{t: .3} atm, I.D. = {igd*1000: .5} ms'
        p_plot.append(t)
        igd_plot_p.append(igd * 1000)
        axs[0].plot(states.time_in_ms, states.T, '.-', label=lbl)

axs[0].set_xlabel(r'Time [ms]$\longrightarrow$')
axs[0].set_ylabel(r"Flame Temperature [K]$\longrightarrow$")
axs[0].grid(linestyle='-.')
axs[0].legend()

axs[1].plot(igd_plot_p, p_plot, '.-')
axs[1].set_xlabel(r'Auto-ignition time [ms]$\longrightarrow$')
axs[1].set_ylabel(r"Pressure [atm]$\longrightarrow$")
axs[1].grid(linestyle='-.')

fig1.suptitle(f"Auto-ignition of hydrogen at 1250 K and pressure variation")


fig2, axs = plt.subplots(2,1, figsize=(10, 10), layout='tight')

T = np.linspace(950, 1450, 9)
t_plot = []
igd_plot_t = []

for t in T:
    _, states, igd = combust(t, 5*ct.one_atm, {'H2':2, 'O2':1, 'N2': 3.76}, 'constant pressure', dt= 1e-6)
    
    
    if igd > 0:
        lbl = f'{t} K, I.D. = {igd*1000: .5} ms'
        t_plot.append(t)
        igd_plot_t.append(igd * 1000)
        axs[0].plot(states.time_in_ms, states.T, '.-', label=lbl)

axs[0].set_xlabel(r'Time [ms]$\longrightarrow$')
axs[0].set_ylabel(r"Flame Temperature [K]$\longrightarrow$")
axs[0].grid(linestyle='-.')
axs[0].legend()

axs[1].plot(igd_plot_t, t_plot, '.-')
axs[1].set_xlabel(r'Auto-ignition time [ms]$\longrightarrow$')
axs[1].set_ylabel(r"Initial temperature [K]$\longrightarrow$")
axs[1].grid(linestyle='-.')

fig2.suptitle(f"Auto-ignition of hydrogen at 5 atm and initial temperature variation")

plt.show()


fig4, axs = plt.subplots(1,1, figsize=(10, 10), layout='tight')

sim_iters = 5000  
times = np.zeros(sim_iters)
gas = ct.Solution('gri30.yaml')
data = np.zeros([sim_iters,4])

def check_time(tab, T_start):
    for i in range(len(tab)):
        if tab[i, 0] >= T_start + 400:
            return i
    return -1

fuel_mole = np.array([0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
P = np.linspace(1,5,9)
temperature = 1250

print('%15s    %15s' % ('H2 [mole]', 'Ignition time [s]'))

for t in P:
    fuel_plot = []
    time_plot = []
    for m in fuel_mole:
        gas.TPX = temperature, t*ct.one_atm, 'H2:%f, O2:1, N2:3.76' % m
        r = ct.IdealGasReactor(gas)
        sim = ct.ReactorNet([r])
        time = 0
        for n in range(sim_iters):
            time += 2.e-6  
            sim.advance(time)
            data[n,0] = r.T
            times[n] = time
            
        idx = check_time(data, temperature)
        if idx != -1:
            fuel_plot.append(m)
            time_plot.append(times[idx])
            if t == P[0]:
                print('%15.2f %15.6f' % (m, times[idx]))
        else:
            if t == P[0]:
                print('%15.2f %15s' % (m, 'No ignition'))
                
    axs.plot(fuel_plot, time_plot, '.-', label= f'{t: .3} atm')

axs.set_xlabel(r'Fuel quantity [mole]$\longrightarrow$')
axs.set_ylabel(r"Ignition time [s]$\longrightarrow$")
axs.grid(linestyle='-.')
axs.legend()
fig4.suptitle(f"Ignition time at %dK with varying fuel quantity and pressure" % temperature)
    
plt.show()