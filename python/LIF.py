import numpy as np
import matplotlib.pyplot as plt

## parameters
T = 50  # msec
dt = 0.125  # simulation time step
time = np.arange(0, T + dt, dt)
t_rest = 0  # initial refractory time

## LIF properties
Vm = np.zeros(len(time))  # potential (V) trace over time
Rm = 1  # resistance (kOhm)
Cm = 10  # capacitance (uF)
tau_m = Rm * Cm  # time constant (msec)
tau_ref = 4  # refractory period (msec)
Vth = 0.5  # spike threshold (V)
V_spike = 1.0  # spike delta (V)

## Input stimulus
I = 1.5  # input current (A)

## Iterate over each time step
for i, t in enumerate(time):
    if t > t_rest:
        Vm[i] = Vm[i - 1] + (-Vm[i - 1] + I * Rm) / tau_m * dt
        if Vm[i] >= Vth:
            Vm[i] += V_spike
            t_rest = t + tau_ref

## Plot membrane potential trace
plt.plot(time, Vm)
plt.title('Leaky Integrate-and-Fire Example')
plt.ylabel('Membrane Potential (V)')
plt.xlabel('Time (msec)')
plt.ylim([0, 2])
plt.show()
