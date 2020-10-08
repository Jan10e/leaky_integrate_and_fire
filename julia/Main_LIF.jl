# Leaky integrate and fire neuron
# (this is same code as in Python to compare)
#
# Name: Jantine Broek
# Date: 13 June 2020
##################################################################
## working directory
cd("/Users/jantinebroek/Documents/03_projects/04_IF/")


## Loads packages
using Plots         # Select the Plots package for figures
pyplot()            # Plots package will use Pyplot (matplotlib needs to be installed).


## parameters
T = 50              # msec
dt = 0.125          # simulation time step
t_rest = 0          # initial refractory time
time = collect(0:dt:T+dt)

## LIF properties
Rm = 1          # resistance (kOhm)
Cm = 10         # capacitance (uF)
tau_m = Rm * Cm  # time constant (msec)
tau_ref = 4     # refractory period (msec)
Vth = 0.5       # spike threshold (V)
V_spike = 1.0   # spike delta (V)


## Input stimulus
I = 1.0         # input current (A)


## Iterate over each time step
function LIFfn(t_rest)
    for (i, t) in enumerate(time)
        if t > t_rest
            Vm[i] = Vm[i-1] + ((-Vm[i-1] + I * Rm) / tau_m) * dt
            if Vm[i] >= Vth
                Vm[i] += V_spike
                t_rest = t + tau_ref
            end
        end
    end
    return Vm
end


## Data structures
Vm = zeros(length(time))  # potential (V) trace over time


## Simulate
V = LIFfn(t_rest)

## Plot membrane potential trace
plot(
    time,
    V,
    title = "Leaky Integrate-and-Fire Example",
    ylabel = "Membrane Potential (V)",
    xlabel = "Time (msec)",
    ylim = [0, 2],
    linewidth = 2,
    color = :green,
    legend = :none,
)
