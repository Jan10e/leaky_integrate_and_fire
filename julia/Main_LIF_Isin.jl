# Leaky integrate and fire neuron
# Sinusoidal input
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
# Vm = zeros(length(time))  # potential (V) trace over time
V_0 = 0.            # init val
Rm = 1              # resistance (kOhm)
Cm = 10             # capacitance (uF)
tau_m = Rm * Cm     # time constant (msec)
tau_ref = 4         # refractory period (msec)
Vth = 0.5           # spike threshold (V)
V_spike = 1.0       # spike delta (V)


## Step function with sinusoidal input
T_step_start = 10
T_step_stop = 40


## Iterate over each time step

function LIFfn(t_rest, I_sin)

    #integrate
    for (i, t) in enumerate(time)
        if t > t_rest
            Vm[i] = Vm[i-1] + ((-Vm[i-1] + I_sin[i-1] * Rm) / tau_m) * dt
            if Vm[i] >= Vth
                Vm[i] += V_spike
                t_rest = t + tau_ref
            end
        end
    end
    return Vm
end


function sinfn(time)
    I_sin = zeros(length(time))

    # step function with sinusoidal input
    for iT = 1:length(time)

        if time[iT] >= T_step_start && time[iT] <= T_step_stop

            # sinusoidal input current (A)
            # I_sin[iT] = 1.0 + sin.(time[iT])
            I_sin[iT] = 1.0 .+ 2.5*sin.(time[iT])     # increase amplitude
            # I_sin[iT] = 1.0 + sin.(3.0*time[iT])    # increase frequency
        end
    end

    return I_sin
end


## Data structures
Vm = zeros(length(time))  # potential (V) trace over time



## Simulate
I_sin = sinfn(time)
Vv = LIFfn(t_rest, I_sin)


## Plot membrane potential trace
p1 = plot(
    time,
    I_sin,
    title = "Sinusoidal step input",
    ylabel = "Input current (A)",
    xlabel = "Time (msec)",
    linewidth = 2,
    legend = :none
)

p2 = plot(
    time,
    Vv,
    title = "LIF with sinusoidal input",
    ylabel = "Membrane Potential (V)",
    xlabel = "Time (msec)",
    ylim = [0, 2],
    linewidth = 2,
    color = :green,
    legend = :none,
)

plot(p1, p2, layout = (2,1))
