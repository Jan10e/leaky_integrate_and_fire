# Multi-Quadratic Integrate and Fire (MQIF) neuron
# based on https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=235138&file=/Van_Pottelbergh_2018/MQIF_bistability.py#tabs-2
#
#
# Name: Jantine Broek
# Date: October 2020
#####################################################################
## working directory
cd("/Users/jantinebroek/Documents/03_projects/04_IF/code/julia")

## Loads packages
using Plots         # Select the Plots package for figures
pyplot()            # Plots package will use Pyplot (matplotlib needs to be installed).

## Parameters
C = 1               # capacitance
g_f = 1             # conductance fast ion channels
g_s = 0.2           # conductance slow ion channels
g_u = 0             # conductance ultra-slow ion channels

# input
I_app = 10

# steady state values
V_f0 = -40
V_s0 = -35
V_u0 = -40


## Resets
ΔV_u = 0            # reset addition for ultraslow gating
V_max = 40          # Voltage threshold
V_spike = 80       # spike delta (V)
V_r = -45 #-70        # reset potential fast gating
V_sr = 50 #-30       # reset potential slow gating

## Membrane time constants
τ_s = 10            # time-scale slow
τ_u = 100           # time-scale ultra-slow

## Initial values and time vector
T_final = 200       # msec
dt = 1e-2           # simulation time step
Tt = collect(0:dt:T_final)

# step function
T_step_start = floor(Int, 0.1 * length(Tt))
T_step_stop = floor(Int, 0.7 * length(Tt))
I_step = zeros(length(Tt))
I_step[T_step_start:T_step_stop] .= I_app

## Functions
dV(V, V_s, V_u, V_f0, V_s0, V_u0, g_f, g_s, g_u, C, I) = (g_f*(V - V_f0)^2 - g_s*(V_s - V_s0)^2 - g_u*(V_u - V_u0)^2 + I) / C
dVs(V, V_s, τ_s) = (V - V_s) / τ_s
dVu(V, V_u, τ_u) = (V - V_u) / τ_u

function MQIFfn(V_max, V_r, V_sr, ΔV_u, I_step)

    # initial value
    V = -30
    Vprev = -30
    V_s = -40
    V_u = -40

    # data structures
    spike_train = zeros(length(Tt))
    Vv = zeros(length(Tt))

    # integrate
    for i = 1: length(Tt)-1
        V += dV(V, V_s, V_u, V_f0, V_s0, V_u0, g_f, g_s, g_u, C, I_step[i]) * dt
        V_s += dVs(Vprev, V_s, τ_s) * dt
        V_u += dVu(Vprev, V_u, τ_u) * dt

        # reset
        if V >= V_max
            Vv[i] = V_spike
            V = V_r
            V_s = V_sr
            V_u = V_u + ΔV_u
            spike_train[i] = 1
        else
            Vv[i] = copy(V)
        end

        Vprev = copy(V)
    end

    return Vv, spike_train
end



## Iterate over each time step
@time (V, spk) = MQIFfn(V_max, V_r, V_sr, ΔV_u, I_step)

## Plot input current and membrane potential trace
p1 = plot(
    Tt,
    I_step,
    title = "input current",
    ylabel = "current (A)",
    linewidth = 2,
    color = :black
)

p2 = plot(
    Tt,
    V,
    title = "time trace",
    ylabel = "Membrane Potential (V)",
    xlabel = "time (msec)",
    linewidth = 2,
    legend = :none,
)

## Get rasterplot
spk_dot1 = findall(x -> x == 1, spk[:, 1])
p3 = scatter(
    Tt[spk_dot1, 1],
    zeros(length(spk_dot1), 1),
    title = "scatterplot",
    legend = :none
)

plot(p1, p2, p3, layout = (3,1))
