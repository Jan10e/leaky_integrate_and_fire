# Multi-Quadratic Integrate and Fire (MQIF) neuron
# based on https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=235138&file=/Van_Pottelbergh_2018/MQIF_bistability.py#tabs-2
#
#
# Name: Jantine Broek
# Date: October 2020
#####################################################################
## working directory
cd("/Users/jantinebroek/Documents/03_projects/04_IF/")

## Loads packages
using Plots         # Select the Plots package for figures
pyplot()            # Plots package will use Pyplot (matplotlib needs to be installed).

## Parameters
C = 1               # capacitance
g_f = 1             # conductance fast ion channels
g_s = 0.2           # conductance slow ion channels
g_u = 0             # conductance ultra-slow ion channels
V_th = -30          # Voltage threshold
# V_spike = 1.0       # spike delta (V)
I = 5               # input

## Resets
V_rest = -70        # resting potential fast gating
V_srest = -30       # resting potential slow gating
ΔV_u = 0            # reset addition for ultraslow gating

## Membrane time constants
τ_s = 10            # time-scale slow
τ_u = 100           # time-scale ultra-slow

## Steady state values and time vector
V_f0 = -70
V_s0 = -35
V_u0 = -40

T_final = 20       # msec
dt = 1e-2           # simulation time step
Tt = collect(0:dt:T_final)

## Functions
dV(V, V_s, V_u, V_f0, V_s0, V_u0, g_f, g_s, g_u, C, I) = (g_f*(V - V_f0)^2 - g_s*(V_s - V_s0)^2 - g_u*(V_u - V_u0)^2 + I) / C
dVs(V, V_s, τ_s) = (V - V_s0) / τ_s
dVu(V, V_u, τ_u) = (V - V_u0) / τ_u

function MQIFfn(V_th, V_rest, V_srest, ΔV_u)

    # initial value
    V = V_rest
    Vprev = V_rest
    V_s = V_srest
    V_u = V_srest

    # data structures
    spike_train = zeros(length(Tt))
    Vv = zeros(length(Tt))

    # integrate
    for i = 1: length(Tt)-1
        V += dV(V, V_s, V_u, V_f0, V_s0, V_u0, g_f, g_s, g_u, C, I) * dt
        V_s += dVs(Vprev, V_s, τ_s) * dt
        V_u += dVu(Vprev, V_u, τ_u) * dt

        # reset
        if V >= V_th
            V = V_rest
            V_s = V_srest
            V_u = V_u + ΔV_u
            spike_train[i] = 1
        end

        Vprev = copy(V)
        Vv[i] = copy(V)
    end

    return Vv, spike_train
end


## Iterate over each time step
@time (V, spk) = MQIFfn(V_th, V_rest, V_srest, ΔV_u)

## Plot membrane potential trace
p1 = plot(
    Tt,
    V,
    title = "Membrane Dynamics",
    ylabel = "Membrane Potential (V)",
    xlabel = "time (msec)",
    linewidth = 2,
    legend = :none,
)

## Get rasterplot
spk_dot1 = findall(x -> x == 1, spk[:, 1])
p2 = scatter(Tt[spk_dot1, 1], zeros(length(spk_dot1), 1), legend = :none)

plot(p1, p2, layout = (2,1))
