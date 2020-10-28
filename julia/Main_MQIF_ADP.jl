# Multi-Quadratic Integrate and Fire (MQIF) neuron
# based on https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=235138&file=/Van_Pottelbergh_2018/MQIF_ADP.py#tabs-2
#
# ADP
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
g_s = 0.5           # conductance slow ion channels
g_u = 0             # conductance ultra-slow ion channels

# input
I_app_up1 = 3
I_app_up2 = 5

# steady state values
V_f0 = -40
V_s0 = -39
V_u0 = -40


## Resets
ΔV_u = 0            # reset addition for ultraslow gating
V_max = -30          # Voltage threshold
V_spike = 0       # spike delta (V)

V_r = -40           # reset potential fast gating
V_sr = -35          # reset potential slow gating

## Membrane time constants
τ_s = 10            # time-scale slow
τ_u = 100           # time-scale ultra-slow

## Initial values and time vector
T_final = 500       # msec
dt = 1e-2           # simulation time step
Tt = collect(0:dt:T_final)

str_input = "Iup1_$I_app_up1 Iup2_$I_app_up2 Vr$V_r Vsr$V_sr"

# step function
T_step_start1 = floor(Int, 0.1 * length(Tt))
T_step_stop1 = floor(Int, 0.13 * length(Tt))
T_step_start2 = floor(Int, 0.3 * length(Tt))
T_step_stop2 = floor(Int, 0.38 * length(Tt))
T_step_start3 = floor(Int, 0.5 * length(Tt))
T_step_stop3 = floor(Int, 0.53 * length(Tt))
T_step_start4 = floor(Int, 0.8 * length(Tt))
T_step_stop4 = floor(Int, 0.88 * length(Tt))

I_step = zeros(length(Tt))
I_step[T_step_start1:T_step_stop1] .= I_app_up1
I_step[T_step_start2:T_step_stop2] .= I_app_up1
I_step[T_step_start3:T_step_stop3] .= I_app_up2
I_step[T_step_start4:T_step_stop4] .= I_app_up2

## Functions
dV(V, V_s, V_u, V_f0, V_s0, V_u0, g_f, g_s, g_u, C, I) = (g_f*(V - V_f0)^2 - g_s*(V_s - V_s0)^2 - g_u*(V_u - V_u0)^2 + I) / C
dVs(V, V_s, τ_s) = (V - V_s) / τ_s
dVu(V, V_u, τ_u) = (V - V_u) / τ_u

function MQIFfn(V_max, V_r, V_sr, ΔV_u, I_step)

    # initial value
    V = -41.5
    Vprev = -41.5
    V_s = -41.5
    V_u = -41.5

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
    color = :black,
    legend = :none
)

p2 = plot(
    Tt,
    V,
    title = "time trace",
    ylabel = "Membrane Potential (V)",
    xlabel = "time (msec)",
    linewidth = 1,
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

#save figures
cd("/Users/jantinebroek/Documents/03_projects/04_IF/code/figures")
savefig("MQIF_ADP_" * str_input * ".eps")
