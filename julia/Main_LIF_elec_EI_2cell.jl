# 2 electrically connected Leaky Integrate and Fire Nn
# Based on this matlab script: https://github.com/faezeamin/NC-Network-of-LIF-Nns-Coupled-by-Electrical-Synapses
#
# It solves the following eqn for the sub-threshold dynamics
# τV̇ = - (V - V_rest) + I_ext - I_syn
#
# V = post-synaptic potential
# V_rest = resting potential
# I_ext = external input to Nn (cumulative activity)
# I_syn = synaptic input from neighbour Nns; where
#   I_syn = (V_post - V_pre) * g_post,pre
#
# Name: Jantine Broek
# Date: Sep 2020
##################################################################
## working directory
cd("/Users/jantinebroek/Documents/03_projects/04_IF/")

## Loads packages
using Plots         # Select the Plots package for figures
pyplot()            # Plots package will use Pyplot (matplotlib needs to be installed).
using LinearAlgebra


## parameters
Nn = 2     # number of cells
V_rest = -70   # resting potential
V_th = -65   # threshold potential


## membrane time constant
tau_E = 1       # excitatory cell Tt constant (msec)
tau_I = 1.4     # inhibitory cell Tt constant (msec)
tau = [tau_E, tau_I]


## external forcing
I_ext = 5.5


## conductivity coefficient of the synapses
g = (rand(Nn, Nn) * 0.1)
g = g - Diagonal(Diagonal(g))


## initial voltage and time vector
T_final = 50     # msec
dt = 1e-2         # simulation Tt step

V_E = -65.1  # spike threshold (V) for excitatory (E) cells
V_I = -70    # spike threshold (V) for inhibitory (I) cells
V_0 = [V_E, V_I]


## Functions
dV(V, V_rest, I_ext, I_syn, tau) = (-(V - V_rest) + I_ext - I_syn) ./ tau

function LIF_ODE(V_th, V_rest, tau, dt, I_ext, I_syn, V_init)

    # initial value
    V = V_init

    # integrate
    V += dV(V, V_rest, I_ext, I_syn, tau) * dt

    # reset
    spk = false
    if V >= V_th
        V = V_rest
        spk = true
    end

    return V, spk
end


function I_syn_elec(iN, iT, g, V)
    I_syn_pre = (V[iT, iN] .- V[iT, :]) .* g[iN, :]
    I_syn_pre = sum(I_syn_pre, dims = 1)
    I_syn_elec = I_syn_pre[1]
    return I_syn_elec
end


## Data structures

Tt = collect(0:dt:T_final)
Vv = zeros(length(Tt), Nn)
Vv[1, :] = V_0
spike_train = zeros(length(Tt), Nn)


## Solve for multiple neurons

for iT = 1:length(Tt)-1
    for iN = 1:Nn
        I_syn = I_syn_elec(iN, iT, g, Vv)
        (Vv[iT+1, iN], spk) =
            LIF_ODE(V_th, V_rest, tau[iN], dt, I_ext, I_syn, Vv[iT, iN])

        if spk == true
            spike_train[iT+1, iN] = 1
        end

    end
end


## Plot membrane potential trace
p1 = plot(
    Tt,
    Vv,
    title = "Membrane Dynamics",
    ylabel = "Membrane Potential (V)",
    xlabel = "time (msec)",
    linewidth = 2,
    legend = :none,
)


## Get rasterplot
spk_E = findall(x -> x == 1, spike_train[:, 1])
p2 = scatter(Tt[spk_E, 1], zeros(length(spk_E), 1), label = "E")
spk_I = findall(x -> x == 1, spike_train[:, 2])
scatter!(
    Tt[spk_I, 1],
    zeros(length(spk_I), 1) .+ 1,
    ylim = [0, 10],
    title = "Raster Plot",
    ylabel = "Nn number",
    xlabel = "time (msec)",
    label = "I",
)

plot(p1, p2, layout = (2,1))
