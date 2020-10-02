# 2 chemically connected Leaky Integrate and Fire neurons
# Based on this matlab script: https://github.com/faezeamin/NC-Two-LIF-Neurons-Chemical-Synapse/blob/master/Two_LIF_chemSyn.m
#
# It solves the following eqn for the sub-threshold dynamics
#
# τV̇ = -(V - V_rest) + I_ext - I_syn,
#
# V = post-synaptic potential
# V_rest = resting potential
# I_ext = external input to Nn (cumulative activity)
# I_syn = synaptic input from neighbour Nns; where
#         negative I_syn corresponds to EPSP and positive I_syn corresponds to IPSP
#
# The input is a rectangular pulse, where alpha and beta and N specify the shape.
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
tau_E = 10       # excitatory cell Tt constant (msec)
tau_I = 12       # inhibitory cell Tt constant (msec)
tau = [tau_E, tau_I]


## external forcing
# I_ext = 5.5
I_ext = [V_th - V_rest + .01, V_th - V_rest + .01]


## conductivity coefficient of the synapses
# g = (rand(Nn, Nn) * 0.1)
# g = g - Diagonal(Diagonal(g))
g = [0 .01; .02 0]


## initial voltage and time vector
T_final = 150     # msec
dt = 1e-2         # simulation Tt step

V_E = -65.1  # spike threshold (V) for excitatory (E) cells
V_I = -70    # spike threshold (V) for inhibitory (I) cells
V_0 = [V_E, V_I]


## E/I characterisation
# Characterisation of post-synaptic potential.
α = [.6, 1.5]  # first element is negative I_syn (EPSP), second element is the positive I_syn (IPSP).
β = [.3, .3]    # first element is negative I_syn (EPSP), second element is the positive I_syn (IPSP).
N_duration = 1     # N is the concentration of neurotransmitters in synaptic cleft.

# Whether a synapse is excitatory or inhibitory is specified using E_syn, the synaptic reversal potential (syn_rev).
# For inhibitory synapses, syn_rev equals to syn_rev_inh, and for excitatory ones the value is set to zero.
# syn_rev = zeros(2,2);
syn_rev_inh = -80
syn_rev = [0 0; syn_rev_inh 0]


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


function I_syn_chem(iN, iT, g, S, N, α, β, dt, E_syn, V, E_syn_inh)
    n = size(g, 1)
    s = zeros(1, n)

    I_syn = g[iN,:] .* S[iN,:,iT] .* (V - E_syn[iN,:])
    I_syn = sum(I_syn, dims = 2)

    for i = 1:n

        if E_syn[iN, i] == E_syn_inh
            α_1 = α[2]
            β_1 = β[2]
        else
            α_1 = α[1]
            β_1 = β[2]
        end

        s[i] = S[iN,i,iT] + dt * (α_1 * N[iT, i] * (1 - S[iN,i,iT]) - β_1 * S[iN,i,iT])
    end

    return I_syn, s
end


## Data structures

Tt = collect(0:dt:T_final)
Vv = zeros(length(Tt), Nn)
Vv[1, :] = V_0
spike_train = zeros(length(Tt), Nn)


## Rectangular pulse
rec_bound = zeros(Nn, Nn, length(Tt))     # rec_bound represents the fraction of bound receptors.
N = zeros(length(Tt), Nn)
I_syn = zeros(length(Tt), Nn)


## Solve for multiple neurons

for iT = 1:length(Tt)-1
    for iN = 1:Nn

        # get chemical input
        (I_chem, s) =
            I_syn_chem(iN, iT, g, rec_bound, N, α, β, dt, syn_rev, Vv, syn_rev_inh)
        I_syn[iT, iN] = I_chem
        rec_bound[iN, :, iT + 1] = s

        # solve with given input
        (Vv[iT+1, iN], spk) =
            LIF_ODE(V_th, V_rest, tau[iN], dt, I_ext[iN], I_chem, Vv[iT, iN])

        if spk == true
            spike_train[iT+1, iN] = 1
            # n = rectpuls[t - T[iT] - .5 * N_duration, N_duration]
            n = rectangle(t - T[iT] - .5 * N_duration, N_duration)
            n = n'
            N[:,iN] = N[:,iN] + n
        end

    end
end


## Plot membrane potential trace
plot(
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
scatter(Tt[spk_E, 1], zeros(length(spk_E), 1), label = "E")
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
