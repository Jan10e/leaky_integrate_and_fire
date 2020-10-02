% Coupled IF network with one excitatory (E) and one inhibitory (I) cell
% Based on: https://github.com/faezeamin/NC-Network-of-LIF-Neurons-Coupled-by-Electrical-Synapses
%
%
% Name: Jantine Broek
% Date: September 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program creates a network of n LIF neurons coupled by
% electrical synapses. It solves the following equation to pursue the dynamics
% of each neuron's sub-threshold membrane potential:
%          tau * v_dot = -(v - v_rest) + I_ext - I_synps,
%
% where v is the potential of post synaptic neuron, v_rest is the resting
% potential, I _ext is external input to the neuron which could be from 
% a sensory stimulus or the cumulative activity of neurons around, and
% I_synps is the synaptic input from the neighbour neurons in the network. It
% is calculated using 
%          I_synps = ( V_post - V_pre ) * g_post,pre;
%
% where g is conductivity coefficient of the synapses.
% Since an LIF neuron just depicts the dynamics of sub-threshold membrane 
% potential, this model cannot show the urge of ions during action
% potential.
% Setting the values of g, tau, and V0 in different fashions or randomly 
% results in various fascinating patterns of activity for even a network 
% made up of three neurons.
% In the output, the membrane potential of all the neurons as well as the
% raster plot of their spiking activity is presented.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc; clear all

% parameters LIF
Nn = 2;         % number of neurons
V_th = -65;     % threshold potential
V_rest = -70;   % resting potential

% external forcing
I_ext = 5.5;

% membrane time constant
tau = [1, 1.4];

% simulation time steps
T_final = 150;    % total run time
dt = .01;
n_step = T_final/dt + 1;

% conductance
g = (rand (Nn,Nn) *.1);
g = g - diag(diag(g));
%g=[0,2; 0,0];

% data structures
Tt = zeros(n_step, Nn);
Vv = zeros(n_step, Nn);
Vv(1,1) = -65.1;     % spike thresholod excitatory cell
Vv(1,2) = -70;       % spike threshold inhibitory cell
spike_train = zeros(n_step, Nn);


%% Solve coupled LIF

for i = 1:n_step-1
    for j=1:Nn
        
        I_syn = functions.I_elec_syn(j, i, g, Vv);
        v_a1 = Vv(i,j);  
        [v_a2,spk] = functions.LIF_ODE(V_th, V_rest, tau(j), dt, I_ext, I_syn, v_a1);
        Tt(i+1,j) = Tt(i,j)+dt;
        Vv(i+1,j) = v_a2;
        
        if spk == true
            spike_train(i+1,j) = 1;
        end
        
        v_a1 = v_a2;
    end    
end


%% Plot
figure(1);   
plot(Tt,Vv)
title('Membrane Dynamics')
xlabel('Time')
ylabel('V')

figure(2)
functions.rasterPlot(spike_train, Tt, Nn)
