% Coupled IF network with one excitatory (E) and one inhibitory (I) cell
% Based on: https://github.com/faezeamin/NC-Network-of-LIF-Neurons-Coupled-by-Electrical-Synapses
%
%
% Name: Jantine Broek
% Date: September 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This program generates a network of n LIF neurons which are connected 
%  through chemical synapses. It solves LIF equation as follows:

%     tau * v_dot = -(v - v_rest) + I_ext - I_syn,

%  where tau is the time constant of the neurons, v_rest is the resting 
%  state membrane potential. I_ext is the external input, which is 
%  equivalent to a sensory stimulus or any input sourced from the activity
%  of population of neurons nearby. I_syn is accumulative synaptic inputs
%  arrived from all the pre-synaptic neurons to the given post-synaptic one
%  following the arrival of action potential. It is calculated through

%      I_syn(t) = sum( g_(ij) * S_(ij)(t) * ( v(t)-E_syn(ij) ) )

%  where i is associated index of pre synaptic neurons, and j is for 
%  post-synaptic neuron. v is the postsynaptic 
%  potential, and E_syn is the synaptic  reversal potential. Whether a 
%  synapse is excitatory or inhibitory is determined by E_syn. For 
%  inhibitory synapses, E_syn equals to E_syn_inh, and for excitatory 
%  ones the value is set to zero.
%  g is conductivity matrix. its elements demonstrate maximal conductance 
%  of the synapse. S represents the fraction of 
%  bound receptors. 
% 
%            S_dot(t) = alpha * N * (1-S(t)) - beta * S(t)
            
%  N is equivalent to the concentration of transmitters released into the 
%  synaptic cleft following the arrival of an action potential at the 
%  pre-synaptic terminal. Since the concentration of neurotransmitters in 
%  the cleft rises and falls very rapidly, it is assumed that N occurs as
%  a rectangular pulse. alpha and beta are the rate constants for 
%  transmitter binding to and unbinding from post-synaptic receptor 
%  respectively. alpha and beta and duration of N characterize the 
%  dynamics of post-synaptic potential (EPSP & IPSP). They directly specify 
%  the shape of S (Destexhe, 1993).  The first elements of alpha and beta
%  correspond to EPSP, and the second ones correspond to IPSP.

%  In this program you can specify the percentage of connectivities in the 
%  neuronal graph. So the network can be fully conected 
%  (Connectivity_rate=1), or each neuron is connected to 10% of the other
%  neurons as it is in the cortical networks. Matrix a determined the
%  synaptic connectivity pattern. Among the connected neurons, the 
%  proportion of inhibitory to excitatory neurons can be set using
%  inh2exc_rate. In cortex this value is said to be 1 to 4. 
 
% This program returns 4 graphs. Figure(1) depicts the sub-threshold
% dynamics of neurons' membrane potential. Figure(2) is raster plot of the
% two neurons. Figure(3) represents the dynamics of N and S,
% which determine the input I_syn. You can produce EPSP and IPSP
% corresponding to NMDA, GABAa, GABAb, etc. by setting relevant values for 
% for alpha and beta and therefore for the shape the S profile.
% Figure(4) demonstrates the accumulative synaptic inputs from pre-synaptic 
% neurons to the post-synaptic one. Note that negative I_syn corresponds 
% to EPSP and positive I_syn corresponds to IPSP (Consider the negative 
% sign before I_syn in the LIF equation). 

% Using the given values, we can see a synchronous activity in a network of 
% 50 neurons. You can see that the frequency of oscillation in this network 
% (given the current set values) is approximately 26 Hz.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc; clear all

% parameters LIF
Nn = 20;        % number of neurons
V_th = -65;     % threshold potential
V_rest = -70;   % resting potential

% membrane time constant
tau = rand (Nn,1)*10+10; 

% simulation time steps
T_final = 500;
dt = .01;
Tt = 0:dt:T_final;
n_tSteps = T_final/dt+1;

% external forcing
I_ext = .1 * abs( randn(Nn,1) ) + (V_th - V_rest);

% other
N_duration = 1;

alpha = [.6,1.5];
beta = [.3,.3];

Connectivity_rate = 1;
a = rand(Nn,Nn);
a (a>Connectivity_rate) = 3;
a (a<Connectivity_rate) = 1;
a (a==3) = 0;

inh2exc_rate = 1/3.5;
inh2all = 1/(1+1/inh2exc_rate);
E_syn_inh = -80;
E_syn = zeros(Nn,Nn);
E_syn(a==1) = rand (sum(sum(a),2),1);
E_syn (E_syn>inh2all) =0;
E_syn (E_syn~=0) = E_syn_inh;


% conductance
g_mean_exc = 0.005;
%g_mean_inh = g_mean_exc / inh2exc_rate ; %To satisfy exc-inh balance of the network 
g_mean_inh = .005;
g = zeros(Nn,Nn);
g(E_syn == E_syn_inh) = ...
    .001 * abs (randn(sum(sum(E_syn==E_syn_inh),2),1) ) + g_mean_inh;
g(a==1 & E_syn==0) = ...
    .001 * abs (randn(sum(sum(a),2) - sum(sum(E_syn==E_syn_inh),2),1) ) + g_mean_exc;
g = g - diag(diag(g));



% data structures
spike_train = zeros(n_tSteps,Nn);
T = zeros(n_tSteps,1);
S = zeros(Nn, Nn,n_tSteps);
N = zeros(n_tSteps,Nn);
I_synps = zeros(n_tSteps,Nn);
V = zeros(n_tSteps,Nn);
V(1,1:Nn) = rand(1,Nn)*(V_th - V_rest) + V_rest;


%% Solve coupled LIF

for tStep=1:n_tSteps -1
    
    for j=1:Nn
        v_a1 = V(tStep,j);
        
        [i_synps,s] = functions.I_chem_synps(j,tStep,g,S,N,alpha,beta,dt,E_syn,v_a1,E_syn_inh);
        I_synps (tStep,j) = i_synps;
        S(j,:,tStep+1) = s;
        
        [v_a2,spk] = functions.LIF_ODE(V_th, V_rest, tau(j), dt, I_ext(j), i_synps, v_a1 );
        V(tStep+1,j) = v_a2;
        
        if spk == true
            spike_train(tStep,j) = 1;
            n = rectpuls(Tt-T(tStep)-.5*N_duration,N_duration);
            n = n';
            N(:,j) = N(:,j)+n;
        end
         
    end    
    T(tStep+1) = T(tStep)+dt;
  
end


%% Plots
figure(1);  
plot(T,V)
title(' Dynamics of Membrane Potential')
xlabel('Time')
ylabel('V')

figure (2)
plot (T,N)
hold on;
ss(:,1) = S(1,2,:);
plot (T,ss)
hold on
ss(:,1) = S(2,1,:);
plot (T,ss)
title('Model Prameters: Concentration of Transmitters, [N], & Fraction of Bound Receptors, S')
xlabel('Time')
ylabel('[N] & S')


figure(3)
plot(T,I_synps)
%hold on 
%plot(T,I_ext)
%plot(T,-I_synps+I_ext)
title('Synaptic Input') 
xlabel('Time')
ylabel('I_{syn}')


figure(4)
functions.rasterPlot(spike_train,T,Nn)

% figure (5)
% plot (T,spike_train)
% title('Spike Train')