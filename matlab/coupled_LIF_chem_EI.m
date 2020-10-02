% Coupled IF network with one excitatory (E) and one inhibitory (I) cell
% Based on: https://github.com/faezeamin/NC-Two-LIF-Neurons-Chemical-Synapse/blob/master/Two_LIF_chemSyn.m
%
%
% Name: Jantine Broek
% Date: September 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is generating two LIF neurons which are connected through a
% chemical synapse. It solves LIF equation as follows:

% tau * v_dot = -(v - v_rest) + I_ext - I_syn,

% where I_ext is the external input, which can be a sensory stimulus or any
% input from the population activity of other neurons. v_rest is the
% membrane potential of resting state. I_syn is accumulative synaptic inputs
% arrived from all the pre-synaptic neurons to the given post-synaptic one.
% v_th is the threshold potential. tau is the time constant of the neurons.

% alpha and beta and N specify the shape of S.

% Note that negative I_syn corresponds to EPSP and positive I_syn corresponds to IPSP
% (Consider the negative sign before I_syn in the LIF equation).
%%
clear; clc

% Parameters
Nn = 2;
V_th = -65;
V_rest = -70;

% Membrane time constant
%tau = rand (no_neurons,1)*10+10;
tau = [10,12];

% External forcing
I_ext = [V_th - V_rest + .01, V_th - V_rest + .01];

% Conductivity matrix.
g = [0, .01; .02, 0]; % its elements demonstrate maximal conductance of the synapses.

%% Initial voltage and time vector
T_final = 150;
dt = .01;
n_tSteps = T_final/dt +1;

V_E = -70;
V_I = -65.1;


%% Rectangular pulse
rec_bound = zeros(Nn, Nn, n_tSteps);     % rec_bound represents the fraction of bound receptors.
N = zeros(n_tSteps,Nn);
I_syn = zeros(n_tSteps,Nn);


%% E/I characterisation
% Characterisation of post-synaptic potential.
alpha = [.6, 1.5];  % first element is negative I_syn (EPSP), second element is the positive I_syn (IPSP).
beta = [.3, .3];    % first element is negative I_syn (EPSP), second element is the positive I_syn (IPSP).
N_duration = 1;     % N is the concentration of neurotransmitters in synaptic cleft.

% Whether a synapse is excitatory or inhibitory is specified using E_syn, the synaptic reversal potential (syn_rev).
% For inhibitory synapses, syn_rev equals to syn_rev_inh, and for excitatory ones the value is set to zero.
% syn_rev = zeros(2,2);
syn_rev_inh = -80;
syn_rev = [0, 0; syn_rev_inh, 0];


%% Data Structures

T = zeros(n_tSteps,1);
t = 0:dt:T_final;
V = zeros(n_tSteps,Nn);
% V(1,1:Nn) = rand(1,Nn)*4 +V_rest;
V(1, 1:2) = [V_E, V_I];
spike_train = zeros(n_tSteps,Nn);


%% Solve for multiple neurons

for iT = 1:n_tSteps -1
    
    for iN = 1:Nn
        v_a1 = V(iT,iN);
        
        % get chemical input
        [I_chem, s] = functions.I_chem_synps(iN,iT,g,rec_bound,N,alpha,beta,dt,syn_rev,v_a1,syn_rev_inh);
        I_syn(iT,iN) = I_chem;
        rec_bound(iN,:,iT+1) = s;
        
        % solve
        [v_a2, spk] = functions.LIF_ODE(V_th, V_rest, tau(iN), dt, I_ext(iN), I_chem, v_a1 );
        V(iT+1,iN) = v_a2;
        
        if spk == true
            spike_train(iT,iN) = 1;
            n = rectpuls(t-T(iT)-.5*N_duration,N_duration);
            n = n';
            N(:,iN) = N(:,iN)+n;
        end
        
    end
    T(iT+1) = T(iT)+dt;
    
end


%% Plots
% Figure(1) depicts the sub-threshold dynamics of neurons' membrane potential.
figure(1);
plot(T,V)
title(' Dynamics of Membrane Potential')
xlabel('Time')
ylabel('V')
legend ('Neuron 1' , 'Neuron 2', 'location', 'east')

% Figure(2) is raster plot of the two neurons.
figure(2)
functions.rasterPlot(spike_train,T,Nn)

% Figure(3) represents the dynamics of N and S, which determine the input I_syn. You can produce EPSP and IPSP
% corresponding to NMDA, GABAa, GABAb, etc. by setting relevant values for alpha and beta and
% therefore for the shape the S profile.
figure (3)
p1 = plot (T,N(:,1));
hold on;
p2 = plot (T,N(:,2));
ss(:,1) = rec_bound(1,2,:);
p3 = plot (T,ss);
hold on
ss(:,1) = rec_bound(2,1,:);
p4 = plot (T,ss);
title('Model Prameters: Concentration of Transmitters, [N], & Fraction of Bound Receptors, S')
xlabel('Time')
ylabel('[N] & S')
legend([p1 p2 p3 p4],{'[N_1]', '[N_2]', 'S_{1,2}', 'S_{2,1}'})

% Figure(4) demonstrates the accumulative synaptic inputs from pre-synaptic neurons to the post-synaptic one
figure(4)
plot(T,I_syn)
hold on
%plot(T,I_ext)
% plot(T,-I_synps+I_ext)
title('Synaptic Input')
xlabel('Time')
ylabel('I_{syn}')
legend('I_{syn,1}', 'I_{syn,2}')



% figure (5)
% plot (T,spike_train)
% title('Spike Train')





