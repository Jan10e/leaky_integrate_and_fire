% Examples using SimLIFNet, to look at architecture of LIF connected cells
%
% Examples from comments in function
%
% Name: Jantine Broek
% Date: October 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% === Examples ===

% 1) Reproduce figure 1 from Lewis & Rinzel 2003:
% Two neurons coupled with inhibition and synaptic strength of -0.2
W = [0 -0.2; -0.2 0];     

% Asynchronous behaviour at low driving currents (1.1) and synaptic
% densities of 3, and initial membrane potentials of 0.4 and 0
[spkFull, NetParams, V] = functions.SimLIFNet(W,'simTime',35,'tstep',1e-2,... 
           'offsetCurrents',1.1*[1; 1],'synapticDensity',[3; 3],...
           'initialConditions',[0.4; 0]);
% Now observe synchronous behavior with stronger forcing, 1.6
[spkFull, NetParams, V] = functions.SimLIFNet(W,'simTime',35,'tstep',1e-2,... 
       'offsetCurrents',1.6*[1; 1],'synapticDensity',[3; 3],...
       'initialConditions',[0.4; 0]);

   
% 2) Investigate random network architectures
% A random sparse network with noise
W = 1*rand(8)-0.5;
W(randperm(numel(W))>round(0.25*numel(W))) = 0;
[spk, NetParams, V] = functions.SimLIFNet(W,'simTime',100,'tstep',1e-2,...
   'offsetCurrents',0.8*ones(length(W),1),...
   'noiseAmplitude',0.2*ones(length(W),1));

% A deterministic medium-sized highly inhibitory network initialized
% randomly
W = log(abs(randn(12)));
[spk, NetParams, V] = functions.SimLIFNet(W,'simTime',35,'tstep',1e-2,...
   'offsetCurrents',1.1*ones(length(W),1));


% 3) An example using a neuron with a refractory period:
% All neurons will drive neuron 3 at high rate, but neuron 3 will have a
% long refractory period
W = [0 0 0.5; 0 0 0.5; 0 0 0.2];
[spk, NetParams, V] = functions.SimLIFNet(W,'simTime',50,'tstep',1e-2,...
     'offsetCurrents',[1.5 1.5 0.5]','refractoryTime',[0 0 5]',...
     'initialConditions',[0.2 0.4 0]');

 
% 4) Examine effects of synaptic density:
% Connect a driver neuron (1) to others via varied synaptic densities
W = [0 ones(1,3)*0.4; zeros(3,4)];
[spk, NetParams, V] = functions.SimLIFNet(W,'simTime',35,'tstep',1e-2,...
      'offsetCurrents',[1.6 0.6 0.6 0.6]',...
      'synapticDensity',[0 1 4 7; zeros(3,4)]);

  
% 5) Apply forcing functions:
% Neuron 1 excites 3, neuron 2 inhibits 3.
W = [0 0 0.5; 0 0 -0.5; 0 0 0];

% Neuron 1 has constant forcing starting at t=10, neuron 2 has sinusoidal
% forcing throughout
Ffcns = {@(t) 1.5*heaviside(t-10), 1; @sin, 2};
[spk, NetParams, V] = functions.SimLIFNet(W,'simTime',35,'tstep',1e-2,...
       'offsetCurrents',[0.8 0.8 0.8]','forcingFunctions',Ffcns);

   
% 6) Forcing with noise:
% Apply a constant current to a single neuron from t=10 to 50 in the
% presence of noise (heaviside is slow for larger networks)
[spk, NetParams, V] = functions.SimLIFNet(0,'forcingFunctions', ... 
      {@(t) 1.6*(heaviside(t-10)-heaviside(t-50)), 1}, ... 
      'noiseAmplitude',0.3);

  
% 7) Building logical operators with LIF networks (exploiting adjustable
% synaptic densities):
% Build a network with modular logical elements. The last (10th) neuron
% will serve as the output. Neurons 1 and 2 are inputs to an OR operator;
% neurons 3 and 4 are inputs to an AND operator; neuron 6 is the input to a
% NOT opperator while input 5 is an interneuron providing chronic
% stimulation to the output neuron, neurons 8 and 9 are inputs to an XOR
% operator with slow synaptic connections while neuron 7 is a fast
% inhibitory interneuron.
W = zeros(10);
W([99 98 97 96 95 94 93 92 91 69 68])=[3 3 -15 -1.5 2 1 1 2 2 1 1];
A = 4*ones(10); A([69 68 99 98 97])=[8 8 1 1 8];
F = {@(u) 1.4*(ge(u,5)-ge(u,15)+ge(u,35)-ge(u,45)),1; ...
      @(u) 1.4*(ge(u,20)-ge(u,30)+ge(u,35)-ge(u,45)),2; ...
      @(u) 1.4*(ge(u,55)-ge(u,65)+ge(u,70)-ge(u,80)),3; ...
      @(u) 1.4*(ge(u,70)-ge(u,80)+ge(u,85)-ge(u,95)),4; ...
      @(u) 1.4*(ge(u,100)-ge(u,140)),5; ...
      @(u) 1.4*(ge(u,120)-ge(u,130)),6; ...
      @(u) 1.4*(ge(u,160)-ge(u,170)+ge(u,190)-ge(u,210)),8; ...
      @(u) 1.4*(ge(u,175)-ge(u,185)+ge(u,190)-ge(u,210)),9};
[spk, NetParams, V] = functions.SimLIFNet(W,'simTime',220,'tstep',1e-2, ... 
    'forcingFunctions',F,'refractoryTime',ones(length(W),1)/2, ...
    'synapticDensity',A);
set(NetParams.handles.sub1,'ylim',[-0.4 1.9])


% 8) Errors due to time-step sizes using multi-part forcing
% Build a network with complex forcing
W = [0 -0.2; -0.1 0];
F = { @(u) (u>3 & u<5)*0.8 + (u>=5 & u<10.5)*sin(2*pi/2*u) + ...
          (u>=10.5)*exp(-(u-10.5)/3),1};

% Run the simulation using a range of time steps.
dt = [logspace(-4,-1,4) 0.5];
V = cell(size(dt));
pr = [1 zeros(1,length(dt)-1)];
for k=1:length(dt)
   [spk, NetParams, V{k}] = functions.SimLIFNet(W,'simTime',25,'tstep',dt(k), ... 
   'forcingFunctions',F,'offsetCurrents',[1; 1.2],'plotResults',pr(k));
end


% Plot the difference between the the smallest time step simulation and
% each of the other larger time step simulations (the error) for the second
% neuron in the network.
figure; hold on
cmp = hsv(length(dt));
for k=1:length(dt)
    plot( linspace(0,25,length(V{k})), V{1}(2,1:dt(k)/dt(1):end)-V{k}(2,:),'color',cmp(k,:),'linewidth',2 )
end
ylabel('Difference from Smallest Timestep'); xlabel('time')
legend(num2str(dt'),'location','SW')      



