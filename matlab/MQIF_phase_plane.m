% Title: Phase plane of MQIF model
%
% Author: Jantine Broek
% Created: October 2020
%% 
clear
clc
dbstop if error

%% Paths
dir_base = '/Users/jantinebroek/Documents/03_projects/04_IF/code';

dir_work = '/matlab';
dir_data = '/data';
dir_fig = '/figures';

cd(fullfile(dir_base, dir_work));

%% Options

options = optimoptions('fsolve','MaxFunEvals',1e16,'MaxIter',...
    1e16,'TolFun',1e-14,'TolX',1e-14,'Display', 'off');

%% Functions and Parameters
C = 1;               % capacitance
g_f = 1;             % conductance fast ion channels

% membrane time constants
tau_s = 10;          % time-scale s  

% steady-state value
v_f0 = -40;         % fast time-scale
v_s0 = -35;         % slow time-scale

% reset values
v_max = -30;         % spike cutoff
v_spike = 80;       % approx of spike


%% Parameters for different kinds of spiking

% bistability
g_s = 0.2;        % conductance slow ion channels
v_r = -40;        % voltage reset
v_sr = -30;       % reset for slow gating  
I_app = 50;       % input

% spike latency
% g_s = 0.5;        % conductance slow ion channels
% v_r = -45;        % voltage reset
% v_sr = 0;         % reset for slow gating  
% I_app = 10;       % input
    

%% Time and Initial values

% time span and step (ms)
T_final = 50; 
dt = 1e-2;                          
Tt = 0:dt:T_final;

% initial values
v_init = -80;         
u_init = -40;


% pulse of input DC current
I=[zeros(1,floor(0.1*length(Tt))), I_app * ones(1,floor(0.7*length(Tt))), ...
    zeros(1,floor(0.2*length(Tt))), 0];    

I_0 = zeros(1, length(Tt)); 


%% Integrate: forward Euler

[Vv, Uu] = functions.MQIF_forward_euler(g_f, g_s, v_f0, ...
    v_s0, tau_s, C, v_r, v_sr, I, v_max, v_spike, v_init, u_init, Tt); 

[Vv0, Uu0] = functions.MQIF_forward_euler(g_f, g_s, v_f0, ...
    v_s0, tau_s, C, v_r, v_sr, I_0, v_max, v_spike, v_init, u_init, Tt); 

%% Plot time trace
f1 = figure(1);

% plot input
subplot(4,1,1)
plot(Tt, I, 'k', 'LineWidth', 3)
title('Input')

% plot time trace
subplot(4,1,2)
plot(Tt, Vv, 'b', 'LineWidth', 2); 
title('Time trace')


%% ____________ Plots Phase Plane ________________%
subplot(4,1,3)

Iapp_0 = 0;

syms vv1 vv2
dot_vv_0 = (g_f*(vv1 - v_f0)^2 - g_s * (vv2 - v_s0)^2 + Iapp_0) / C;
dot_uu = (vv1 - vv2) / tau_s;

hold on

% plot V-nullcline for I =0
nullcline_v_0 = ezplot(dot_vv_0, [-80 80]);
set(nullcline_v_0, 'color',[0, 0.6, 0.748], 'Linewidth', 2);

% plot u-nullcline
nullcline_u = ezplot(dot_uu, [-80 80]);
set(nullcline_u, 'color', [0.4660, 0.6740, 0.1880], 'Linewidth', 2);

% trajectory for no input
plot(Vv0, Uu0, 'color',[1 0.08 0.57],'LineWidth', 1) 
hold off

xlabel('membrane potential V', 'interpreter','latex','fontsize',11);
ylabel('recovery variable u', 'interpreter','latex','fontsize',11);
title('Phase Plane, with I = 0', 'fontsize', 12);


subplot(4,1,4)
Iapp_up = I_app;

hold on
dot_vv_Iapp = (g_f*(vv1 - v_f0)^2 - g_s * (vv2 - v_s0)^2 + Iapp_up) / C;

% plot V-nullcline for I = I_app
nullcline_v_Iapp = ezplot(dot_vv_Iapp, [-80 80]);
set(nullcline_v_Iapp, 'color',[0, 0.6, 0.748], 'Linewidth', 2);

% plot u-nullcline
nullcline_u = ezplot(dot_uu, [-80 80]);
set(nullcline_u, 'color', [0.4660, 0.6740, 0.1880], 'Linewidth', 2);

% trajectory
plot(Vv, Uu, 'color',[1 0.08 0.57],'LineWidth', 1) 

xlabel('membrane potential V', 'interpreter','latex','fontsize',11);
ylabel('recovery variable u', 'interpreter','latex','fontsize',11);
title('Phase Plane, with I = I_{app}', 'fontsize', 12);




%% Export/Save
outfile = ['MQIF_phase_plane', '_Vr', num2str(v_r), '_Vsr',...
           num2str(v_sr), '_I', num2str(I_app)];
       
suffix_fig = '';
suffix_data = '';       

out_mat = [outfile, suffix_data, '.mat'];
out_fig_png = [outfile, suffix_fig, '.png'];
out_fig_eps = [outfile, suffix_fig, '.eps'];

outpath_data = fullfile(dir_base, dir_data, out_mat);
outpath_fig_png = fullfile(dir_base, dir_fig, out_fig_png);
outpath_fig_eps = fullfile(dir_base, dir_fig, out_fig_eps);

% figures
% saveas(f1, outpath_fig_eps,'eps')
% saveas(f1, outpath_fig_png,'png')

