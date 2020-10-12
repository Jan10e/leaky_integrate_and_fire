% Title: Phase plane of MQIF model
%
% Author: Jantine Broek
% Created: October 2020
%% 
clear
clc
dbstop if error

%% Paths
dir_base = '/Users/jantinebroek/Documents/03_projects/04_IF';

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
g_s = 0.2;           % conductance slow ion channels

% membrane time constants
tau_s = 10;          % time-scale s  

% steady-state value
v_f0 = -70;         % fast time-scale
v_s0 = -35;         % slow time-scale

% resets
v_max = -30;        % spike cutoff
v_spike = 50;       % approx of spike

v_r = -70;          % voltage reset
v_sr = -30;         % reset for slow gating  

% input
I_app = 10;
    

%% Time and Initial values

% time span and step (ms)
T = 10; 
dt = 1e-2;                          
n = round(T/dt);        % number of simulation steps

% initial values
v = v_f0 * ones(1,n); 
% u=0*v;
u = v_s0 * ones(1,n);

% pulse of input DC current
I=[zeros(1,0.1*n), I_app * ones(1,0.9*n)];    


%% Integrate: forward Euler

% forward Euler method
for i = 1:n-1                         
    v(i+1) = v(i) + dt * ((g_f * (v(i) - v_f0)^2 - g_s * (u(i)-v_s0)^2 + I(i)) / C);
    u(i+1) = u(i) + dt * ((v(i) - u(i)) / tau_s);
    
    if v(i+1) >= v_max          % a spike is fired!
        v(i) = v_spike;         % padding the spike amplitude
        v(i+1) = v_r;           % membrane voltage reset      
        u(i+1) = v_sr;          % recovery variable update
    end
end       


figure(1);

% plot input
subplot(4,1,1)
plot(dt*(1:n), I, 'b', 'LineWidth', 4)
title('Input')

% plot time trace
subplot(4,1,2)
plot(dt*(1:n), v, 'k', 'LineWidth', 3); 
title('Time trace')


%% ____________ Plots Phase Plane ________________%
subplot(4,1,3)

Iapp_0 = I(1);

syms vv1 vv2
dot_vv_0 = (g_f*(vv1 - v_f0)^2 - g_s * (vv2 - v_s0)^2 + Iapp_0) / C;
dot_uu = (vv1 - vv2) / tau_s;

hold on

% plot V-nullcline for I =0
nullcline_v_0 = ezplot(dot_vv_0, [-80 80]);
set(nullcline_v_0, 'color',[0, 0.6, 0.748], 'Linewidth', 3);

% plot u-nullcline
nullcline_u = ezplot(dot_uu, [-80 80]);
set(nullcline_u, 'color', [0.4660, 0.6740, 0.1880], 'Linewidth', 3);

% trajectory
plot(v, u, 'color',[1 0.08 0.57],'LineWidth', 2) 
hold off

xlabel('membrane potential V', 'interpreter','latex','fontsize',11);
ylabel('recovery variable u', 'interpreter','latex','fontsize',11);
title('Phase Plane, with I = 0', 'fontsize', 12);


subplot(4,1,4)
Iapp_up = I(500);

hold on
dot_vv_Iapp = (g_f*(vv1 - v_f0)^2 - g_s * (vv2 - v_s0)^2 + Iapp_up) / C;

% plot V-nullcline for I = I_app
nullcline_v_Iapp = ezplot(dot_vv_Iapp, [-80 80]);
set(nullcline_v_Iapp, 'color',[0, 0.6, 0.748], 'Linewidth', 3);

% plot u-nullcline
nullcline_u = ezplot(dot_uu, [-80 80]);
set(nullcline_u, 'color', [0.4660, 0.6740, 0.1880], 'Linewidth', 3);

% trajectory
plot(v, u, 'color',[1 0.08 0.57],'LineWidth', 2) 

xlabel('membrane potential V', 'interpreter','latex','fontsize',11);
ylabel('recovery variable u', 'interpreter','latex','fontsize',11);
title('Phase Plane, with I = I_{app}', 'fontsize', 12);




%% Export/Save
% outfile = ['MQIF_phase_plane', num2str(v0), '_w', num2str(w0), '_e',...
%            num2str(epsilon), '_I', num2str(I_app)];
%        
% suffix_fig = '';
% suffix_data = '';       
% 
% out_mat = [outfile, suffix_data, '.mat'];
% out_fig_png = [outfile, suffix_fig, '.png'];
% out_fig_eps = [outfile, suffix_fig, '.eps'];
% 
% outpath_data = fullfile(dir_base, dir_data, out_mat);
% outpath_fig_png = fullfile(dir_base, dir_fig, out_fig_png);
% outpath_fig_eps = fullfile(dir_base, dir_fig, out_fig_eps);
% 
% % figures
% % saveas(f1, outpath_fig_eps,'eps')
% saveas(f1, outpath_fig_png,'png')

