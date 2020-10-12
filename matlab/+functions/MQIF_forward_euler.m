function [V, V_s] = MQIF_forward_euler(g_f, g_s, V_f0, V_s0, tau_s, C, V_init, Vs_init, I, V_th, Tt)
%simulation of mirrored Fitzhugh Nagumo model

dt = Tt(2)-Tt(1);

%% Resets
V_rest = -70;        % resting potential fast gating
V_srest = -30;       % resting potential slow gating

%% Euler solve info
V = V_init;
V_s = Vs_init;

%% Integrate using Forward Euler
for i = 1:length(Tt)-1
    V(i + 1)    = V(i) + dt .* ((g_f*(V(i) - V_f0)^2 - g_s*(V_s(i) - V_s0)^2 + I) / C);
    V_s(i + 1)  = V_s(i) + dt .* ((V(i) - V_s0) / tau_s);
    
    %reset
    if V >= V_th
        V = V_rest;
        V_s = V_srest;
    end
    
    
    
end

end