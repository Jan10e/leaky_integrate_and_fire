function [v, u] = MQIF_forward_euler(g_f, g_s, v_f0, v_s0, tau_s, C, v_r, v_sr, I, v_max, v_spike, v_init, u_init, Tt)

%simulation of mirrored Fitzhugh Nagumo model
dt = Tt(2)-Tt(1);

%% Euler solve info
v = v_init; 
u = u_init;

%% Integrate using Forward Euler

for i = 1:length(Tt)-1
    
    v(i+1) = v(i) + dt * ((g_f * (v(i) - v_f0)^2 - g_s * (u(i)-v_s0)^2 + I(i)) / C);
    u(i+1) = u(i) + dt * ((v(i) - u(i)) / tau_s);
    
    if v(i+1) >= v_max          % a spike is fired!
        v(i) = v_spike;         % padding the spike amplitude
        v(i+1) = v_r;           % membrane voltage reset      
        u(i+1) = v_sr;          % recovery variable update
    end
    
end   


end