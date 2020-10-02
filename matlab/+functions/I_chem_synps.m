function [I_syn,s] = I_chem_synps(iN,iT,g,S,N,alpha,beta,dt,E_syn,v_a1,E_syn_inh)

n = size(g,1);
s = zeros(1,n);

I_syn = g(iN,:) .* S(iN,:,iT) .* (v_a1 - E_syn(iN,:));
I_syn = sum(I_syn,2);

% g_syn = g(j,:) .* S(j,:,tStep);
% g_syn = sum(g_syn,2);
% i_synps = g_syn * ( v_a1 - E_syn(j) );


for i = 1:n
    
    if E_syn(iN,i) == E_syn_inh
        alpha_1 = alpha(2);
        beta_1 = beta(2);
    else
        alpha_1 = alpha(1);
        beta_1 = beta(2);
    end
    
    s(i) = S(iN,i,iT) + dt * ( alpha_1 * N(iT,i) *(1-S(iN,i,iT)) - ...
        beta_1 * S(iN,i,iT));
end

end

