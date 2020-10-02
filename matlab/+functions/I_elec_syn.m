function I_syn = I_elec_syn(neuron, t_step, g, V)
    
%     n = size(g,1);
%     I_synps =0;
    
%     for k = 1:n
%         I_synps = I_synps + g(j,k) * ( V(tStep,j) - V(tStep,k) );
%     end

    I_syn = ( V(t_step, neuron) - V(t_step, :) ) .* g(neuron, :);
    I_syn = sum(I_syn,2);
end

