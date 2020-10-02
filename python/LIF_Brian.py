from brian2 import *
from numpy import *

if __name__ == '__main__':

    N = 1   # number of neurons

    tau_m   = 10 * ms   # membrane time constant
    v_r     = 0 * mV    # reset potential
    v_th    = 15 * mV   # threshold potential
    I_c     = 20 * mV   # constant input current

    eqs = '''
    dv/dt = -(v-I)/tau_m : volt
    I : volt
    '''

    lif = NeuronGroup(N, model=eqs, threshold=v_th, reset=v_r)

    # You can add randomness in initial membrane potential by changing the following line
    lif.v = v_r * mV + 0 * mV * rand(len(lif))
    lif.I = I_c

    spikes = SpikeMonitor(lif)
    v_trace = StateMonitor(lif, 'v', record=True)

    run(0.1 * second)

    figure(1)
    plot(v_trace.times/ms, v_trace[0]/mV)
    xlabel('Time (ms)', fontsize = 24)
    ylabel('v (mV)', fontsize = 24)
    yticks([0, 4, 8, 12, 16])

    show()

    # plot the f-I curve for a LIF neuron with constant input
    delta_abs = 5
    tau_m = 10
    v_th = 15
    Is = linspace(0, 50, 51)

    f1 = 1.0 / (0 + tau_m * log(Is/(Is - v_th)))
    f2 = 1.0 / (delta_abs + tau_m * log(Is/(Is - v_th)))

    figure(2)
    plot(Is, f1, 'k--')
    plot(Is, f2, 'k--')
    plot([20, 20], [0, 0.3], 'r--')
    xlabel('I', fontsize = 24)
    ylabel('f (kHz)', fontsize = 24)
    yticks([0, .1, .2, .3])

    show()