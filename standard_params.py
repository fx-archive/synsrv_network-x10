
from brian2.units import *

N_e = 400
N_i = int(0.2*N_e)

tau = 20*ms                 # membrane time constant
tau_e = 3*ms                # EPSP time constant
tau_i = 5*ms                # IPSP time constant
El = -60*mV                 # resting value
Ee = 0*mV                   # reversal potential Excitation 
Ei = -80*mV                 # reversal potential Inhibition
sigma = 5.0**0.5 *mV        # noise amplitude

Vr_e = -70*mV
Vr_i = -60*mV
Vt_e = -54*mV
Vt_i = -58*mV

ascale = 1.0
a_ee = 1.5 
a_ie = 1.5 
a_ei = -1.5 
a_ii = -1.5 

p_ee = 0.1 
p_ie = 0.1
p_ei = 0.1
p_ii = 0.5

taupre = 15*ms
taupost = 30*ms
Aplus = 15*0.001
Aminus = -7.5*0.001
amax = 40*0.001

synEE_rec = 1

ATotalMax = 40.*0.001

# scaling
scl_active = 0
dt_synEE_scaling = 10*ms

# intrinsic plasticity
it_active = 0
eta_ip = 0.2*mV*ms
it_dt = 10*ms
h_ip = 3*Hz

# structural plasticity
strct_active = 0
prn_thrshld = 0.001 * ATotalMax
insert_P = 0.001
strct_dt = 10*ms
a_insert = 0.01 * ATotalMax


#preT  = 100*second
T  = 10*second
netw_dt = 0.1*ms


# neuron_method = 'euler'
# synEE_method = 'euler'
