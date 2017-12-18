
# import logging?

import standard_params as prm
import models as mod

import numpy as np

from brian2.units import ms,mV
from pypet.brian2.parameter import Brian2Parameter, Brian2MonitorResult

from brian2 import NeuronGroup, StateMonitor, SpikeMonitor, run, \
                   PoissonGroup, Synapses, set_device, device, Clock, \
                   defaultclock, prefs, network_operation

def add_params(tr):

    tr.v_standard_parameter=Brian2Parameter
    tr.v_fast_access=True

    tr.f_add_parameter('netw.N_e', prm.N_e)
    tr.f_add_parameter('netw.N_i', prm.N_i)
    
    tr.f_add_parameter('netw.tau',   prm.tau)
    tr.f_add_parameter('netw.tau_e', prm.tau_e)
    tr.f_add_parameter('netw.tau_i', prm.tau_i)
    tr.f_add_parameter('netw.El',    prm.El)
    tr.f_add_parameter('netw.Ee',    prm.Ee)
    tr.f_add_parameter('netw.Ei',    prm.Ei)
    tr.f_add_parameter('netw.sigma', prm.sigma)
    
    tr.f_add_parameter('netw.Vr_e',  prm.Vr_e)
    tr.f_add_parameter('netw.Vr_i',  prm.Vr_i)
    tr.f_add_parameter('netw.Vt_e',  prm.Vt_e)
    tr.f_add_parameter('netw.Vt_i',  prm.Vt_i)
    
    tr.f_add_parameter('netw.a_ee',  prm.a_ee)
    tr.f_add_parameter('netw.a_ie',  prm.a_ie)
    tr.f_add_parameter('netw.a_ei',  prm.a_ei)
    tr.f_add_parameter('netw.a_ii',  prm.a_ii)

    tr.f_add_parameter('netw.p_ee',  prm.p_ee)
    tr.f_add_parameter('netw.p_ie',  prm.p_ie)
    tr.f_add_parameter('netw.p_ei',  prm.p_ei)
    tr.f_add_parameter('netw.p_ii',  prm.p_ii)

    # STDP
    tr.f_add_parameter('netw.taupre',    prm.taupre)
    tr.f_add_parameter('netw.taupost',   prm.taupost)
    tr.f_add_parameter('netw.Aplus',     prm.Aplus)
    tr.f_add_parameter('netw.Aminus',    prm.Aminus)
    tr.f_add_parameter('netw.amax',      prm.amax)

    # scaling
    tr.f_add_parameter('netw.config.scl_active', prm.scl_active)
    tr.f_add_parameter('netw.ATotalMax',        prm.ATotalMax)
    tr.f_add_parameter('netw.dt_synEE_scaling', prm.dt_synEE_scaling)

    # intrinsic plasticity
    tr.f_add_parameter('netw.config.it_active', prm.it_active)
    tr.f_add_parameter('netw.eta_ip', prm.eta_ip)
    tr.f_add_parameter('netw.it_dt',  prm.it_dt)
    tr.f_add_parameter('netw.h_ip',   prm.h_ip)

    # structural plasticity
    tr.f_add_parameter('netw.prn_thrshld', prm.prn_thrshld)
    tr.f_add_parameter('netw.insert_P',    prm.insert_P)
    tr.f_add_parameter('netw.a_insert',    prm.a_insert)
    tr.f_add_parameter('netw.strct_dt',    prm.strct_dt)
    
    tr.f_add_parameter('netw.mod.condlif_sig',   mod.condlif_sig)
    tr.f_add_parameter('netw.mod.nrnEE_thrshld', mod.nrnEE_thrshld)
    tr.f_add_parameter('netw.mod.nrnEE_reset',   mod.nrnEE_reset)
    tr.f_add_parameter('netw.mod.synEE_mod',     mod.synEE_mod)
    tr.f_add_parameter('netw.mod.synEE_pre',     mod.synEE_pre)
    tr.f_add_parameter('netw.mod.synEE_post',    mod.synEE_post)
    tr.f_add_parameter('netw.mod.synEE_scaling', mod.synEE_scaling)
    tr.f_add_parameter('netw.mod.intrinsic_mod', mod.intrinsic_mod)
    tr.f_add_parameter('netw.mod.strct_mod',     mod.strct_mod)
    
    # tr.f_add_parameter('netw.mod.neuron_method', prm.neuron_method)
    # tr.f_add_parameter('netw.mod.synEE_method',  prm.synEE_method)

    #tr.f_add_parameter('netw.sim.preT',  prm.T)
    tr.f_add_parameter('netw.sim.T',  prm.T)
    tr.f_add_parameter('netw.sim.dt', prm.netw_dt)

    tr.f_add_parameter('netw.config.strct_active', prm.strct_active)
   

    
def run_net(tr):

    prefs.codegen.target = 'numpy'
    # prefs.codegen.target = 'cython'
    # set_device('cpp_standalone', directory='./build', build_on_run=False)

    namespace = tr.netw.f_to_dict(short_names=True, fast_access=True)

    defaultclock.dt = tr.netw.sim.dt

    GExc = NeuronGroup(N=tr.N_e, model=tr.condlif_sig, threshold=tr.nrnEE_thrshld,
                       reset=tr.nrnEE_reset, #method=tr.neuron_method,
                       namespace=namespace)
    GInh = NeuronGroup(N=tr.N_i, model=tr.condlif_sig, threshold ='V > Vt',
                       reset='V=Vr_i', #method=tr.neuron_method,
                       namespace=namespace)
    GExc.Vt, GInh.Vt = tr.Vt_e, tr.Vt_i
    GExc.V , GInh.V  = np.random.uniform(tr.Vr_e/mV, tr.Vt_e/mV, size=tr.N_e)*mV, \
                       np.random.uniform(tr.Vr_i/mV, tr.Vt_i/mV, size=tr.N_i)*mV

    SynEE = Synapses(target=GExc, source=GExc, model=tr.synEE_mod,
                     on_pre=tr.synEE_pre, on_post=tr.synEE_post,
                     #method=tr.synEE_method,
                     namespace=namespace)
    SynIE = Synapses(target=GInh, source=GExc, on_pre='ge_post += a_ie',
                     namespace=namespace)
    SynEI = Synapses(target=GExc, source=GInh, on_pre='gi_post += a_ei',
                     namespace=namespace)
    SynII = Synapses(target=GInh, source=GInh, on_pre='gi_post += a_ii',
                     namespace=namespace)

    def generate_connections(N_tar, N_src, p, same=False):
        nums = np.random.binomial(N_tar-1, p, N_src)
        i = np.repeat(np.arange(N_src), nums)
        j = []
        if same:
            for k,n in enumerate(nums):
                j+=list(np.random.choice([*range(k-1)]+[*range(k+1,N_tar)],
                                         size=n, replace=False))
        else:
            for k,n in enumerate(nums):
                j+=list(np.random.choice([*range(N_tar)],
                                         size=n, replace=False))

        return i, np.array(j)

    if not tr.strct_active:
        sEE_src, sEE_tar = generate_connections(tr.N_e, tr.N_e, tr.p_ee, same=True) 

    sIE_src, sIE_tar = generate_connections(tr.N_i, tr.N_e, tr.p_ie)
    sEI_src, sEI_tar = generate_connections(tr.N_e, tr.N_i, tr.p_ei)
    sII_src, sII_tar = generate_connections(tr.N_i, tr.N_i, tr.p_ii, same=True)

    if tr.strct_active:
        SynEE.connect(True)        
    else:
        SynEE.connect(i=sEE_src, j=sEE_tar)
        
    SynIE.connect(i=sIE_src, j=sIE_tar)
    SynEI.connect(i=sEI_src, j=sEI_tar)
    SynII.connect(i=sII_src, j=sII_tar)

    if not tr.strct_active:
        tr.f_add_result('sEE_src', sEE_src)
        tr.f_add_result('sEE_tar', sEE_tar)
        
    tr.f_add_result('sIE_src', sIE_src)
    tr.f_add_result('sIE_tar', sIE_tar)
    tr.f_add_result('sEI_src', sEI_src)
    tr.f_add_result('sEI_tar', sEI_tar)
    tr.f_add_result('sII_src', sII_src)
    tr.f_add_result('sII_tar', sII_tar)

    SynEE.a = tr.a_ee
    SynEE.syn_active = 0
    SynEE.insert_P = tr.insert_P

    # synaptic scaling
    if tr.netw.config.scl_active:
        SynEE.summed_updaters['Asum_post']._clock = Clock(
            dt=tr.dt_synEE_scaling)
        SynEE.run_regularly(tr.synEE_scaling, dt = tr.dt_synEE_scaling,
                            when='end')

    # intrinsic plasticity
    if tr.netw.config.it_active:
        GExc.h_ip = tr.h_ip
        GExc.run_regularly(tr.intrinsic_mod, dt = tr.it_dt, when='end')

    # structural plasticity
    if tr.netw.config.strct_active:
        print("activated strct_mod")
        SynEE.run_regularly(tr.strct_mod, dt = tr.strct_dt, when='end')

        @network_operation(dt=10*ms, when='end')
        def f():
            print("Hello World")

    #run(tr.sim.preT)
    
    GExc_stat = StateMonitor(GExc, ['V', 'Vt', 'ge', 'gi'], record=[0,1,2])
    SynEE_stat = StateMonitor(SynEE, ['a','Apre', 'Apost'], record=[0,1,2])

    GExc_spks = SpikeMonitor(GExc)
    
    GInh_stat = StateMonitor(GInh, ['V', 'Vt', 'ge', 'gi'], record=[0,1,2])
    GInh_spks = SpikeMonitor(GInh)

    GExc_vts = StateMonitor(GExc, ['Vt'], record=True, dt=tr.sim.T/2.)
    SynEE_a = StateMonitor(SynEE, ['a','syn_active'], record=True, dt=tr.sim.T/2.)

    run(tr.sim.T)
    #device.build(directory='./build')

    GExc_vts.record_single_timestep()
    SynEE_a.record_single_timestep()

    tr.v_standard_result = Brian2MonitorResult
    tr.f_add_result('GExc_stat', GExc_stat)
    tr.f_add_result('SynEE_stat', SynEE_stat)
    print("Saving exc spikes...   ", GExc_spks.get_states()['N'])
    tr.f_add_result('GExc_spks', GExc_spks)
    tr.f_add_result('GInh_stat', GInh_stat)
    print("Saving inh spikes...   ", GInh_spks.get_states()['N'])
    tr.f_add_result('GInh_spks', GInh_spks)
    tr.f_add_result('SynEE_a', SynEE_a)
    tr.f_add_result('GExc_vts', GExc_vts)

