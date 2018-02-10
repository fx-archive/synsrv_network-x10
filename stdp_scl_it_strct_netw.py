
import os

import standard_params as prm
import models as mod
from utils import generate_connections, generate_full_connectivity

import numpy as np

from brian2.units import ms,mV,second
from pypet.brian2.parameter import Brian2Parameter, Brian2MonitorResult

from brian2 import NeuronGroup, StateMonitor, SpikeMonitor, run, \
                   PoissonGroup, Synapses, set_device, device, Clock, \
                   defaultclock, prefs, network_operation, implementation, \
                   check_units


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

    tr.f_add_parameter('netw.ascale', prm.ascale)
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
    tr.f_add_parameter('netw.synEE_rec',      prm.synEE_rec)
    

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
    # tr.f_add_parameter('netw.mod.synEE_mod',     mod.synEE_mod)
    # gets added later taking into accout tr.synEE_rec!
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


# @implementation('cpp', code=r'''
# #include <fstream>

# double monitor_now(double t, int i, int j, double w) {
#     std::ofstream outfile;

#     outfile.open("../spk_vals", std::ios_base::app);
#     outfile << t << "," << i << "," << j << "," << w << "\n";

#     return 0.0; // we need to return a dummy value
# } '''
#                 )
# @check_units(t=second, i=1, j=1, w=1, result=1)
# def record_at_spike_time(t, i, j, w):
#     return 0.0

# need to use vATotalMax here instead of ATotalMax, as otherwise
# the float of ATotalMax gets inserted into the function
# parameters:
#
#   double syn_scale(double a, double 0.121, ...)
#
# which returns an error
@implementation('cpp', code=r'''
   
    double syn_scale(double a, double vATotalMax, double Asum_post) {
      
      double a_out;

      if (Asum_post==0.){
          a_out = 0.;
      }
      else{
          a_out = a*(vATotalMax/Asum_post);
      }

      return a_out;
    } ''')
@check_units(a=1, vATotalMax=1, Asum_post=1, result=1)
def syn_scale(a, vATotalMax, Asum_post):
    return -1.


@implementation('cpp', code=r'''
    #include <fstream>
    
    double record_turnover(double t, int was_active_before, int should_become_active, int should_stay_active, int syn_active, int i, int j) {

      if (int(was_active_before==0)*should_become_active==1){
          std::ofstream outfile;          
          outfile.open("turnover", std::ios_base::app);
          outfile << 1 << "," << t << "," << i << "," << j << "\n";
      }
      else if (was_active_before*int(should_stay_active==0)){
           std::ofstream outfile;     
           outfile.open("turnover", std::ios_base::app);
           outfile << 0 << "," << t << "," << i << "," << j << "\n";
      }

      return 0.0; // we need to return a dummy value
    } ''')

@check_units(t=second, was_active_before=1, should_become_active=1,
             should_stay_active=1, syn_active=1, i=1, j=1, result=1)
def record_turnover(t, was_active_before, should_become_active,
                    should_stay_active, syn_active, i, j):
    return 0.0

@implementation('cpp', code=r'''
    #include <fstream>
    
    double record_spk(double t, int i, int j, double a, double Apre, double Apost, int syn_active, int preorpost) {

       if (syn_active > 0){
          std::ofstream outfile;          
          outfile.open("spk_register", std::ios_base::app);
          outfile << t << "," << i << "," << j << "," << a << "," << Apre << "," << Apost << "," << preorpost << "\n";
       }

      return 0.0; // we need to return a dummy value
    } ''')

@check_units(t=second, i=1, j=1, a=1, Apre=1, Apost=1, syn_active=1,
             preorpost=1, result=1)
def record_spk(t, i, j, a, Apre, Apost, syn_active, preorpost):
    return 0.0

    
def run_net(tr):

    # prefs.codegen.target = 'numpy'
    # prefs.codegen.target = 'cython'
    set_device('cpp_standalone', directory='./builds/%.4d'%(tr.v_idx),
               build_on_run=False)

    print("Started process with id ", str(tr.v_idx))

    namespace = tr.netw.f_to_dict(short_names=True, fast_access=True)
    namespace['idx'] = tr.v_idx

    defaultclock.dt = tr.netw.sim.dt

    GExc = NeuronGroup(N=tr.N_e, model=tr.condlif_sig,
                       threshold=tr.nrnEE_thrshld,
                       reset=tr.nrnEE_reset, #method=tr.neuron_method,
                       namespace=namespace)
    GInh = NeuronGroup(N=tr.N_i, model=tr.condlif_sig,
                       threshold ='V > Vt',
                       reset='V=Vr_i', #method=tr.neuron_method,
                       namespace=namespace)

    # set initial thresholds fixed, init. potentials uniformly distrib.
    GExc.Vt, GInh.Vt = tr.Vt_e, tr.Vt_i
    GExc.V , GInh.V  = np.random.uniform(tr.Vr_e/mV, tr.Vt_e/mV,
                                         size=tr.N_e)*mV, \
                       np.random.uniform(tr.Vr_i/mV, tr.Vt_i/mV,
                                         size=tr.N_i)*mV

    if tr.synEE_rec:
        tr.f_add_parameter('synEE_pre', mod.syn_EE_pre + mod.synEE_pre_rec)
        tr.f_add_parameter('synEE_post', mod.syn_EE_post + mod.synEE_post_rec)
    
    # E<-E advanced synapse model, rest simple
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

    

    if tr.strct_active:
        sEE_src, sEE_tar = generate_full_connectivity(tr.N_e, same=True)
        SynEE.connect(i=sEE_src, j=sEE_tar)
        SynEE.syn_active = 0
    else:
        sEE_src, sEE_tar = generate_connections(tr.N_e, tr.N_e, tr.p_ee,
                                                same=True)
        SynEE.connect(i=sEE_src, j=sEE_tar)

        tr.f_add_result('sEE_src', sEE_src)
        tr.f_add_result('sEE_tar', sEE_tar)

        SynEE.syn_active = 1


    sIE_src, sIE_tar = generate_connections(tr.N_i, tr.N_e, tr.p_ie)
    sEI_src, sEI_tar = generate_connections(tr.N_e, tr.N_i, tr.p_ei)
    sII_src, sII_tar = generate_connections(tr.N_i, tr.N_i, tr.p_ii,
                                            same=True)

    SynIE.connect(i=sIE_src, j=sIE_tar)
    SynEI.connect(i=sEI_src, j=sEI_tar)
    SynII.connect(i=sII_src, j=sII_tar)
        
    tr.f_add_result('sIE_src', sIE_src)
    tr.f_add_result('sIE_tar', sIE_tar)
    tr.f_add_result('sEI_src', sEI_src)
    tr.f_add_result('sEI_tar', sEI_tar)
    tr.f_add_result('sII_src', sII_src)
    tr.f_add_result('sII_tar', sII_tar)

    if tr.strct_active:
        SynEE.a = 0
    else:
        SynEE.a = tr.a_ee
        
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
    number_active_synapses = []

    if tr.netw.config.strct_active:
        SynEE.run_regularly(tr.strct_mod, dt = tr.strct_dt, when='end')

    # -------------- recording ------------------        

    #run(tr.sim.preT)

    
    # GExc_stat = StateMonitor(GExc, ['V', 'Vt', 'ge', 'gi'], record=[0,1,2])
    # SynEE_stat = StateMonitor(SynEE, ['a','Apre', 'Apost'], record=range(20))

    # GExc_spks = SpikeMonitor(GExc)
    
    # GInh_stat = StateMonitor(GInh, ['V', 'Vt', 'ge', 'gi'], record=[0,1,2])
    # GInh_spks = SpikeMonitor(GInh)

    # GExc_vts = StateMonitor(GExc, ['Vt'], record=True, dt=tr.sim.T/2.)
    SynEE_a = StateMonitor(SynEE, ['a','syn_active'],
                           record=range(tr.N_e*(tr.N_e-1)), dt=tr.sim.T/10.)
    

    GExc_stat = []
    SynEE_stat = []

    GExc_spks = []
    
    GInh_stat = []
    GInh_spks = []

    GExc_vts = []
    SynEE_a = []
    
    run(tr.sim.T, report='text')
    device.build(directory='./builds/%.4d'%(tr.v_idx))

    #GExc_vts.record_single_timestep()
    SynEE_a.record_single_timestep()

    # it looks like only pure numpy arrays can be stored as results
    number_active_synapses = np.array(number_active_synapses)
    tr.f_add_result('SynAct_stat', number_active_synapses)

    # it looks like only pure numpy arrays can be stored as results
    # tr.f_add_result('dead_times', np.array(dead_times))
    # tr.f_add_result('life_times', np.array(life_times))
 
    
    tr.v_standard_result = Brian2MonitorResult

    tr.f_add_result('GExc_stat', GExc_stat)
    tr.f_add_result('SynEE_stat', SynEE_stat)
    #print("Saving exc spikes...   ", GExc_spks.get_states()['N'])
    tr.f_add_result('GExc_spks', GExc_spks)
    tr.f_add_result('GInh_stat', GInh_stat)
    #print("Saving inh spikes...   ", GInh_spks.get_states()['N'])
    tr.f_add_result('GInh_spks', GInh_spks)
    tr.f_add_result('SynEE_a', SynEE_a)

    tr.f_add_result('GExc_vts', GExc_vts)

    # ----------------- add raw data ------------------------
    fpath = './builds/%.4d/'%(tr.v_idx)

    from pathlib import Path

    Path(fpath+'turnover').touch()
    turnover_data = np.genfromtxt(fpath+'turnover',delimiter=',')
    tr.f_add_result('turnover', turnover_data)
    os.remove(fpath+'turnover')

    Path(fpath+'spk_register').touch()
    spk_register_data = np.genfromtxt(fpath+'spk_register',delimiter=',')
    tr.f_add_result('spk_register', spk_register_data)
    os.remove(fpath+'spk_register')
    
