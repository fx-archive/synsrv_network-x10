
import sys, os, shutil, pickle, neo, scipy

from . import models as mod
from .utils import generate_connections, generate_full_connectivity, \
                   generate_N_connections

import numpy as np

from brian2.units import ms,mV,second,Hz
from pypet.brian2.parameter import Brian2MonitorResult

from brian2 import NeuronGroup, StateMonitor, SpikeMonitor, run, \
                   PoissonGroup, Synapses, set_device, device, Clock, \
                   defaultclock, prefs, network_operation, Network, \
                   PoissonGroup, PopulationRateMonitor, profiling_summary

from elephant.conversion import BinnedSpikeTrain
from elephant.spike_train_correlation import corrcoef, cch
import quantities as pq


from .cpp_methods import syn_scale, syn_EI_scale, \
                         record_turnover, record_turnover_EI, \
                         record_spk, record_spk_EI


    
def run_net(tr):

    # prefs.codegen.target = 'numpy'
    # prefs.codegen.target = 'cython'
    if tr.n_threads > 1:
        prefs.devices.cpp_standalone.openmp_threads = tr.n_threads
        
    set_device('cpp_standalone', directory='./builds/%.4d'%(tr.v_idx),
               build_on_run=False)

    print("Started process with id ", str(tr.v_idx))

    T = tr.T1 + tr.T2 + tr.T3 + tr.T4 + tr.T5

    namespace = tr.netw.f_to_dict(short_names=True, fast_access=True)
    namespace['idx'] = tr.v_idx

    defaultclock.dt = tr.netw.sim.dt

    # collect all network components dependent on configuration
    # (e.g. poisson vs. memnoise) and add them to the Brian 2
    # network object later
    netw_objects = []

    if tr.external_mode=='memnoise':
        neuron_model = tr.condlif_memnoise
    elif tr.external_mode=='poisson':
        neuron_model = tr.condlif_poisson

    GExc = NeuronGroup(N=tr.N_e, model=neuron_model,
                       threshold=tr.nrnEE_thrshld,
                       reset=tr.nrnEE_reset, #method=tr.neuron_method,
                       namespace=namespace)
    GInh = NeuronGroup(N=tr.N_i, model=neuron_model,
                       threshold ='V > Vt',
                       reset='V=Vr_i', #method=tr.neuron_method,
                       namespace=namespace)

    if tr.external_mode=='memnoise':
        GExc.mu, GInh.mu = tr.mu_e, tr.mu_i
        GExc.sigma, GInh.sigma = tr.sigma_e, tr.sigma_i

  
    GExc.Vt, GInh.Vt = tr.Vt_e, tr.Vt_i
    GExc.V , GInh.V  = np.random.uniform(tr.Vr_e/mV, tr.Vt_e/mV,
                                         size=tr.N_e)*mV, \
                       np.random.uniform(tr.Vr_i/mV, tr.Vt_i/mV,
                                         size=tr.N_i)*mV

    netw_objects.extend([GExc,GInh])


    if tr.external_mode=='poisson':
    
        if tr.PInp_mode == 'pool':
            PInp = PoissonGroup(tr.NPInp, rates=tr.PInp_rate,
                                namespace=namespace, name='poissongroup_exc')
            sPN = Synapses(target=GExc, source=PInp, model=tr.poisson_mod,
                           on_pre='gfwd_post += a_EPoi',
                           namespace=namespace, name='synPInpExc')

            sPN_src, sPN_tar = generate_N_connections(N_tar=tr.N_e,
                                                      N_src=tr.NPInp,
                                                      N=tr.NPInp_1n)

        elif tr.PInp_mode == 'indep':
            PInp = PoissonGroup(tr.N_e, rates=tr.PInp_rate,
                                namespace=namespace)
            sPN = Synapses(target=GExc, source=PInp, model=tr.poisson_mod,
                           on_pre='gfwd_post += a_EPoi',
                           namespace=namespace, name='synPInp_inhInh')
            sPN_src, sPN_tar = range(tr.N_e), range(tr.N_e)


        sPN.connect(i=sPN_src, j=sPN_tar)



        if tr.PInp_mode == 'pool':
            PInp_inh = PoissonGroup(tr.NPInp_inh, rates=tr.PInp_inh_rate,
                                    namespace=namespace, name='poissongroup_inh')
            sPNInh = Synapses(target=GInh, source=PInp_inh, model=tr.poisson_mod,
                               on_pre='gfwd_post += a_EPoi',
                               namespace=namespace)
            sPNInh_src, sPNInh_tar = generate_N_connections(N_tar=tr.N_i,
                                                            N_src=tr.NPInp_inh,
                                                            N=tr.NPInp_inh_1n)


        elif tr.PInp_mode == 'indep':

            PInp_inh = PoissonGroup(tr.N_i, rates=tr.PInp_inh_rate,
                                    namespace=namespace)
            sPNInh = Synapses(target=GInh, source=PInp_inh, model=tr.poisson_mod,
                              on_pre='gfwd_post += a_EPoi',
                              namespace=namespace)
            sPNInh_src, sPNInh_tar = range(tr.N_i), range(tr.N_i)


        sPNInh.connect(i=sPNInh_src, j=sPNInh_tar)

        netw_objects.extend([PInp, sPN, PInp_inh, sPNInh])
    

        
    if tr.syn_noise:
        synEE_mod = '''%s 
                       %s
                       %s''' %(tr.synEE_noise, tr.synEE_mod, tr.synEE_scl_mod)

        synEI_mod = '''%s 
                       %s
                       %s''' %(tr.synEE_noise, tr.synEE_mod, tr.synEI_scl_mod)

    else:
        synEE_mod = '''%s 
                       %s
                       %s''' %(tr.synEE_static, tr.synEE_mod, tr.synEE_scl_mod)

        synEI_mod = '''%s 
                       %s
                       %s''' %(tr.synEE_static, tr.synEE_mod, tr.synEI_scl_mod)


    synEE_pre_mod = mod.synEE_pre
    synEE_post_mod = mod.syn_post
    
    if tr.stdp_active:
        synEE_pre_mod  = '''%s 
                            %s''' %(synEE_pre_mod, mod.syn_pre_STDP)
        synEE_post_mod = '''%s 
                            %s''' %(synEE_post_mod, mod.syn_post_STDP)

    if tr.synEE_rec:
        synEE_pre_mod  = '''%s 
                            %s''' %(synEE_pre_mod, mod.synEE_pre_rec)
        synEE_post_mod = '''%s 
                            %s''' %(synEE_post_mod, mod.synEE_post_rec)


        
    # E<-E advanced synapse model
    SynEE = Synapses(target=GExc, source=GExc, model=synEE_mod,
                     on_pre=synEE_pre_mod, on_post=synEE_post_mod,
                     namespace=namespace, dt=tr.synEE_mod_dt)

    if tr.istdp_active and tr.istdp_type=='dbexp':

        synEI_pre_mod  = '''%s 
                            %s''' %(mod.synEI_pre, mod.syn_pre_STDP)
        synEI_post_mod = '''%s 
                            %s''' %(mod.syn_post, mod.syn_post_STDP)

    elif tr.istdp_active and tr.istdp_type=='sym':

        synEI_pre_mod  = '''%s 
                            %s''' %(mod.synEI_pre_sym, mod.syn_pre_STDP)
        synEI_post_mod = '''%s 
                            %s''' %(mod.synEI_post_sym, mod.syn_post_STDP)

    if tr.istdp_active and tr.synEI_rec:

            synEI_pre_mod  = '''%s 
                                %s''' %(synEI_pre_mod, mod.synEI_pre_rec)
            synEI_post_mod = '''%s 
                                %s''' %(synEI_post_mod, mod.synEI_post_rec)
            
    if tr.istdp_active:        
        SynEI = Synapses(target=GExc, source=GInh, model=synEI_mod,
                         on_pre=synEI_pre_mod, on_post=synEI_post_mod,
                         namespace=namespace, dt=tr.synEE_mod_dt)
    else:
        SynEI = Synapses(target=GExc, source=GInh, on_pre='gi_post += a_ei',
                         namespace=namespace)

    #other simple  
    SynIE = Synapses(target=GInh, source=GExc, on_pre='ge_post += a_ie',
                     namespace=namespace)

    SynII = Synapses(target=GInh, source=GInh, on_pre='gi_post += a_ii',
                     namespace=namespace)

    if tr.strct_active:
        sEE_src, sEE_tar = generate_full_connectivity(tr.N_e, same=True)
        SynEE.connect(i=sEE_src, j=sEE_tar)
        SynEE.syn_active = 0

    else:
        sEE_src, sEE_tar = generate_connections(tr.N_e, tr.N_e, tr.p_ee, same=True)
        SynEE.connect(i=sEE_src, j=sEE_tar)
        SynEE.syn_active = 1

    if tr.istrct_active:
        print('istrct active')
        sEI_src, sEI_tar = generate_full_connectivity(Nsrc=tr.N_i, Ntar=tr.N_e, same=False)
        SynEI.connect(i=sEI_src, j=sEI_tar)
        SynEI.syn_active = 0

    else:
        print('istrct not active')
        sEI_src, sEI_tar = generate_connections(tr.N_e, tr.N_i, tr.p_ei)
        print(len(sEI_src))
        SynEI.connect(i=sEI_src, j=sEI_tar)

        # initial values, as they are not later set
        # by istrct initialization
        SynEI.a = tr.a_ei
        SynEI.syn_active = 1

        

    sIE_src, sIE_tar = generate_connections(tr.N_i, tr.N_e, tr.p_ie)
    sII_src, sII_tar = generate_connections(tr.N_i, tr.N_i, tr.p_ii,
                                            same=True)

    SynIE.connect(i=sIE_src, j=sIE_tar)
    SynII.connect(i=sII_src, j=sII_tar)

    tr.f_add_result('sEE_src', sEE_src)
    tr.f_add_result('sEE_tar', sEE_tar)
    tr.f_add_result('sIE_src', sIE_src)
    tr.f_add_result('sIE_tar', sIE_tar)
    tr.f_add_result('sEI_src', sEI_src)
    tr.f_add_result('sEI_tar', sEI_tar)
    tr.f_add_result('sII_src', sII_src)
    tr.f_add_result('sII_tar', sII_tar)


    if tr.syn_noise:
        SynEE.syn_sigma = tr.syn_sigma

    if tr.syn_noise and tr.istdp_active:
        SynEI.syn_sigma = tr.syn_sigma

    SynEE.insert_P = tr.insert_P
    SynEE.p_inactivate = tr.p_inactivate
    SynEE.stdp_active=1

    if tr.istdp_active:
        SynEI.insert_P = tr.insert_P_ei
        SynEI.p_inactivate = tr.p_inactivate_ei
        SynEI.stdp_active=1
     

    # make randomly chosen synapses active at beginning
    rs = np.random.uniform(size=tr.N_e*(tr.N_e-1))
    initial_active = (rs < tr.p_ee).astype('int')
    initial_a = initial_active * tr.a_ee
    SynEE.syn_active = initial_active
    SynEE.a = initial_a

    if tr.istrct_active:
        rs = np.random.uniform(size=tr.N_i*tr.N_e)
        initial_active = (rs < tr.p_ei).astype('int')
        initial_a = initial_active * tr.a_ei
        SynEI.syn_active = initial_active
        SynEI.a = initial_a



    # recording of stdp in T4
    SynEE.stdp_rec_start = tr.T1+tr.T2+tr.T3
    SynEE.stdp_rec_max = tr.T1+tr.T2+tr.T3 + tr.stdp_rec_T

    if tr.istdp_active:
        SynEI.stdp_rec_start = tr.T1+tr.T2+tr.T3
        SynEI.stdp_rec_max = tr.T1+tr.T2+tr.T3 + tr.stdp_rec_T

  
       
    # synaptic scaling
    if tr.netw.config.scl_active:

        if tr.syn_scl_rec:
            SynEE.scl_rec_start = tr.T1+tr.T2+tr.T3
            SynEE.scl_rec_max = tr.T1+tr.T2+tr.T3 + tr.scl_rec_T
        else:
            SynEE.scl_rec_start = T+10*second
            SynEE.scl_rec_max = T
        
        SynEE.summed_updaters['AsumEE_post']._clock = Clock(
            dt=tr.dt_synEE_scaling)
        synee_scaling = SynEE.run_regularly(tr.synEE_scaling,
                                            dt=tr.dt_synEE_scaling,
                                            when='end',
                                            name='synEE_scaling')

    if tr.netw.config.iscl_active:

        if tr.syn_iscl_rec:
            SynEI.scl_rec_start = tr.T1+tr.T2+tr.T3
            SynEI.scl_rec_max = tr.T1+tr.T2+tr.T3 + tr.scl_rec_T
        else:
            SynEI.scl_rec_start = T+10*second
            SynEI.scl_rec_max = T
        
        SynEI.summed_updaters['AsumEI_post']._clock = Clock(
            dt=tr.dt_synEE_scaling)

        synei_scaling = SynEI.run_regularly(tr.synEI_scaling,
                                            dt=tr.dt_synEE_scaling,
                                            when='end',
                                            name='synEI_scaling')


    # # intrinsic plasticity
    # if tr.netw.config.it_active:
    #     GExc.h_ip = tr.h_ip
    #     GExc.run_regularly(tr.intrinsic_mod, dt = tr.it_dt, when='end')

    # structural plasticity
    if tr.netw.config.strct_active:
        if tr.strct_mode == 'zero':    
            if tr.turnover_rec:
                strct_mod  = '''%s 
                                %s''' %(tr.strct_mod, tr.turnover_rec_mod)
            else:
                strct_mod = tr.strct_mod
                
            strctplst = SynEE.run_regularly(strct_mod, dt=tr.strct_dt,
                                            when='end', name='strct_plst_zero')
           
        elif tr.strct_mode == 'thrs':
            if tr.turnover_rec:
                strct_mod_thrs  = '''%s 
                                %s''' %(tr.strct_mod_thrs, tr.turnover_rec_mod)
            else:
                strct_mod_thrs = tr.strct_mod_thrs
                
            strctplst = SynEE.run_regularly(strct_mod_thrs,
                                            dt=tr.strct_dt,
                                            when='end',
                                            name='strct_plst_thrs')



    if tr.netw.config.istrct_active:
        if tr.strct_mode == 'zero':    
            if tr.turnover_rec:
                strct_mod_EI  = '''%s 
                                   %s''' %(tr.strct_mod, tr.turnoverEI_rec_mod)
            else:
                strct_mod_EI = tr.strct_mod
                
            strctplst_EI = SynEI.run_regularly(strct_mod_EI, dt=tr.strct_dt,
                                               when='end', name='strct_plst_EI')
           
        elif tr.strct_mode == 'thrs':
            raise NotImplementedError

    netw_objects.extend([SynEE, SynEI, SynIE, SynII])


    # keep track of the number of active synapses
    sum_target = NeuronGroup(1, 'c : 1 (shared)', dt=tr.csample_dt)

    sum_model = '''NSyn : 1 (constant)
                   c_post = (1.0*syn_active_pre)/NSyn : 1 (summed)'''
    sum_connection = Synapses(target=sum_target, source=SynEE,
                              model=sum_model, dt=tr.csample_dt,
                              name='get_active_synapse_count')
    sum_connection.connect()
    sum_connection.NSyn = tr.N_e * (tr.N_e-1)
    

    if tr.adjust_insertP:
        # homeostatically adjust growth rate
        growth_updater = Synapses(sum_target, SynEE)
        growth_updater.run_regularly('insert_P_post *= 0.1/c_pre',
                                     when='after_groups', dt=tr.csample_dt,
                                     name='update_insP')
        growth_updater.connect(j='0')

        netw_objects.extend([sum_target, sum_connection, growth_updater])


    if tr.istrct_active:

        # keep track of the number of active synapses
        sum_target_EI = NeuronGroup(1, 'c : 1 (shared)', dt=tr.csample_dt)

        sum_model_EI = '''NSyn : 1 (constant)
                          c_post = (1.0*syn_active_pre)/NSyn : 1 (summed)'''
        sum_connection_EI = Synapses(target=sum_target_EI, source=SynEI,
                                     model=sum_model_EI, dt=tr.csample_dt,
                                     name='get_active_synapse_count_EI')
        sum_connection_EI.connect()
        sum_connection_EI.NSyn = tr.N_e * tr.N_i



        if tr.adjust_EI_insertP:
            # homeostatically adjust growth rate
            growth_updater_EI = Synapses(sum_target_EI, SynEI)
            growth_updater_EI.run_regularly('insert_P_post *= 0.1/c_pre',
                                            when='after_groups', dt=tr.csample_dt,
                                            name='update_insP_EI')
            growth_updater_EI.connect(j='0')

            netw_objects.extend([sum_target_EI, sum_connection_EI, growth_updater_EI])



            
    # -------------- recording ------------------        

    GExc_recvars = []
    if tr.memtraces_rec:
        GExc_recvars.append('V')
    if tr.vttraces_rec:
        GExc_recvars.append('Vt')
    if tr.getraces_rec:
        GExc_recvars.append('ge')
    if tr.gitraces_rec:
        GExc_recvars.append('gi')
    if tr.gfwdtraces_rec and tr.external_mode=='poisson':
        GExc_recvars.append('gfwd')

    GInh_recvars = GExc_recvars
    
    GExc_stat = StateMonitor(GExc, GExc_recvars, record=[0,1,2],
                             dt=tr.GExc_stat_dt)
    GInh_stat = StateMonitor(GInh, GInh_recvars, record=[0,1,2],
                             dt=tr.GInh_stat_dt)
    
    # SynEE stat
    SynEE_recvars = []
    if tr.synee_atraces_rec:
        SynEE_recvars.append('a')
    if tr.synee_activetraces_rec:
        SynEE_recvars.append('syn_active')
    if tr.synee_Apretraces_rec:
        SynEE_recvars.append('Apre')
    if tr.synee_Aposttraces_rec:
        SynEE_recvars.append('Apost')

    SynEE_stat = StateMonitor(SynEE, SynEE_recvars,
                              record=range(tr.n_synee_traces_rec),
                              when='end', dt=tr.synEE_stat_dt)

    if tr.istdp_active:
        # SynEI stat
        SynEI_recvars = []
        if tr.synei_atraces_rec:
            SynEI_recvars.append('a')
        if tr.synei_activetraces_rec:
            SynEI_recvars.append('syn_active')
        if tr.synei_Apretraces_rec:
            SynEI_recvars.append('Apre')
        if tr.synei_Aposttraces_rec:
            SynEI_recvars.append('Apost')

        SynEI_stat = StateMonitor(SynEI, SynEI_recvars,
                                  record=range(tr.n_synei_traces_rec),
                                  when='end', dt=tr.synEI_stat_dt)
        netw_objects.append(SynEI_stat)
        

    if tr.adjust_insertP:

        C_stat = StateMonitor(sum_target, 'c', dt=tr.csample_dt,
                              record=[0], when='end')
        insP_stat = StateMonitor(SynEE, 'insert_P', dt=tr.csample_dt,
                                 record=[0], when='end')
        netw_objects.extend([C_stat, insP_stat])

    if tr.adjust_EI_insertP:

        C_EI_stat = StateMonitor(sum_target_EI, 'c', dt=tr.csample_dt,
                                 record=[0], when='end')
        insP_EI_stat = StateMonitor(SynEI, 'insert_P', dt=tr.csample_dt,
                                    record=[0], when='end')
        netw_objects.extend([C_EI_stat, insP_EI_stat])


    
    GExc_spks = SpikeMonitor(GExc)    
    GInh_spks = SpikeMonitor(GInh)

    GExc_rate = PopulationRateMonitor(GExc)
    GInh_rate = PopulationRateMonitor(GInh)

    if tr.external_mode=='poisson':
        PInp_spks = SpikeMonitor(PInp)
        PInp_rate = PopulationRateMonitor(PInp)
        netw_objects.extend([PInp_spks,PInp_rate])


    if tr.synee_a_nrecpoints==0:
        SynEE_a_dt = 10*tr.sim.T2
    else:
        SynEE_a_dt = tr.sim.T2/tr.synee_a_nrecpoints
    SynEE_a = StateMonitor(SynEE, ['a','syn_active'],
                           record=range(tr.N_e*(tr.N_e-1)),
                           dt=SynEE_a_dt,
                           when='end', order=100)

    if tr.istdp_active and tr.synei_a_nrecpoints>0:
        SynEI_a_dt = tr.sim.T2/tr.synei_a_nrecpoints

        if tr.istrct_active:
            record_range = range(tr.N_e*tr.N_i)
        else:
            record_range = range(len(sEI_src))

        SynEI_a = StateMonitor(SynEI, ['a','syn_active'],
                               record=record_range,
                               dt=SynEI_a_dt,
                               when='end', order=100)

        netw_objects.append(SynEI_a)

        

    netw_objects.extend([GExc_stat, GInh_stat,
                         SynEE_stat, SynEE_a, 
                         GExc_spks, GInh_spks,
                         GExc_rate, GInh_rate])
    

    net = Network(*netw_objects)

    
    def set_active(*argv):
        for net_object in argv:
            net_object.active=True

    def set_inactive(*argv):
        for net_object in argv:
            net_object.active=False



    ### Simulation periods
            

    # --------- T1 ---------
    # initial recording period,
    # all recorders active

    T1T3_recorders = [GExc_spks, GInh_spks, 
                      SynEE_stat, 
                      GExc_stat, GInh_stat,
                      GExc_rate, GInh_rate]

    if tr.istdp_active:
        T1T3_recorders.append(SynEI_stat)
    

    set_active(*T1T3_recorders)

    if tr.external_mode=='poisson':
        set_active(PInp_spks, PInp_rate)
       
    net.run(tr.sim.T1, report='text',
            report_period=300*second, profile=True)

    # --------- T2 ---------
    # main simulation period
    # only active recordings are:
    #   1) turnover 2) C_stat 3) SynEE_a
    
    set_inactive(*T1T3_recorders)
    
    if tr.external_mode=='poisson':
        set_inactive(PInp_spks, PInp_rate)

    net.run(tr.sim.T2, report='text',
            report_period=300*second, profile=True)

    # --------- T3 ---------
    # second recording period,
    # all recorders active

    set_active(*T1T3_recorders)
    
    if tr.external_mode=='poisson':
        set_active(PInp_spks, PInp_rate)
    
    net.run(tr.sim.T3, report='text',
            report_period=300*second, profile=True)

    # --------- T4 ---------
    # record STDP and scaling weight changes to file
    # through the cpp models
    
    set_inactive(*T1T3_recorders)

    if tr.external_mode=='poisson':
        set_inactive(PInp_spks, PInp_rate)
    
    net.run(tr.sim.T4, report='text',
            report_period=300*second, profile=True)


    # --------- T5 ---------
    # freeze network and record Exc spikes
    # for cross correlations

    synee_scaling.active=False
    if tr.netw.config.iscl_active:
        synei_scaling.active=False
    strctplst.active=False
    strctplst_EI.active=False
    SynEE.stdp_active=0
    if tr.netw.config.istdp_active:
        SynEI.stdp_active=0

    set_active(GExc_spks)
    set_active(GInh_spks)


    net.run(tr.sim.T5, report='text',
            report_period=300*second, profile=True)
        
    SynEE_a.record_single_timestep()
    if tr.istdp_active and tr.synei_a_nrecpoints>0:
        SynEI_a.record_single_timestep()

    device.build(directory='builds/%.4d'%(tr.v_idx), clean=True,
                 compile=True, run=True, debug=False)


    # -----------------------------------------
    
    # save monitors as raws in build directory
    raw_dir = 'builds/%.4d/raw/'%(tr.v_idx)
    
    if not os.path.exists(raw_dir):
        os.makedirs(raw_dir)

    with open(raw_dir+'namespace.p','wb') as pfile:
        pickle.dump(namespace,pfile)   

    with open(raw_dir+'gexc_stat.p','wb') as pfile:
        pickle.dump(GExc_stat.get_states(),pfile)   
    with open(raw_dir+'ginh_stat.p','wb') as pfile:
        pickle.dump(GInh_stat.get_states(),pfile)   
        
    with open(raw_dir+'synee_stat.p','wb') as pfile:
        pickle.dump(SynEE_stat.get_states(),pfile)

    if tr.istdp_active:
        with open(raw_dir+'synei_stat.p','wb') as pfile:
            pickle.dump(SynEI_stat.get_states(),pfile)   
        
    with open(raw_dir+'synee_a.p','wb') as pfile:
        SynEE_a_states = SynEE_a.get_states()
        if tr.crs_crrs_rec:
            SynEE_a_states['i'] = list(SynEE.i)
            SynEE_a_states['j'] = list(SynEE.j)
        pickle.dump(SynEE_a_states,pfile)

    if tr.istdp_active and tr.synei_a_nrecpoints>0:
        with open(raw_dir+'synei_a.p','wb') as pfile:
            SynEI_a_states = SynEI_a.get_states()
            if tr.crs_crrs_rec:
                SynEI_a_states['i'] = list(SynEI.i)
                SynEI_a_states['j'] = list(SynEI.j)
            pickle.dump(SynEI_a_states,pfile)
        

    if tr.adjust_insertP:
        with open(raw_dir+'c_stat.p','wb') as pfile:
            pickle.dump(C_stat.get_states(),pfile)   

        with open(raw_dir+'insP_stat.p','wb') as pfile:
            pickle.dump(insP_stat.get_states(),pfile)

    if tr.adjust_EI_insertP:
        with open(raw_dir+'c_EI_stat.p','wb') as pfile:
            pickle.dump(C_EI_stat.get_states(),pfile)   

        with open(raw_dir+'insP_EI_stat.p','wb') as pfile:
            pickle.dump(insP_EI_stat.get_states(),pfile)   


    with open(raw_dir+'gexc_spks.p','wb') as pfile:
        pickle.dump(GExc_spks.get_states(),pfile)   
    with open(raw_dir+'ginh_spks.p','wb') as pfile:
        pickle.dump(GInh_spks.get_states(),pfile)

    if tr.external_mode=='poisson':
        with open(raw_dir+'pinp_spks.p','wb') as pfile:
            pickle.dump(PInp_spks.get_states(),pfile)

    with open(raw_dir+'gexc_rate.p','wb') as pfile:
        pickle.dump(GExc_rate.get_states(),pfile)
        if tr.rates_rec:
            pickle.dump(GExc_rate.smooth_rate(width=25*ms),pfile)   
    with open(raw_dir+'ginh_rate.p','wb') as pfile:
        pickle.dump(GInh_rate.get_states(),pfile)
        if tr.rates_rec:
            pickle.dump(GInh_rate.smooth_rate(width=25*ms),pfile)

    if tr.external_mode=='poisson':
        with open(raw_dir+'pinp_rate.p','wb') as pfile:
            pickle.dump(PInp_rate.get_states(),pfile)
            if tr.rates_rec:
                pickle.dump(PInp_rate.smooth_rate(width=25*ms),pfile)   


    # ----------------- add raw data ------------------------
    fpath = 'builds/%.4d/'%(tr.v_idx)

    from pathlib import Path

    Path(fpath+'turnover').touch()
    turnover_data = np.genfromtxt(fpath+'turnover',delimiter=',')    
    os.remove(fpath+'turnover')

    with open(raw_dir+'turnover.p','wb') as pfile:
        pickle.dump(turnover_data,pfile)


    Path(fpath+'turnover_EI').touch()
    turnover_EI_data = np.genfromtxt(fpath+'turnover_EI',delimiter=',')    
    os.remove(fpath+'turnover_EI')

    with open(raw_dir+'turnover_EI.p','wb') as pfile:
        pickle.dump(turnover_EI_data,pfile)   

        
    Path(fpath+'spk_register').touch()
    spk_register_data = np.genfromtxt(fpath+'spk_register',delimiter=',')
    os.remove(fpath+'spk_register')
    
    with open(raw_dir+'spk_register.p','wb') as pfile:
        pickle.dump(spk_register_data,pfile)

        
    Path(fpath+'spk_register_EI').touch()
    spk_register_EI_data = np.genfromtxt(fpath+'spk_register_EI',delimiter=',')
    os.remove(fpath+'spk_register_EI')
    
    with open(raw_dir+'spk_register_EI.p','wb') as pfile:
        pickle.dump(spk_register_EI_data,pfile)


    Path(fpath+'scaling_deltas').touch()
    scaling_deltas_data = np.genfromtxt(fpath+'scaling_deltas',delimiter=',')
    os.remove(fpath+'scaling_deltas')
    
    with open(raw_dir+'scaling_deltas.p','wb') as pfile:
        pickle.dump(scaling_deltas_data,pfile)

        
    Path(fpath+'scaling_deltas_EI').touch()
    scaling_deltas_data = np.genfromtxt(fpath+'scaling_deltas_EI',delimiter=',')
    os.remove(fpath+'scaling_deltas_EI')
    
    with open(raw_dir+'scaling_deltas_EI.p','wb') as pfile:
        pickle.dump(scaling_deltas_data,pfile)

               
     
    with open(raw_dir+'profiling_summary.txt', 'w+') as tfile:
        tfile.write(str(profiling_summary(net)))



    # --------------- cross-correlations ---------------------

    if tr.crs_crrs_rec:

        GExc_spks = GExc_spks.get_states()
        synee_a = SynEE_a_states
        wsize = 100*pq.ms

        for binsize in [1*pq.ms, 2*pq.ms, 5*pq.ms]: 

            wlen = int(wsize/binsize)

            ts, idxs = GExc_spks['t'], GExc_spks['i']
            idxs = idxs[ts>tr.T1+tr.T2+tr.T3+tr.T4]
            ts = ts[ts>tr.T1+tr.T2+tr.T3+tr.T4]
            ts = ts - (tr.T1+tr.T2+tr.T3+tr.T4)

            sts = [neo.SpikeTrain(ts[idxs==i]/second*pq.s,
                                  t_stop=tr.T5/second*pq.s) for i in
                   range(tr.N_e)]

            crs_crrs, syn_a = [], []

            for f,(i,j) in enumerate(zip(synee_a['i'], synee_a['j'])):
                if synee_a['syn_active'][-1][f]==1:

                    crs_crr, cbin = cch(BinnedSpikeTrain(sts[i],
                                                         binsize=binsize),
                                        BinnedSpikeTrain(sts[j],
                                                         binsize=binsize),
                                        cross_corr_coef=True,
                                        border_correction=True,
                                        window=(-1*wlen,wlen))

                    crs_crrs.append(list(np.array(crs_crr).T[0]))
                    syn_a.append(synee_a['a'][-1][f])


            fname = 'crs_crrs_wsize%dms_binsize%fms_full' %(wsize/pq.ms,
                                                            binsize/pq.ms)

            df = {'cbin': cbin, 'crs_crrs': np.array(crs_crrs),
                  'syn_a': np.array(syn_a), 'binsize': binsize,
                  'wsize': wsize, 'wlen': wlen}


            with open('builds/%.4d/raw/'%(tr.v_idx)+fname+'.p', 'wb') as pfile:
                pickle.dump(df, pfile)


        GInh_spks = GInh_spks.get_states()
        synei_a = SynEI_a_states
        wsize = 100*pq.ms

        for binsize in [1*pq.ms, 2*pq.ms, 5*pq.ms]: 

            wlen = int(wsize/binsize)

            ts_E, idxs_E = GExc_spks['t'], GExc_spks['i']
            idxs_E = idxs_E[ts_E>tr.T1+tr.T2+tr.T3+tr.T4]
            ts_E = ts_E[ts_E>tr.T1+tr.T2+tr.T3+tr.T4]
            ts_E = ts_E - (tr.T1+tr.T2+tr.T3+tr.T4)

            ts_I, idxs_I = GInh_spks['t'], GInh_spks['i']
            idxs_I = idxs_I[ts_I>tr.T1+tr.T2+tr.T3+tr.T4]
            ts_I = ts_I[ts_I>tr.T1+tr.T2+tr.T3+tr.T4]
            ts_I = ts_I - (tr.T1+tr.T2+tr.T3+tr.T4)

            sts_E = [neo.SpikeTrain(ts_E[idxs_E==i]/second*pq.s,
                                    t_stop=tr.T5/second*pq.s) for i in
                     range(tr.N_e)]

            sts_I = [neo.SpikeTrain(ts_I[idxs_I==i]/second*pq.s,
                                    t_stop=tr.T5/second*pq.s) for i in
                     range(tr.N_i)]

            crs_crrs, syn_a = [], []

            for f,(i,j) in enumerate(zip(synei_a['i'], synei_a['j'])):
                if synei_a['syn_active'][-1][f]==1:

                    crs_crr, cbin = cch(BinnedSpikeTrain(sts_I[i],
                                                         binsize=binsize),
                                        BinnedSpikeTrain(sts_E[j],
                                                         binsize=binsize),
                                        cross_corr_coef=True,
                                        border_correction=True,
                                        window=(-1*wlen,wlen))

                    crs_crrs.append(list(np.array(crs_crr).T[0]))
                    syn_a.append(synei_a['a'][-1][f])


            fname = 'EI_crrs_wsize%dms_binsize%fms_full' %(wsize/pq.ms,
                                                            binsize/pq.ms)

            df = {'cbin': cbin, 'crs_crrs': np.array(crs_crrs),
                  'syn_a': np.array(syn_a), 'binsize': binsize,
                  'wsize': wsize, 'wlen': wlen}


            with open('builds/%.4d/raw/'%(tr.v_idx)+fname+'.p', 'wb') as pfile:
                pickle.dump(df, pfile)


    # -----------------  clean up  ---------------------------
    shutil.rmtree('builds/%.4d/results/'%(tr.v_idx))
            

    # ---------------- plot results --------------------------

    #os.chdir('./analysis/file_based/')

    from analysis.overview_winh import overview_figure
    overview_figure('builds/%.4d'%(tr.v_idx), namespace)

    from analysis.synw_fb import synw_figure
    synw_figure('builds/%.4d'%(tr.v_idx), namespace)

    from analysis.synw_log_fb import synw_log_figure
    synw_log_figure('builds/%.4d'%(tr.v_idx), namespace)
    
    # from analysis.turnover_fb import turnover_figure
    # turnover_figure('builds/%.4d'%(tr.v_idx), namespace, fit=False)

    # from analysis.turnover_fb import turnover_figure
    # turnover_figure('builds/%.4d'%(tr.v_idx), namespace, fit=True)

          
