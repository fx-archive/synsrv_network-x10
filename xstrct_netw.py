
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


from .cpp_methods import syn_scale, record_turnover, record_spk


    
def run_net(tr):

    # prefs.codegen.target = 'numpy'
    # prefs.codegen.target = 'cython'
    if tr.n_threads > 1:
        prefs.devices.cpp_standalone.openmp_threads = tr.n_threads
        
    set_device('cpp_standalone', directory='./builds/%.4d'%(tr.v_idx),
               build_on_run=False)

    print("Started process with id ", str(tr.v_idx))

    T = tr.T1 + tr.T2 + tr.T3 + tr.T4

    namespace = tr.netw.f_to_dict(short_names=True, fast_access=True)
    namespace['idx'] = tr.v_idx

    defaultclock.dt = tr.netw.sim.dt

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

    # set initial thresholds fixed, init. potentials uniformly distrib.

    if tr.external_mode=='memnoise':
        GExc.mu, GInh.mu = tr.mu_e, tr.mu_i
        GExc.sigma, GInh.sigma = tr.sigma_e, tr.sigma_i
        
    GExc.Vt, GInh.Vt = tr.Vt_e, tr.Vt_i
    GExc.V , GInh.V  = np.random.uniform(tr.Vr_e/mV, tr.Vt_e/mV,
                                         size=tr.N_e)*mV, \
                       np.random.uniform(tr.Vr_i/mV, tr.Vt_i/mV,
                                         size=tr.N_i)*mV


    synEE_pre_mod = mod.synEE_pre
    synEE_post_mod = mod.synEE_post

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
    
    

    if tr.stdp_active:
        synEE_pre_mod  = '''%s 
                            %s''' %(synEE_pre_mod, mod.synEE_pre_STDP)
        synEE_post_mod = '''%s 
                            %s''' %(synEE_post_mod, mod.synEE_post_STDP)

    
    if tr.synEE_rec:
        synEE_pre_mod  = '''%s 
                            %s''' %(synEE_pre_mod, mod.synEE_pre_rec)
        synEE_post_mod = '''%s 
                            %s''' %(synEE_post_mod, mod.synEE_post_rec)

        
    # E<-E advanced synapse model, rest simple
    SynEE = Synapses(target=GExc, source=GExc, model=tr.synEE_mod,
                     on_pre=synEE_pre_mod, on_post=synEE_post_mod,
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
        srcs_full, tars_full = generate_full_connectivity(tr.N_e, same=True)
        SynEE.connect(i=srcs_full, j=tars_full)
        SynEE.syn_active = 0


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

        
    SynEE.insert_P = tr.insert_P
    SynEE.p_inactivate = tr.p_inactivate
    SynEE.stdp_active=1

    # make randomly chosen synapses active at beginning
    rs = np.random.uniform(size=tr.N_e*(tr.N_e-1))
    initial_active = (rs < tr.p_ee).astype('int')
    initial_a = initial_active * tr.a_ee
    SynEE.syn_active = initial_active
    SynEE.a = initial_a
    
        
    # synaptic scaling
    if tr.netw.config.scl_active:
        SynEE.summed_updaters['Asum_post']._clock = Clock(
            dt=tr.dt_synEE_scaling)
        synscaling = SynEE.run_regularly(tr.synEE_scaling,
                                         dt=tr.dt_synEE_scaling, when='end')

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


            
    # -------------- recording ------------------        

    #run(tr.sim.preT)

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

    
    GExc_spks = SpikeMonitor(GExc)    
    GInh_spks = SpikeMonitor(GInh)
    
    if tr.external_mode=='poisson':
        PInp_spks = SpikeMonitor(PInp)

    GExc_rate = PopulationRateMonitor(GExc)
    GInh_rate = PopulationRateMonitor(GInh)

    if tr.external_mode=='poisson':
        PInp_rate = PopulationRateMonitor(PInp)


    if tr.synee_a_nrecpoints==0:
        SynEE_a_dt = 10*tr.sim.T2
    else:
        SynEE_a_dt = tr.sim.T2/tr.synee_a_nrecpoints
    SynEE_a = StateMonitor(SynEE, ['a','syn_active'],
                           record=range(tr.N_e*(tr.N_e-1)),
                           dt=SynEE_a_dt,
                           when='end', order=100)

    if tr.external_mode=='poisson':
        net = Network(GExc, GInh, PInp, sPN, sPNInh, SynEE, SynEI, SynIE, SynII,
                      GExc_stat, GInh_stat, SynEE_stat, SynEE_a,
                      GExc_spks, GInh_spks, PInp_spks, GExc_rate, GInh_rate,
                      PInp_rate, PInp_inh)
    else:
        net = Network(GExc, GInh, SynEE, SynEI, SynIE, SynII,
                      GExc_stat, GInh_stat, SynEE_stat, SynEE_a,
                      GExc_spks, GInh_spks, GExc_rate, GInh_rate)


    spks_recorders = [GExc_spks, GInh_spks]
    if tr.external_mode=='poisson':
        spks_recorders.append(PInp_spks)
        
    stat_recorders = [SynEE_stat, GExc_stat, GInh_stat]
    
    rate_recorders = [GExc_rate, GInh_rate]
    if tr.external_mode=='poisson':
        rate_recorders.append(PInp_rate)

    for rcc in spks_recorders:
        rcc.active=False
    for rcc in stat_recorders:
        rcc.active=False
    for rcc in rate_recorders:
        rcc.active=False
    
    for rcc in stat_recorders:
        rcc.active=True
    if tr.rates_rec:
        for rcc in rate_recorders:
            rcc.active=True
    if tr.spks_rec:
        for spr in spks_recorders:
            spr.active=True

    device.insert_code('main', '''
    cout << "Testing direct insertion of code." << endl;
    ''')
       
    net.run(tr.sim.T1, report='text',
            report_period=300*second, profile=True)

    for rcc in spks_recorders:
        rcc.active=False
    for rcc in stat_recorders:
        rcc.active=False
    for rcc in rate_recorders:
        rcc.active=False


    # print('''Hack solution: Simulate single timestep 
    #          to avoid missing simulation chunks''')
    # net.run(tr.dt)

    # for time_step in range(int(tr.sim.T2/(10000*second))):
    #     net.run(10000*second, report='text')

    net.run(tr.sim.T2, report='text', report_period=300*second,
            profile=True)
        

    for rcc in stat_recorders:
        rcc.active=True
    if tr.rates_rec:
        for rcc in rate_recorders:
            rcc.active=True
    if tr.spks_rec:
        for spr in spks_recorders:
            spr.active=True

    net.run(tr.sim.T3, report='text', report_period=300*second,
            profile=True)


    # freeze network
    synscaling.active=False
    strctplst.active=False
    SynEE.stdp_active=0
    
    for rcc in stat_recorders:
        rcc.active=False
    for rcc in rate_recorders:
        rcc.active=False
    for spr in spks_recorders:
        spr.active=True

    if tr.external_mode=='poisson':
        PInp_rate.active=False        
    
    net.run(tr.sim.T4, report='text', report_period=300*second)
        
    SynEE_a.record_single_timestep()
 

    device.build(directory='builds/%.4d'%(tr.v_idx), clean=True,
                 compile=True, run=True, debug=False)
    
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
    with open(raw_dir+'synee_a.p','wb') as pfile:
        SynEE_a_states = SynEE_a.get_states()
        if tr.crs_crrs_rec:
            SynEE_a_states['i'] = list(SynEE.i)
            SynEE_a_states['j'] = list(SynEE.j)
        pickle.dump(SynEE_a_states,pfile)

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
    
    Path(fpath+'spk_register').touch()
    spk_register_data = np.genfromtxt(fpath+'spk_register',delimiter=',')
    os.remove(fpath+'spk_register')
    
    with open(raw_dir+'spk_register.p','wb') as pfile:
        pickle.dump(spk_register_data,pfile)
     
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
            idxs = idxs[ts>tr.T1+tr.T2+tr.T3]
            ts = ts[ts>tr.T1+tr.T2+tr.T3]
            ts = ts - (tr.T1+tr.T2+tr.T3)

            sts = [neo.SpikeTrain(ts[idxs==i]/second*pq.s,
                                  t_stop=tr.T4/second*pq.s) for i in
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


    # -----------------  clean up  ---------------------------
    shutil.rmtree('builds/%.4d/results/'%(tr.v_idx))
            

    # ---------------- plot results --------------------------

    #os.chdir('./analysis/file_based/')

    from analysis.overview_fb import overview_figure
    overview_figure('builds/%.4d'%(tr.v_idx), namespace)

    from analysis.synw_fb import synw_figure
    synw_figure('builds/%.4d'%(tr.v_idx), namespace)

    from analysis.synw_log_fb import synw_log_figure
    synw_log_figure('builds/%.4d'%(tr.v_idx), namespace)

    from analysis.turnover_fb import turnover_figure
    turnover_figure('builds/%.4d'%(tr.v_idx), namespace, fit=False)

    # from analysis.turnover_fb import turnover_figure
    # turnover_figure('builds/%.4d'%(tr.v_idx), namespace, fit=True)

          

