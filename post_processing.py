
import pickle
import numpy as np

from brian2.units import ms,mV,second,Hz

from analysis.methods.process_survival import extract_survival
from analysis.methods.process_turnover_pd import extract_lifetimes
from analysis.methods.resample_dA import resample_spk_register, \
                                         resample_scl_deltas

def post_process(tr):

    t_cut = tr.pp_tcut
    t_split = (tr.T2-t_cut)/2.

    if t_split/second > 0:

        with open('builds/%.4d/raw/turnover.p' %(tr.v_idx), 'rb') as pfile:
            turnover = pickle.load(pfile)

        full_t, ex_ids = extract_survival(turnover,
                                          tr.N_e,
                                          t_split=t_split,
                                          t_cut=t_cut)

        fpath = 'builds/%.4d/raw/survival_full_t.p' %(tr.v_idx)
        with open(fpath, 'wb') as pfile:
            out = {'t_split': t_split, 't_cut': t_cut,
                   'full_t': full_t, 'excluded_ids': ex_ids}
            pickle.dump(out, pfile)



    Tmax = tr.sim.T1+tr.sim.T2+tr.sim.T3

    lts_wthsrv, lts_dthonly, ex_ids = extract_lifetimes(turnover,
                                                        tr.N_e,
                                                        t_cut=t_cut,
                                                        Tmax=Tmax)

    bpath = 'builds/%.4d' %(tr.v_idx)
    for bin_w in [0.1*second, 0.5*second, 1*second, 10*second, 100*second]:

        bins = np.arange(bin_w/second,
                         Tmax/second+2*bin_w/second,
                         bin_w/second)

        f_add = 'bin%dcs' %(int(bin_w/second*10.))

        counts, edges = np.histogram(lts_wthsrv, bins=bins,
                                     density=True)
        centers = (edges[:-1] + edges[1:])/2.            
        
        with open(bpath+'/raw/lts_wthsrv_'+f_add+'.p', 'wb') as pfile:
            out = {'Tmax': Tmax, 't_cut': t_cut,
                   'counts': counts, 'excluded_ids': ex_ids,
                   'bins': bins, 'centers': centers, 'bin_w': bin_w}
            pickle.dump(out, pfile)


        counts, edges = np.histogram(lts_dthonly, bins=bins,
                                     density=True)
        centers = (edges[:-1] + edges[1:])/2.            

        with open(bpath+'/raw/lts_dthonly_'+f_add+'.p', 'wb') as pfile:
            out = {'Tmax': Tmax, 't_cut': t_cut,
                   'counts': counts, 'excluded_ids': ex_ids,
                   'bins': bins, 'centers': centers, 'bin_w': bin_w}
            pickle.dump(out, pfile)


    for nbins in [25,50,100,250,500,1000,2500,5000]:

        bins = np.logspace(np.log10(1),
                           np.log10((Tmax-t_cut)/second+0.5),
                           num=nbins)

        f_add = 'lognbin%d' %(nbins)

        counts, edges = np.histogram(lts_wthsrv, bins=bins,
                                     density=True)
        centers = (edges[:-1] + edges[1:])/2.            

        with open(bpath+'/raw/lts_wthsrv_'+f_add+'.p', 'wb') as pfile:
            out = {'Tmax': Tmax, 't_cut': t_cut,
                   'counts': counts, 'excluded_ids': ex_ids,
                   'bins': bins, 'centers': centers, 'nbins': nbins}
            pickle.dump(out, pfile)


        counts, edges = np.histogram(lts_dthonly, bins=bins,
                                     density=True)
        centers = (edges[:-1] + edges[1:])/2.            

        with open(bpath+'/raw/lts_dthonly_'+f_add+'.p', 'wb') as pfile:
            out = {'Tmax': Tmax, 't_cut': t_cut,
                   'counts': counts, 'excluded_ids': ex_ids,
                   'bins': bins, 'centers': centers, 'nbins': nbins}
            pickle.dump(out, pfile)



    from analysis.srvprb_all import srvprb_all_figure
    srvprb_all_figure('builds/%.4d'%(tr.v_idx))


    if tr.synEE_rec:
    
        with open(bpath+'/raw/spk_register.p', 'rb') as pfile:
            spkr = pickle.load(pfile)

        for resamp_dt in [1*second, 2*second, 5*second, 10*second]:

            dA_v, a_v = resample_spk_register(spkr, tr, resamp_dt)

            fname = 'stdp_dA_T%ds_bin%dms.p' %(int(tr.stdp_rec_T/second),
                                               int(resamp_dt/second*1000))

            with open(bpath+'/raw/'+fname, 'wb') as pfile:
                pickle.dump((dA_v, a_v), pfile)

                
    if tr.syn_scl_rec:
        
        with open(bpath+'/raw/scl_deltas.p', 'rb') as pfile:
            scl_deltas = pickle.load(pfile)

        for resamp_dt in [1*second, 2*second, 5*second, 10*second]:

            dA_v, a_v = resample_scl_deltas(scl_deltas, tr, resamp_dt)

            fname = 'scl_dA_T%ds_bin%dms.p' %(int(tr.scl_rec_T/second),
                                               int(resamp_dt/second*1000))

            with open(bpath+'/raw/'+fname, 'wb') as pfile:
                pickle.dump((dA_v, a_v), pfile)
