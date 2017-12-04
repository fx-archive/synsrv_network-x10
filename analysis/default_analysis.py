
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from matplotlib import rc

rc('text', usetex=True)
pl.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',    # helvetica font
    r'\usepackage{sansmath}',   # math-font matching helvetica
    r'\sansmath'                # actually tell tex to use it!
    r'\usepackage{siunitx}',    # micro symbols
    r'\sisetup{detect-all}',    # force siunitx to use the fonts
]  

import numpy as np

from brian2.units import mV, ms, second

from pypet.trajectory import Trajectory
from pypet.brian2.parameter import Brian2Parameter, Brian2MonitorResult

tr = Trajectory(#name='stdp_scl_it_strct',
                name='tr1',
                add_time=False, 
                #filename='../data/stdp_scl_it_strct.hdf5',
                filename='../data/ns2.hdf5',
                dynamic_imports=[Brian2MonitorResult, Brian2Parameter])


# pypet.pypetconstants.LOAD_NOTHING  --> 0
# pypet.pypetconstants.LOAD_SKELETON --> 1
# pypet.pypetconstants.LOAD_DATA     --> 2
tr.f_load(load_parameters=2, load_derived_parameters=2, load_results=1)
tr.v_auto_load = True

def raster_plot(ax, tr, crun='run_00000000'):
    df_e, df_i = tr.crun.GExc_spks, tr.crun.GInh_spks
    ax.plot(df_e.t/ms, df_e.i, marker='.', color='blue',
            markersize=.5, linestyle='None')
    ax.plot(df_i.t/ms, df_i.i+tr.N_e, marker='.', color='red',
            markersize=.5, linestyle='None')
    ax.set_title('Spike Raster')
    ax.set_xlabel('time [ms]')

def firing_rate_distribution(ax, tr, crun='run_00000000', bins=25):
    df_e, df_i = tr.crun.GExc_spks, tr.crun.GInh_spks
    ax.hist(np.bincount(df_e.i)/(tr.T/second), bins=bins)
    ax.set_title('Firing Rate Distribution')

def voltage_traces(ax, tr, crun='run_00000000'):
    df = tr.crun.GExc_stat 
    for i in df.record:
        ax.plot(df.V[i]/mV)
    ax.set_ylim(tr.Vr_e/mV-5, tr.Vt_e/mV+5)
    ax.set_title('Membrane Voltage Traces')

def synaptic_weight_distribution(ax, tr, crun='run_00000000', bins=50):
    df = tr.crun.SynEE_a
    bins = np.linspace(0, tr.amax, num=bins)
    ax.hist(np.nan_to_num(df.a.flatten()), bins=bins)
    ax.set_title('Synaptic Weight Distribution')

def single_synapse_weight_distribution(ax, tr, crun='run_00000000', bins=50):
    df = tr.crun.SynEE_a
    bins = np.linspace(0, tr.amax, num=bins)
    ax.hist(np.nan_to_num(df.a.flatten()), bins=bins)

def membrane_threshold_traces(ax, tr, crun='run_00000000'):
    df = tr.crun.GExc_stat
    for i in df.record:
        ax.plot(df.Vt[i]/mV)
    ax.set_title('{:.3f} ms mV, {:.0f} ms'.format(tr.eta_ip/(ms*mV),tr.it_dt/ms))
    #ax.set_ylim(tr.Vr_e/mV-5, tr.Vt_e/mV+5)

def membrane_threshold_distribution(ax, tr, crun='run_00000000'):
    df = tr.crun.GExc_vts
    for i,t in enumerate(df.t):
        ax.hist(df.Vt[:,i]/mV)
    ax.set_title('Membrane Threshold Distribution')
    ax.set_ylabel('threshold [mV]')


def default_analysis_figure(tr, crun='run_00000000'):
    
    pl.close()
    fig = pl.figure()

    # ----- 3x4 -----
    s = 150
    fig.set_size_inches(1920/s,1080/s)
    ax1 = pl.subplot2grid((3, 4), (0, 0))
    ax2 = pl.subplot2grid((3, 4), (0, 1))
    ax3 = pl.subplot2grid((3, 4), (0, 2))
    ax4 = pl.subplot2grid((3, 4), (0, 3))
    ax5 = pl.subplot2grid((3, 4), (1, 0))
    ax6 = pl.subplot2grid((3, 4), (1, 1))
    ax7 = pl.subplot2grid((3, 4), (1, 2))
    ax8 = pl.subplot2grid((3, 4), (1, 3))
    ax9 = pl.subplot2grid((3, 4), (2, 0))
    ax10 = pl.subplot2grid((3, 4), (2, 1))
    ax11 = pl.subplot2grid((3, 4), (2, 2))
    ax12 = pl.subplot2grid((3, 4), (2, 3))
    axs=[ax1,ax2,ax3,ax4,ax5,ax6, ax7, ax8, \
         ax9,ax10,ax11,ax12]

    assert(tr.v_crun == crun)

    raster_plot(ax1, tr, crun=crun)
    firing_rate_distribution(ax2, tr, crun=crun, bins=25)
    voltage_traces(ax3, tr, crun=crun)
    synaptic_weight_distribution(ax5, tr, crun=crun)
    membrane_threshold_distribution(ax9, tr, crun=crun)
    membrane_threshold_traces(ax10, tr, crun=crun)
            
    pl.tight_layout()
    pl.savefig("default_analysis_output/50sec_dfa_{:s}.png".format(crun), dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    # tr.v_idx = 2
    # crun = tr.v_crun
    # default_analysis_figure(tr,crun)
    
    for crun in tr.f_iter_runs():
        default_analysis_figure(tr, crun)
