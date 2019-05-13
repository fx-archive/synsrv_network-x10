import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl

import os, git

from pypet import Environment, Trajectory
from pypet.brian2.network import NetworkManager
from pypet.brian2.parameter import Brian2Parameter, \
                                   Brian2MonitorResult

from .explored_params import explore_dict, name

from .add_parameters import add_params
from .xstrct_netw import run_net
from .post_processing import post_process


# control the number of cores to be used for computation
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--ncores", "-c", help="No. cores", nargs=1)
args = parser.parse_args()
ncores = int(args.ncores[0])
print("Using {:d} cores".format(ncores))

# get info on git repository
repo = git.Repo('./code/')
commit = repo.commit(None)

filename = os.path.join('data', name+'_'+str(commit)[:6]+'.hdf5')

# if not the first run, tr2 will be merged later
label = 'tr1'
    
env = Environment(trajectory=label,
                  add_time=False,
                  filename=filename,
                  continuable=False, # ??
                  lazy_debug=False,  # ??
                  multiproc=True,     
                  ncores=ncores,
                  use_pool=False, # likely not working w/ brian2
                  wrap_mode='QUEUE', # ??
                  overwrite_file=False, # ??
                  git_repository='./code/',
                  git_fail = True)


tr = env.trajectory

add_params(tr)

tr.f_add_parameter('mconfig.git.sha1', str(commit))
tr.f_add_parameter('mconfig.git.message', commit.message)

tr.f_explore(explore_dict)


def run_sim(tr):
    run_net(tr)
    post_process(tr)
    
env.run(run_sim)



                  


