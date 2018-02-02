
import numpy as np


def generate_connections(N_tar, N_src, p, same=False):
    ''' 
    connect source to target with probability p
      - if populations SAME, avoid self-connection
      - if not SAME, connect any to any avoididing multiple
    return list of sources i and targets j
    '''
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

def generate_full_connectivity(N, same):
    if not same==True:
        raise NotImplementedError
    i = []
    j = []
    for k in range(N):
        i.extend([k]*(N-1))
        targets = list(range(N))
        del targets[k]
        j.extend(targets)

    assert len(i)==len(j)
    return np.array(i), np.array(j)
