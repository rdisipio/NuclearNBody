#!/usr/bin/env python3

import numpy as np

import dimod
import dwave_networkx as dnx
import neal

from dwave.system import EmbeddingComposite, FixedEmbeddingComposite, TilingComposite, DWaveSampler
from dwave_tools import get_embedding_with_short_chain, get_energy, anneal_sched_custom, qubo_quadratic_terms_from_np_array

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def make_hamiltonian( n_sp = 4 ):
    '''Creates the Hamiltonian'''

    # one-body hamiltonian
    # kinetic energy and external field
    # T = sum_ij T_ij a_i^† a_j
    #T = np.zeros( [n_sp,n_sp] )
    T = { (0,1):1, (1,3):-1, (2,4):2, (2,3):-1.5 }

    # two-body hamiltonian
    # conservation laws/selection rules 
    # strongly restrict the elements
    # V = sum_ijkl V_ijkl a_i^† a_j^† a_k a_l
    #V = np.zeros( [n_sp,n_sp,n_sp,n_sp] )
    V = { (0,1,2,3):1, (0,1,0,3):1, (1,0,2,3):-2 }

    # let's put something in by hand according to
    # H = sum_ij T_ij a_i^† a_j + sum_ijkl V_ijkl a_i^† a_j^† a_k a_l

    H = {**T, **V}

    return H

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def ground_state(x):
    q = np.array( [ str(x) for x in list( x.sample.values() ) ] )
    s = "".join(q)
    s = "|ψ> = %s|0>" % s
    return s

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

    # number of single particles (nucleons)
    n_sp = 4
    num_reads = 1000

    Q = make_hamiltonian( n_sp )

    # Strength of the reduction constraint. 
    # Insufficient strength can result in the binary quadratic model
    # not having the same minimizations as the polynomial.
    strength = 5.0 

    # see: https://docs.ocean.dwavesys.com/projects/dimod/en/latest/reference/generated/dimod.higherorder.utils.make_quadratic.html#dimod.higherorder.utils.make_quadratic
    bqm = dimod.make_quadratic(Q, strength, dimod.BINARY)

    print(bqm)

    sampler = neal.SimulatedAnnealingSampler()

    results = sampler.sample( bqm, num_reads=num_reads).aggregate()

    best_fit = results.first
    gs = ground_state(best_fit)

    print("INFO: ground state:")
    print(gs)

    