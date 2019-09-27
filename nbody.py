#!/usr/bin/env python3

import numpy as np
from decimal import Decimal

import dimod
import dwave_networkx as dnx
import neal

from dwave.system import EmbeddingComposite, FixedEmbeddingComposite, TilingComposite, DWaveSampler
from dwave_tools import get_embedding_with_short_chain, get_energy, anneal_sched_custom, qubo_quadratic_terms_from_np_array

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def onebody(i, n, l):
    '''Expectation value for the one body part, 
    Harmonic oscillator in three dimensions'''

    omega = 10.0

    return omega * ( 2*n[i] + l[i] + 1.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def make_hamiltonian( n_sp, quantum_numbers, nninteraction ):
    '''Creates the Hamiltonian'''

    N = quantum_numbers['n']
    L = quantum_numbers['l']
    J = quantum_numbers['j']
    MJ = quantum_numbers['mj']
    TZ = quantum_numbers['tz']

    # one-body hamiltonian
    # kinetic energy and external field
    # T = sum_ij T_ij a_i^† a_j
    #T = np.zeros( [n_sp,n_sp] )
    #T = { (0,1):1, (1,3):-1, (2,4):2, (2,3):-1.5 }
    T = {}
    # no self-loops allowed ???
    for i in range(n_sp):
        idx = (i,i)
        T[idx] = onebody( i, N, L )
        print( idx, T[idx])

    # two-body hamiltonian
    # conservation laws/selection rules 
    # strongly restrict the elements
    # V = sum_ijkl V_ijkl a_i^† a_j^† a_k a_l
    #V = np.zeros( [n_sp,n_sp,n_sp,n_sp] )
    #V = { (0,1,2,3):1, (0,1,0,3):1, (1,0,2,3):-2 }

    V = {}
    for i in range(n_sp):
        for j in range(n_sp):
            if L[i] != L[j] and J[i] != J[j] and MJ[i] != MJ[j] and TZ[i] != TZ[j]: continue

            for k in range(n_sp):
                for l in range(n_sp):

                    if (MJ[i]+MJ[k]) != (MJ[j]+MJ[l]) and (TZ[i]+TZ[k]) != (TZ[j]+TZ[l]): continue

                    if nninteraction[i][j][k][l] == 0.: continue

                    idx = (i,j,k,l)
                    V[idx] = nninteraction[i][j][k][l]

                    print(idx, ':', V[idx], ',')

    # let's put something in by hand according to
    # H = sum_ij T_ij a_i^† a_j + sum_ijkl V_ijkl a_i^† a_j^† a_k a_l

    #V = {
    #    (0,1,0,2): -5,
    #}

    H = {**T, **V}
    
    print("INFO: Hamiltonian:")
    print(H)

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
    n_sp = 3
    num_reads = 1000

    # see:
    # https://github.com/ManyBodyPhysics/LectureNotesPhysics/blob/master/Programs/Chapter8-programs/python/hfnuclei.py
    quantum_numbers = {
        'index' : [],
        'n' : [],
        'l' : [],
        'j' : [],
        'mj' : [],
        'tz' : [],
    }
    spOrbitals = 0
    with open("qnumbers.dat", "r") as qnumfile:
        for line in qnumfile:
            nums = line.split()
            if len(nums) != 0:
                quantum_numbers['index'].append(int(nums[0]))
                quantum_numbers['n'].append(int(nums[1]))
                quantum_numbers['l'].append(int(nums[2]))
                quantum_numbers['j'].append(int(nums[3]))
                quantum_numbers['mj'].append(int(nums[4]))
                quantum_numbers['tz'].append(int(nums[5]))
                spOrbitals += 1

    nninteraction = np.zeros([spOrbitals, spOrbitals, spOrbitals, spOrbitals])
    with open("nucleitwobody.dat", "r") as infile:
        for line in infile:
            number = line.split()
            a = int(number[0]) - 1
            b = int(number[1]) - 1
            c = int(number[2]) - 1
            d = int(number[3]) - 1
			#print a, b, c, d, float(l[4])
            nninteraction[a][b][c][d] = Decimal(number[4])

    Q = make_hamiltonian( n_sp, quantum_numbers, nninteraction )

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

    