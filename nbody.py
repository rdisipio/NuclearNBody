#!/usr/bin/env python3

import random
import itertools
import numpy as np
from decimal import Decimal
from collections import Counter

import dimod
import dwave_networkx as dnx
import minorminer
import neal

from dwave.system import EmbeddingComposite, FixedEmbeddingComposite, AutoEmbeddingComposite, DWaveSampler
from dimod.reference.composites import HigherOrderComposite
from dwave_tools import get_embedding_with_short_chain, get_energy, anneal_sched_custom, qubo_quadratic_terms_from_np_array

# see Reverse Annealing:
# https://docs.dwavesys.com/docs/latest/c_fd_ra.html


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def onebody(i, n, l):
    '''Expectation value for the one body part, 
    Harmonic oscillator in three dimensions'''

    omega = 10.0

    return omega * ( 2*n[i] + l[i] + 1.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def make_hamiltonian( n_so, quantum_numbers, nninteraction ):
    '''Creates the Hamiltonian'''

    N = quantum_numbers['n']
    L = quantum_numbers['l']
    J = quantum_numbers['j']
    MJ = quantum_numbers['mj']
    TZ = quantum_numbers['tz']

    # one-body hamiltonian
    # kinetic energy and external field
    # T = sum_ij T_ij a_i^† a_j
    #T = np.zeros( [n_so,n_so] )
    #T = { (0,1):1, (1,3):-1, (2,4):2, (2,3):-1.5 }
    T = {}
    # no self-loops allowed ???
    for i in range(n_so):
        idx = ( str(i),str(i) )
        T[idx] = onebody( i, N, L )
        print( idx, T[idx])

    # two-body hamiltonian
    # conservation laws/selection rules 
    # strongly restrict the elements
    # V = sum_ijkl V_ijkl a_i^† a_j^† a_k a_l
    #V = np.zeros( [n_so,n_so,n_so,n_so] )
    #V = { (0,1,2,3):1, (0,1,0,3):1, (1,0,2,3):-2 }

    V = {}
    for i in range(n_so):
        for j in range(n_so):
            if L[i] != L[j] and J[i] != J[j] and MJ[i] != MJ[j] and TZ[i] != TZ[j]: continue

            for k in range(n_so):
                for l in range(n_so):

                    if (MJ[i]+MJ[k]) != (MJ[j]+MJ[l]) and (TZ[i]+TZ[k]) != (TZ[j]+TZ[l]): continue

                    if nninteraction[i][j][k][l] == 0.: continue

                    idx = ( str(i),str(j),str(k),str(l))
                    V[idx] = nninteraction[i][j][k][l]

                    print(idx, ':', V[idx], ',')

    # let's put something in by hand according to
    # H = sum_ij T_ij a_i^† a_j + sum_ijkl V_ijkl a_i^† a_j^† a_k a_l

    #V = {
    #    (0,1,0,2): -5,
    #}

    H = {**T, **V}

    #H = { 
    #    ('0','0'):-13.6, ('1','1'):-3.4, ('2','2'):-1.5,
    #    ('1','0'):10.2, ('2','0'):12.1, ('2','1'):1.9, 
    #    ('0','1','1','2'):5.0,
    #}
    
    print("INFO: Hamiltonian:")
    print(H)

    return H

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def dirac(x):
    sx = [ str(a) for a in x ]
    s = "".join(sx)
    s = "|%s>" % s
    return s

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if __name__ == '__main__':

    n_so = 5 # number of spin orbitals
    n_p = 2  # number of particles actually present
    num_reads = 100
    max_evals = 2

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

    Q = make_hamiltonian( n_so, quantum_numbers, nninteraction )

    # Strength of the reduction constraint. 
    # Insufficient strength can result in the binary quadratic model
    # not having the same minimizations as the polynomial.
    strength = 5.0 

    # see: https://docs.ocean.dwavesys.com/projects/dimod/en/latest/reference/generated/dimod.higherorder.utils.make_quadratic.html#dimod.higherorder.utils.make_quadratic
    bqm = dimod.make_quadratic(Q, strength, dimod.BINARY)
    print("INFO: simplified binary quadratic model:")
    print(bqm)

    schedule = [(0.0, 1.0), (2.0, 0.7), (98.0, 0.7), (100.0, 1.0)]
    print("INFO: reverse annealing schedule:")
    print(schedule)
    
    neal_sampler = neal.SimulatedAnnealingSampler()
    hardware_sampler = DWaveSampler()
    source = bqm.quadratic
    embedding = minorminer.find_embedding(source, hardware_sampler.edgelist)
    #qpu_sampler = FixedEmbeddingComposite(hardware_sampler, embedding)
    qpu_sampler = dimod.HigherOrderComposite(EmbeddingComposite(DWaveSampler()))


    #processor = dnx.chimera_graph(16, 16, 4).edges()
    #embedding = minorminer.find_embedding(source, processor)

    print("INFO: embedding:")
    print(embedding)

    # set aside the indices of the additional qubits generated by make_quadratic:
    qubit_indices = []
    embedding_keys = list(embedding.keys())
    print(embedding_keys)
    for i in range(len(embedding_keys)):
        if len(embedding_keys[i]) == 1:
            qubit_indices += [i]
    print("INFO: indices of valid qubits:")
    print(qubit_indices)

    #qpu_sampler = FixedEmbeddingComposite(hardware_sampler, embedding)


    initial_states = set( itertools.permutations( [1]*n_p + [0]*(n_so-n_p), n_so )  )
    initial_states = list( zip( np.arange(len(initial_states)), initial_states) )
    print("INFO: possible initial states with %i particles in %i spin orbitals:" % (n_p, n_so) )
    for s in initial_states:
        print( "%-3i ="%s[0], dirac(s[1]))

    for initial_state in initial_states:
        counter = Counter()
        
        #print("DEBUG: initial state:", initial_state[1])
        s0 = {}
        for j in range(n_so):
            s0[str(j)] = initial_state[1][j]
        for v in bqm:
            if v not in initial_state:
                s0[v] = random.choice((0, 1))
        #print("DEBUG: s0=",s0)

        for itrial in range(max_evals):

            reverse_anneal_params = dict(anneal_schedule=schedule,
                                 initial_state=s0,
                                 reinitialize_state=True)

            solver_parameters = {
                'num_reads': num_reads,
                'auto_scale': True,
                **reverse_anneal_params,
            }

            # CPU:
            results = neal_sampler.sample(bqm, **solver_parameters).aggregate()

            # QPU:
            #results = qpu_sampler.sample_hubo(Q, **solver_parameters )
            
            #this_energy = this_result.energy
            #this_q = np.array(list(this_result.sample.values()))

            #print("Results:")
            #print(results)

            # this also contains the state of ancilla qubits
            best_fit = list( results.first.sample.values() )
            #print("DEBUG: best-fit:", best_fit)
            gs = dirac(best_fit)

            # CPU/neal: remove elements corresponding to ancilla qbits
            best_fit_x = [ best_fit[j] for j in qubit_indices ]
            gs = dirac(best_fit_x)

            counter[gs] += 1

            #print("INFO: ground state (%i/%i):" % (itrial, max_evals) )
            #print(gs)

        #x0 = [ a[1] for a in initial_state ]
        #s0 = dirac( x0 )
        s0 = dirac( initial_state[1] )
        print("INFO: counters for initial state", s0)
        for state, c in counter.items():
            print( state, ":", c)
        print("-------")
    