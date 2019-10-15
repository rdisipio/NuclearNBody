#!/usr/bin/env python3

import numpy as np

import dimod
import dwave_networkx as dnx
import neal

from dwave.system import EmbeddingComposite, FixedEmbeddingComposite, TilingComposite, DWaveSampler
from dwave_tools import get_embedding_with_short_chain, get_energy, anneal_sched_custom, qubo_quadratic_terms_from_np_array

num_reads = 1000
n = 3
strength = 5.

# linear
A = [ 1., 1., -1.]
# quadratic
B = [ 
    [1,-1,1], 
    [2,0,1], 
    [-1,2,0] 
    ]
# cubic
C = [ 
    [
        [1,1,1],
        [1,1,1],
        [1,1,1]
    ], 
    [
        [1,1,1],
        [1,1,1],
        [1,1,1]
    ], 
    [
        [1,1,1],
        [1,1,1],
        [1,1,1]
    ]
]

A = np.array(A)
B = np.array(B)
C = np.array(C)

#L = {} #{ 0:1}
#Q = { (0,1):1, (2,3):1 }

H = {}
#for i in range(n):
#    idx=(i,i)
#    H[idx] = A[i]

for i in range(n):
    for j in range(i+1,n):
        idx=(i,j)
        H[idx] = B[i][j]

for i in range(n):
    for j in range(i+1,n):
        for k in range(j+1,n):
            idx=(i,j,k)
            H[idx] = C[i][j][k]

print(H)
#bqm = dimod.BinaryQuadraticModel({0: 1, 1: -1, 2: .5},{(0, 1): .5, (1, 2): 1.5}, 1.4,dimod.SPIN)
#bqm = dimod.BinaryQuadraticModel( L, Q, 1.4,dimod.SPIN)

#Q = { (0,1):1, (2,3):1, (0,1,3):-1, (2,3,4):-1, (3,4):-1 }
bqm = dimod.make_quadratic(H, strength, dimod.BINARY)

#bqm = dimod.BinaryQuadraticModel( H, 0., dimod.SPIN)
print(bqm)

sampler = neal.SimulatedAnnealingSampler()

results = sampler.sample( bqm, num_reads=num_reads).aggregate()
print(results)