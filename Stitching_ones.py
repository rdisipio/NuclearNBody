# -*- coding: utf-8 -*-

#import neal
import dwavebinarycsp
import operator as op
import pyqubo
import dimod
import networkx
import neal
from helper_functions import *
from datetime import datetime
import dwave_networkx as dnx
import matplotlib.pyplot as plt
import numpy as np
startTime = datetime.now()


####################################    
lin_one = {} 

#i1, j1 reprresent succsesful transformations on |Ψ_1> = |1>

num_ones_direct_product = 100 #Number of times (n) in ∏_[i = 2]^[n] (⨂_i |Ψ_i>) 
var_list = []
for num in range(num_ones_direct_product):
    num += 1 #run from 1 to 2
    
    name = 'i' + str(num)
    var_list += [name]
    lin_one.update({name : -2}) #Create i1 : -2
    
    name = 'j' + str(num)
    var_list += [name]
    lin_one.update({'j' + str(num) : -2}) #Create j1 : -2

quad_one = {}
str_one = 2.0

bqm_one = dimod.BinaryQuadraticModel(lin_one, quad_one,
                                      str_one, 'BINARY')

#Solvers
#
#sampler = dimod.ExactSolver()
sampler = neal.SimulatedAnnealingSampler()


resp = sampler.sample(bqm_one, num_reads=1)
result_list = []
for sample, energy in resp.data(['sample', 'energy']):
    result_list += [sample]

transformation_matrix = []

for var in var_list:
    transformation_matrix += [sample[var]]
    
transformation_matrix = np.array(transformation_matrix)
transformation_matrix = np.reshape(transformation_matrix,
                                 (2, num_ones_direct_product))
mult_matrix = np.diagflat(np.arange(1,
                                    num_ones_direct_product + 1, 1))
#transformation_matrix = np.matmul(transformation_matrix, mult_matrix)
print(transformation_matrix)