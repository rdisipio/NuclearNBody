# -*- coding: utf-8 -*-

#import neal
import dwavebinarycsp
import operator as op
import pyqubo
import dimod
import networkx
from helper_functions import *
from datetime import datetime
import dwave_networkx as dnx
import matplotlib.pyplot as plt
import numpy as np
startTime = datetime.now()


####################################    
lin_zero = {'i1': -6.0, 'i2': 2.0, 'j1': 2.0, 'j2': -6.0}
quad_zero = {('i1', 'i2'): 4.0, ('i1', 'j1'): 4.0, ('i1', 'j2'): -4.0,
             ('i2', 'j1'): -4.0, ('i2', 'j2'): 4.0, ('j1', 'j2'): 4.0}
str_zero = 16.0

lin_one = {'i1': -2, 'j1': -2}
quad_one = {}
str_one = 2.0

bqm_zero = dimod.BinaryQuadraticModel(lin_zero, quad_zero,
                                      str_zero, 'BINARY')

bqm_one = dimod.BinaryQuadraticModel(lin_one, quad_one,
                                      str_one, 'BINARY')

#Solvers
#
resp = dimod.ExactSolver().sample(bqm_one)

for sample, energy in resp.data(['sample', 'energy']):
    print(sample, energy, )
