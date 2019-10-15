#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 10:28:45 2019

@author: sahibjeetsingh
"""

import dimod
import dwave_networkx as dnx
import neal
import pyqubo

lin = {'i1': 0, 'i2': 0, 'j1': 0, 'j2': 0}
quad = { ('i1', 'i2'): 1.0, ('i1', 'j1'): 0.0, ('i1', 'j2'): -1.0,
        ('i2', 'j1'): 0.0, ('i2', 'j2'): -1.0, ('j1', 'j2'): 1.0}

dimod.BinaryQuadraticModel(lin, quad, 0, 'BINARY')
resp = dimod.ExactSolver().sample(bqm)
print(resp)