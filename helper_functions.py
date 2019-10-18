#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 17:15:01 2019

@author: sahibjeetsingh
"""

# Helper functions for nbody
####################################
# Create list of varriables like i1, i2, ... i_n_so, j1, j2...
# This must be fed to the add_constraint and represents a qubit
# i1, i2, ... i_n_so are the index of the creation
# operator in binnary

def init_orbital_list(n_so): # n_so is # of spin orbitals
    idx = 0
    so_list = []
    idx2str = { 0: 'i', 1: 'j'}
    while idx != 2:
        for i in range(n_so):
            so_list += [idx2str[idx] + str(i+1)]
        idx += 1
    return so_list

def init_orbital_list_2nd(n_so): # n_so is # of spin orbitals
    idx = 0
    so_list = []
    idx2str = { 0: 'i', 1: 'j', 2: 'k', 3: 'l'}
    while idx != 4:
        for i in range(n_so):
            so_list += [idx2str[idx] + str(i+1)]
        idx += 1
    return so_list

####################################
# Function to represent the anhilation operator
# Fed a state, if input is == False, returns False,
# not necesary here but needed in Create().

def Anhilate(state_input, anhil):
    if state_input == False:
        return False
    clone_input = state_input.copy() # copy input state so we can change it and output it
    if clone_input[anhil - 1] == 0: # a|0> = 0, anhil - 1 since anhil = j and j runs from 1 ro n_so
        return False
    else:
        clone_input[anhil - 1] = 0 # a|1> = |0>
        
        return clone_input  # clone_input = a|Ψ_i>
#
def Create(state_input, create):
    if state_input == False: # since we anhilate first then create (a†a|Ψ>), if a|Ψ> = 0 = False then 
                             # we simply say that a†0 = 0 -> return False
        return False
    
    clone_input = state_input.copy() # copy input state so we can change it and output it
    if clone_input[create - 1] == 1: # a†|1> = 0
        return False
    else:
        clone_input[create - 1] = 1 # a†|0> = |1>
        return clone_input # clone_input = a†|Ψ_i>
