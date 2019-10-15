# -*- coding: utf-8 -*-

#import neal
import dwavebinarycsp
import operator as op
import pyqubo
import dimod
import networkx
from helper_functions import *

####################################    
####################################
# Converts binnary i1 ... i_n_so and j1 ... j_n_so strings to int
# Input arg is the entire list of i1 ... j_n_so
# Returns a tuple of (i, j)
 
def Get_bin(*args):
    bin_i, bin_j = '', ''
    count = 0
    for arg in args:
        if count < n_so: 
            bin_i += str(arg)
            count += 1
            continue
        else:
            bin_j += str(arg)
            count += 1
    i, j =  int(bin_i, 2), int(bin_j, 2)
    if i > n_so: i = 0 # out of range of orbitals and so = 0
    if j > n_so: j = 0 # out of range of orbitals and so = 0
# i, j = 0 is nonsensical, since i, j = 1 ... n_so 
    return (i, j)

####################################

def First_term(*args):
    i, j = Get_bin(*args) # *args are al the qubits from init_orbital_list
    if (i == 0) or (j == 0): return False # Since i, j = 0 is nonsensical
    transformed_state = Create(Anhilate(state_i, j), i) # transformed_state = a†_i a_j |Ψ>

# To enforce a†_i a_j |Ψ> = C * |Ψ>,  uncomment
#    if transformed_state == input_state:
#        return True
#    else: return False
 
    
# All valid transofrmations will return true   
    if transformed_state == False:
        return False
    else:
        print("Final transformed states are:%s" % transformed_state)
        return True  
    
####################################
####################################
        
global state_i # Define the initial state |Ψ_i> as a list
state_i = [1, 0, 0, 0, 0, 0] # Each index corresponds to a location for electron
                 # 0 -> 1S1, 1-> 1S2 3 -> 2S1, 4 -> 2S2, 5-> 2Px1 ... 
n_so = len(state_i) # Number of location for electron

var_list = init_orbital_list(n_so) # List of qubits

# Creates a BQM for a constraint problem
# For details, look at: 
#https://docs.ocean.dwavesys.com/en/latest/examples/scheduling.html

csp = dwavebinarycsp.ConstraintSatisfactionProblem('BINARY')
csp.add_constraint(First_term, var_list)
bqm = dwavebinarycsp.stitch(csp, max_graph_size=20) # Change max_graph_size > 2*n_so


m = bqm.to_networkx_graph() # Convert BQM to a graph, to see if patterns exist 
networkx.draw(m) # Draw graph of BQM
print("DONE")


#Solvers

#resp = dimod.ExactSolver().sample(bqm)
#
#for sample, energy in resp.data(['sample', 'energy']):
#    print(sample, csp.check(sample))#, energy)
