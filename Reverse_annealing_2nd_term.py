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


def Get_bin_2nd_term(*args):
    bin_i, bin_j, bin_k, bin_l = '', '', '', ''
    count = 0
    for arg in args:
        if count < n_so: 
            bin_i += str(arg)
            count += 1
            continue
        if count < 2*n_so: 
            bin_j += str(arg)
            count += 1
            continue
        if count < 3*n_so: 
            bin_k += str(arg)
            count += 1
            continue
        else:
            bin_l += str(arg)
    i, j, k, l =  int(bin_i, 2), int(bin_j, 2), int(bin_k, 2), int(bin_l, 2)
    if i > n_so: i = 0 # out of range of orbitals and so = 0
    if j > n_so: j = 0 # out of range of orbitals and so = 0
    if k > n_so: k = 0
    if l > n_so: l = 0
# i, j = 0 is nonsensical, since i, j = 1 ... n_so 
    return (i, j, k, l)
####################################

def First_term(*args):
    i, j = Get_bin(*args) # *args are al the qubits from init_orbital_list
    if (i == 0) or (j == 0): return True # Since i, j = 0 is nonsensical
    transformed_state = Create(Anhilate(state_i, j), i) # transformed_state = a†_i a_j |Ψ>

# To enforce a†_i a_j |Ψ> = C * |Ψ>,  uncomment
#    if transformed_state == input_state:
#        return True
#    else: return False
 
    
# All valid transofrmations will return true   
    if transformed_state == False:
        return False
    else:
        print("Final transformations in binary are:%s with transformation %s" 
              % ((format(i, 'b'), format(j, 'b')), (i, j)))
        return True  
 
    
def Second_term(*args):
    i, j, k, l = Get_bin_2nd_term(*args)
    if (i == 0) or (j == 0) or (k == 0) or (l == 0): return False
    anhil_state = Anhilate(Anhilate(state_i, l), k)
    transformed_state = Create(Create(anhil_state, j), i)
    
    if transformed_state == False:
        return False
    else:
        print("Final transformations are: %s, %s, %s, %s" % (i, j, k, l))
        return True  

####################################
####################################
        
global state_i # Define the initial state |Ψ_i> as a list
state_i = [1,] # Each index corresponds to a location for electron
                 # 0 -> 1S1, 1-> 1S2 3 -> 2S1, 4 -> 2S2, 5-> 2Px1 ... 
n_so = len(state_i) # Number of location for electron

var_list = init_orbital_list_2nd(n_so) # List of qubits

# Creates a BQM for a constraint problem
# For details, look at: 
#https://docs.ocean.dwavesys.com/en/latest/examples/scheduling.html

csp = dwavebinarycsp.ConstraintSatisfactionProblem('BINARY')
csp.add_constraint(Second_term, var_list)
#print(datetime.now() - startTime)
bqm = dwavebinarycsp.stitch(csp, max_graph_size=20, min_classical_gap = .01) # Change max_graph_size > 2*n_so
print(datetime.now() - startTime)

print(bqm)
m = bqm.to_networkx_graph() # Convert BQM to a graph, to see if patterns exist 

edge_color = []
edge_list = []
node_color = []
for q in m.node:
    if m._node[q]['bias'] > 0:
        node_color += ['y']
        
    if m._node[q]['bias'] == 0:
        node_color += ['k']
        
    if m._node[q]['bias'] < 0:
        node_color += ['r']

for q in m.edges:
    edge_list += [q]
    if m[q[0]][q[1]]['bias'] > 0:
        edge_color += ['y']
        
    if m[q[0]][q[1]]['bias'] == 0:
        edge_color += ['k']
        
    if m[q[0]][q[1]]['bias'] < 0:
        edge_color += ['r']
matrix = bqm.to_numpy_matrix(var_list)
matrix_squared = np.matmul(matrix, matrix)


simplified_matrix = np.zeros(np.shape(matrix))

for i in range(np.shape(matrix)[0]):
    for j in range(np.shape(matrix)[1]):
        if matrix[i][j] > 0:
            simplified_matrix[i][j] = 1
            
        if matrix[i][j] < 0:
            simplified_matrix[i][j] = -1
simplified_matrix_squared = np.matmul(simplified_matrix, simplified_matrix)
eig, eigvect = np.linalg.eig(matrix)
det = np.linalg.det(matrix)

#plt.figure()
#pos = networkx.spring_layout(m)
#networkx.draw_networkx_nodes(m, pos, node_color = node_color)
#networkx.draw_networkx_edges(m, pos, edgelist = edge_list,
#                             edge_color = edge_color)
#networkx.draw_networkx_labels(m, pos, font_color = 'w')
#plt.text(1.15, 0.8, " Black: = 0 \n Red: < 0 \n Yellow: > 0")
#networkx.draw_networkx_edge_labels(m, pos)

print("DONE")

#Solvers
#
#resp = dimod.ExactSolver().sample(bqm)
#
#for sample, energy in resp.data(['sample', 'energy']):
#    print(sample, energy)
