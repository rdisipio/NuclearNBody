# -*- coding: utf-8 -*-

#import neal
import dwavebinarycsp
import operator as op
import pyqubo
import dimod
import networkx

n_so = 3

def init_state(n_p, n_so):
    iter_np = 0
    init_state_dict = {}
    
    for i in range(n_so):
        if iter_np < n_p:
            init_state_dict.update({i+1: 1 })
        else:
            init_state_dict.update({i+1: 0 })
        iter_np += 1
    return init_state_dict

def init_orbital_list(n_so):
    idx = 0
    so_list = []
    idx2str = {0: 'a', 1: 'i', 2: 'j'}
    while idx != 3:
        for i in range(n_so):
            so_list += [idx2str[idx] + str(i+1)]
        idx += 1
    return so_list

def Get_bin(*args):
    state, bin_i, bin_j = '', '', ''
    count = 0
    for arg in args:
        if count < n_so: 
            state += str(arg)
            count += 1
            continue
        elif count < 2*n_so: 
            bin_i += str(arg)
            count += 1
            continue
        else:
            bin_j += str(arg)
            count += 1
    i, j =  int(bin_i, 2), int(bin_j, 2)
    if i > n_so: i = 0
    if j > n_so: j = 0
    return (list(state), i, j)



def Anhilate(state, anhil):
    if state[anhil - 1] == 0:
        return False
    else:
        state[anhil - 1] = 0
        return state 

def Create(state, create):
    if state == False:
        return False
    if state[create - 1] == 1: return False
    else:
        state[create - 1] = 1
        return state
    

def Operate(initial_state, create, anhil):
    transformed_state = Create(Anhilate(initial_state, anhil), create)
    if transformed_state == initial_state:
        return True
    else:
        return False

def First_term(*args):
    state, i, j = Get_bin(*args)
    if (i == 0) or (j == 0): return False
    return Operate(state, i, j)


var_list = init_orbital_list(n_so)

csp = dwavebinarycsp.ConstraintSatisfactionProblem('BINARY')


csp.add_constraint(First_term, var_list)


bqm = dwavebinarycsp.stitch(csp, max_graph_size=40)
m = bqm.to_numpy_matrix(var_list)
print(m)
#resp = dimod.ExactSolver().sample(bqm)
#
#for sample, energy in resp.data(['sample', 'energy']):
#    print(sample, csp.check(sample), energy)
