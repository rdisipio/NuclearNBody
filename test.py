#!/usr/bin/env python3

import numpy as np
from scipy.linalg import block_diag

from pauli import *

def to_string( ψ ):
    ket = "|"
    n = len(ψ)
    #print(ψ.reshape( n//2, 2))

    s = 1
    for p in ψ.reshape( n//2, 2):
        if np.array_equal(p, zero):
            ket += "0"
        elif np.array_equal(p, one): 
            ket += "1"
        elif np.array_equal(p, -zero):
            ket += "0"
            s *= -1
        elif np.array_equal(p, -one): 
            ket += "1"
            s *= -1
    ket += ">"
    if s < 0:
        ket = "-" + ket
    return ket

def prepare(n):
    ψ = np.array( zero )
    for i in range(n-1):
        ψ = np.hstack( [ψ, zero ] )

    a_dagger = [ np.zeros([2*n,2*n]) for i in range(n) ]

    # init a†1, a†2, a†n
    for i in range(n):

        n_I = n-i-1
        n_Z = i

        j = 0
        for _ in range(n_Z):
            a_dagger[i][2*j:2*j+2,2*j:2*j+2] = Z
            j += 1
        
        a_dagger[i][2*j:2*j+2,2*j:2*j+2] = P
        j += 1

        for _ in range(n_I):
            a_dagger[i][2*j:2*j+2,2*j:2*j+2] = I
            j += 1

        #print(a_dagger[i])

        #a_dagger[i] = block_diag(a_dagger[i], P)
        #for _ in range(n_I):
        #    a_dagger[i] = block_diag(a_dagger[i], I)

    return a_dagger, ψ

if __name__ == "__main__":

    n = 3
    r = 2

    a_dagger, ψ = prepare(n)

    print("ψ =", ψ, "=", to_string(ψ))

    for i in range(n):
        print("a†%i =" % (i+1))
        print(a_dagger[i])

    print("ψ =", ψ, "=", to_string(ψ))

    for i in range(n):
        s = np.dot(a_dagger[i],ψ)
        print("a†%i |ψ> = "%(i+1), s, "=", to_string(s) )

    x21 = np.dot( a_dagger[0], np.dot(a_dagger[1],ψ))
    x12 = np.dot( a_dagger[1], np.dot(a_dagger[0],ψ))

    print("a†2 a†1 |ψ> =", x21, "=", to_string(x21) )
    print("a†1 a†2 |ψ> =", x12, "=", to_string(x12) )