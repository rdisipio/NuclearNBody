#!/usr/bin/env python3

import numpy as np

from pauli import *

def to_dirac( ψ ):
    ket = "|"
    n = len(ψ) // 2
    #print(ψ.reshape( n//2, 2))

    s = 1
    for p in ψ.reshape( n, 2):
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
        ψ = np.kron( ψ, zero )

    n_basis = len(ψ)
    basis = [ None for _ in range(n_basis) ]
    for i in range(n_basis):
        basis[i] = np.zeros(n_basis)
        basis[i][i] = 1

    a_dagger = [ None for _ in range(n) ]

    # init a†1, a†2, ...,a†n
    for i in range(n):

        n_I = n-i-1
        n_Z = i

        a_dagger[i] = P
        for _ in range(n_I):
            a_dagger[i] = np.kron( a_dagger[i], I )

        for _ in range(n_Z):
            a_dagger[i] = np.kron( Z, a_dagger[i] )
        
        #print(a_dagger[i])

    return a_dagger, ψ, basis

if __name__ == "__main__":

    n = 3
    r = 2

    a_dagger, ψ, basis = prepare(n)

    print("ψ =", ψ, "=", to_dirac(ψ))

    for i in range(n):
        print("a†%i =" % (i+1))
        print(a_dagger[i])

    print("ψ =", ψ, "=", to_dirac(ψ))

    for i in range(n):
        s = np.dot(a_dagger[i],ψ)
        print("a†%i |ψ> = "%(i+1), s, "=", to_dirac(s) )

    x21 = np.dot( a_dagger[0], np.dot(a_dagger[1],ψ))
    x12 = np.dot( a_dagger[1], np.dot(a_dagger[0],ψ))

    print("a†2 a†1 |ψ> =", x21, "=", to_dirac(x21) )
    print("a†1 a†2 |ψ> =", x12, "=", to_dirac(x12) )