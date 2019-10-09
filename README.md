# NuclearNBody
Nuclear N-body calculation using DWave quantum annealer.

Requires: 
```
dimod
dwave-system
dwave-neal
dwave_networkx
numpy
pyqubo
minorminer
```

To run test example:

```
./nbody.py -s n_spin_orbitals -p n_particles -r num_reads -m max_evals
```

Change ```h``` and ```V``` in ```make_hamiltonian()``` according to the problem.

