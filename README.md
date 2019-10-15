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
networkx
```

To run test example:

```
./nbody.py -s n_spin_orbitals -f n_fermions -r num_reads -m max_evals
```

Change ```h``` and ```V``` in ```make_hamiltonian()``` according to the problem.


```Reverse_anealing.py``` creates a BQM and a graph to represent it. This code requires an initial state and then outputs a BQM where all possible single electron transformations have the lowest energy. Read the comments for more details. To run, simply change varriable ```state_i``` with desired initial state and then run.

```All_permutations_algorithm.py``` has extra qubits to represent the various initial states. This has a longer computation time but it allows for all permutations of the initial state.
