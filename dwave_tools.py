import dwave_networkx as dnx
import minorminer
import numpy as np
import sys
from dwave.system.composites import FixedEmbeddingComposite
from dwave.system.samplers import DWaveSampler

#########################################


def anneal_sched_custom(id=0):
    if id == 0:
        return (
            (0.0, 0.0),  # Start everything at 0
            (10.0, 0.50),
            (110, 0.50),
            (120.0, 1.0))
    elif id == 1:
        return (
            (0.0, 0.0),  # Start everything at 0
            (1.0, 0.25),  # Quickly ramp up to 0.25 anneal at 1us
            (19.0,
             0.75),  # Hold the anneal at 0.25 until 19us, then go up to 0.75
            (20.0, 1.0)  # End the full anneal at 20us
        )
    elif id == 2:
        return (
            (0.0, 0.0),  # Start everything at 0
            (20.0, 0.50),
            (280, 0.50),
            (300.0, 1.0))
    elif id == 3:
        return ((0.0, 0.0), (10.0, 0.5), (12.0, 1.0))  # with quench
    else:
        return None


def max_chain_length(embedding: dict) -> int:
    max_ = 0
    for _, chain in embedding.items():
        if len(chain) > max_:
            max_ = len(chain)
    return max_


#########################################


def get_embedding_with_short_chain(S: dict,
                                   tries: int = 5,
                                   processor: list = None,
                                   verbose=True) -> dict:
    '''Try a few probabilistic embeddings and return the one with the shortest
    chain length
    :param S: Source couplings
    :param tries: Number of probabilistic embeddings
    :param verbose: Whether to print out diagnostic information
    :return: Returns the minor embedding
    '''
    if processor is None:
        # The hardware topology: 16 by 16 pieces of K_4,4 unit cells
        processor = dnx.chimera_graph(16, 16, 4).edges()
    # Try a few embeddings
    embedding = None
    best_chain_length = sys.maxsize
    source = list(S.keys())
    for itry in range(tries):
        #        print(".",  end=' ')
        try:
            emb = minorminer.find_embedding(source, processor, verbose=0)
            chain_length = max_chain_length(emb)
            if (chain_length > 0) and (chain_length < best_chain_length):
                embedding = emb
                best_chain_length = chain_length
                if verbose:
                    print("DEBUG: found embedding at attempt", itry,
                          "max length:", chain_length)
        except:
            pass
    if embedding == None:
        print("WARNING: could not find minor embedding")
        return embedding

    if verbose:
        max_length = max_chain_length(embedding)
        print("INFO: max chain length in best embedding:", max_length)
        print("INFO: best embedding:")
        print(embedding)

    if best_chain_length == sys.maxsize:
        raise Exception("Cannot find embedding")
    return embedding


#########################################


def get_energy(bqm, sample):
    # see https://docs.ocean.dwavesys.com/projects/dimod/en/latest/_modules/dimod/reference/samplers/exact_solver.html
    M = bqm.binary.to_numpy_matrix()
    off = bqm.binary.offset
    E = sample.dot(M).dot(sample.transpose())

    return float(E) + off


#########################################


# define a qbsolv-like workflow
def merge_substates(_, substates):
    a, b = substates
    return a.updated(
        subsamples=hybrid.hstack_samplesets(a.subsamples, b.subsamples))


#########################################


def make_reverse_anneal_schedule(s_target=0.0,
                                 hold_time=10.0,
                                 ramp_back_slope=0.2,
                                 ramp_up_time=0.0201,
                                 ramp_up_slope=None):
    """Build annealing waveform pattern for reverse anneal feature.
    Waveform starts and ends at s=1.0, descending to a constant value
    s_target in between, following a linear ramp.
      s_target:   s-parameter to descend to (between 0 and 1)
      hold_time:  amount of time (in us) to spend at s_target (must be >= 2.0us)
      ramp_slope: slope of transition region, in units 1/us
    """
    # validate parameters
    if s_target < 0.0 or s_target > 1.0:
        raise ValueError("s_target must be between 0 and 1")
    if hold_time < 0.0:
        raise ValueError("hold_time must be >= 0")
    if ramp_back_slope > 0.2:
        raise ValueError("ramp_back_slope must be < 0.2")

    ramp_time = (1.0 - s_target) / ramp_back_slope

    initial_s = 1.0
    pattern = [[0.0, initial_s]]

    # don't add new points if s_target == 1.0
    if s_target < 1.0:
        pattern.append([round(ramp_time, 4), round(s_target, 4)])
        if hold_time != 0:
            pattern.append(
                [round(ramp_time + hold_time, 4),
                 round(s_target, 4)])

    # add last point
    if ramp_up_slope is not None:
        ramp_up_time = (1.0 - s_target) / ramp_up_slope
        pattern.append(
            [round(ramp_time + hold_time + ramp_up_time, 4),
             round(1.0, 4)])
    else:
        pattern.append(
            [round(ramp_time + hold_time + ramp_up_time, 4),
             round(1.0, 4)])

    return pattern


#########################################


def qubo_quadratic_terms_from_np_array(Q):
    J = {}
    nx = Q.shape[0]
    ny = Q.shape[1]
    for a in range(nx):
        for b in range(a + 1, ny):
            idx = (a, b)
            J[idx] = Q[a][b]
    return J
