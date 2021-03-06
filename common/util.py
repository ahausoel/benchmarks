import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from pytriqs.gf import *
from pytriqs.operators import c, c_dag, n, dagger
from pytriqs.utility import mpi

# Get a list of all annihilation operators from a many-body operators
def get_fundamental_operators(op):
    idx_lst = []
    for term, val in op:
        for has_dagger, (bl, orb) in term:
            if not idx_lst.count([bl, orb]):
                idx_lst.append([bl,orb])
    return [c(bl, orb) for bl, orb in idx_lst]

# print on master node
def mpi_print(arg):
    if mpi.is_master_node():
        print arg
