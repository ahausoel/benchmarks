import sys, os
sys.path.append(os.getcwd() + "/..")
sys.path.append(os.getcwd() + "/../../common")
from model import *

from pytriqs.archive import HDFArchive
from pytriqs.utility import mpi
from pytriqs.cthyb import Solver, version

# --------- Construct the CTHYB solver ----------
constr_params = {
        'beta' : beta,
        'gf_struct' : dict(gf_struct),
        'n_iw' : n_iw,
        'n_tau' : 100001
        }
S = Solver(**constr_params)

# --------- Initialize G0_iw ----------
S.G0_iw << G0_iw

# --------- Solve! ----------
solve_params = {
        'h_int' : h_int,
        'n_warmup_cycles' : 1000,
        'n_cycles' : 10000000,
        'length_cycle' : 200
        }
S.solve(**solve_params)

# -------- Save in archive ---------
if mpi.is_master_node():
    with HDFArchive("../results/cthyb.h5",'w') as results:
        results["G"] = S.G_iw

        import inspect
        import __main__
        results.create_group("Solver_Info")
        info_grp = results["Solver_Info"]
        info_grp["solver_name"] = "triqs_cthyb"
        info_grp["version"] = version.version
        info_grp["git_hash"] = version.cthyb_hash
        info_grp["triqs_git_hash"] = version.triqs_hash
        info_grp["script"] = inspect.getsource(__main__)
        info_grp["constr_params"] = constr_params
        info_grp["solve_params"] = solve_params