import sys, os
sys.path.append(os.getcwd() + "/..")
sys.path.append(os.getcwd() + "/../../common")
from model import *

from pytriqs.gf import Fourier, Idx
from pytriqs.archive import HDFArchive
from math import exp

# -------- States ---------
states = ['0', 'up', 'dn', 'updn']

# -------- Eigenenergies ---------
E = {}
E['0'] = 0.
E['up'] = -mu -h
E['dn'] = -mu +h
E['updn'] = U - 2.*mu

# -------- Partition function ---------
Z = sum( [exp(-beta*E[i]) for i in states] )

# -------- Density-Matrix ---------
rho = {}
rho['0'] = exp(-beta*E['0']) / Z
rho['up'] = exp(-beta*E['up']) / Z
rho['dn'] = exp(-beta*E['dn']) / Z
rho['updn'] = exp(-beta*E['updn']) / Z

# -------- Occupation ---------
n = {}
n['up'] = rho['up'] + rho['updn']
n['dn'] = rho['dn'] + rho['updn']


# ================ Green Function G_iw ================

# Exact expressions for Green function according to
# Eq. (B.3) and (B.4) of S. Pairault, et al., EPJ B 16, 85-105 (2000)

G_iw = BlockGf_from_struct(mesh=iw_mesh, struct=gf_struct)

G_iw['up'] << ( 1 - n['dn'] ) * inverse( iOmega_n + mu + h ) + \
                    n['dn']   * inverse( iOmega_n + mu + h  - U )
G_iw['dn'] << ( 1 - n['up'] ) * inverse( iOmega_n + mu - h ) + \
                    n['up']   * inverse( iOmega_n + mu - h  - U )

# Green function in Lambda form
G = { 'up': lambda iw: ( 1 - n['dn'] ) / ( iw + mu + h ) + n['dn'] / ( iw + mu + h - U ),
      'dn': lambda iw: ( 1 - n['up'] ) / ( iw + mu - h ) + n['up'] / ( iw + mu - h - U ) }

# -------- Save in archive ---------
res = HDFArchive("../results/exact.h5",'w')
res["G"] = G_iw


# ================ chi3pp_tau ==================

# Exact expressions for Equal-time two-particle Green function
# chi3_pp(tau,tau') =   <T cdag_up(tau^{+}) c_up(0^{+}) cdag_dn(tau') c_dn(0)>
#                   = - <T cdag_up(tau^{+}) cdag_dn(tau') c_up(0^{+}) c_dn(0)>

def chi3pp(block, tau, taup):

    # beta anti-periodicity
    if tau > beta:
        return -chi3pp(block, tau - beta, taup)
    elif tau < 0:
        return -chi3pp(block, tau + beta, taup)
    elif taup > beta:
        return -chi3pp(block, tau, taup - beta)
    elif taup < 0:
        return -chi3pp(block, tau, taup + beta)

    if tau >= taup:
        if block == ('up','dn'):
            # chi3pp_updn = - rho_updn * <updn| cdag_up(tau^{+}) cdag_dn(tau') c_up(0^{+}) c_dn(0) |updn>
            # Take care of state ordering inside |updn>, i.e. c_dn |updn> = -|up>
            return rho['updn'] * exp(E['updn']*tau + E['dn']*(taup - tau))
        elif block == ('dn','up'):
            # chi3pp_dnup = - rho_updn * <updn| cdag_dn(tau^{+}) cdag_up(tau') c_dn(0^{+}) c_up(0) |updn>
            # Take care of state ordering inside |updn>, i.e. c_dn |updn> = -|up>
            return rho['updn'] * exp(E['updn']*tau + E['up']*(taup - tau))
        elif block == ('up','up'):
            # chi3pp_upup = 0, There is no state for which <m| cdag_up(tau^{+}) cdag_up(tau') c_up(0^{+}) c_up(0) |m> is finite
            return 0.
        elif block == ('dn','dn'):
            # chi3pp_dndn = 0, There is no state for which <m| cdag_dn(tau^{+}) cdag_dn(tau') c_dn(0^{+}) c_dn(0) |m> is finite
            return 0.
        else:
            print "Invalid Block"
            raise
    else: # ===  tau < taup , Take care of overall Minus Sign
        if block == ('up','dn'):
            # chi3pp_updn = rho_updn * <updn| cdag_dn(tau') cdag_up(tau) c_up(0^{+}) c_dn(0) |updn>
            return rho['updn'] * exp(E['updn']*taup + E['up']*(tau - taup))
        elif block == ('dn','up'):
            # chi3pp_dnup = rho_updn * <updn| cdag_up(tau') cdag_dn(tau) c_dn(0^{+}) c_up(0) |updn>
            return rho['updn'] * exp(E['updn']*taup + E['dn']*(tau - taup))
        elif block == ('up','up'):
            # chi3pp_upup = 0, There is no state for which <m| cdag_up(tau') cdag_up(tau) c_up(0^{+}) c_up(0) |m> is finite
            return 0.
        elif block == ('dn','dn'):
            # chi3pp_dndn = 0, There is no state for which <m| cdag_dn(tau') cdag_dn(tau) c_dn(0^{+}) c_dn(0) |m> is finite
            return 0.
        else:
            print "Invalid Block"
            raise


# -------- Initialize chi3pp_tau Container ---------

tau_mesh = MeshImTime(beta, 'Fermion', 2*n_iw+1)
chi3pp_tau = Block2Gf_from_struct(mesh=MeshProduct(tau_mesh, tau_mesh), struct=gf_struct)

for tau, taup in chi3pp_tau[('up','dn')].mesh:
    for block in [('up','up'), ('up','dn'), ('dn','dn'), ('dn','up')]:
        chi3pp_tau[block][tau, taup] = chi3pp(block, tau.value, taup.value)

# -------- Save in archive ---------
res["chi3pp_tau"] = chi3pp_tau


# ================ chi3ph_tau ==================

# Exact expressions for Equal-time two-particle Green function
# chi3_ph(tau,tau') = <T cdag(tau^{+}) c(tau') n(0)>

def chi3ph(block, tau, taup):

    # beta anti-periodicity
    if tau > beta:
        return -chi3ph(block, tau - beta, taup)
    elif tau < 0:
        return -chi3ph(block, tau + beta, taup)
    elif taup > beta:
        return -chi3ph(block, tau, taup - beta)
    elif taup < 0:
        return -chi3ph(block, tau, taup + beta)

    if tau >= taup:
        if block == ('up','dn'):
            # chi3ph_updn = rho_updn * <updn| cdag_up(tau^{+}) c_up(tau') n_dn(0) |updn>
            #             = rho_updn * <updn| cdag_up(tau^{+}) |dn><dn| c_up(tau') |updn><updn| n_dn(0) |updn>
            return rho['updn'] * exp(E['dn']*(taup - tau) + E['updn']*(tau - taup))
        elif block == ('dn','up'):
            # chi3ph_dnup = rho_updn * <updn| cdag_dn(tau^{+}) c_dn(tau') n_up(0) |updn>
            #             = rho_updn * <updn| cdag_up(tau^{+}) |up><up| c_up(tau') |updn><updn| n_dn(0) |updn>
            return rho['updn'] * exp(E['up']*(taup - tau) + E['updn']*(tau - taup))
        elif block == ('up','up'):
            # chi3ph_upup =   rho_up   * <up|   cdag_up(tau^{+}) c_up(tau') n_up(0) |up>
            #               + rho_updn * <updn| cdag_up(tau^{+}) c_up(tau') n_up(0) |updn>
            return rho['up'] * exp(E['up']*(tau - taup)) + rho['updn'] * exp(E['dn']*(taup - tau) + E['updn']*(tau - taup))
        elif block == ('dn','dn'):
            # chi3ph_dndn =   rho_dn   * <dn|   cdag_dn(tau^{+}) c_dn(tau') n_dn(0) |dn>
            #               + rho_updn * <updn| cdag_dn(tau^{+}) c_dn(tau') n_dn(0) |updn>
            return rho['dn'] * exp(E['dn']*(tau - taup)) + rho['updn'] * exp(E['up']*(taup - tau) + E['updn']*(tau - taup))
        else:
            print "Invalid Block"
            raise
    else: # ===  tau < taup , Take care of overall Minus Sign
        if block == ('up','dn'):
            # chi3ph_updn = rho_dn * <dn| c_up(tau') cdag_up(tau) n_dn(0) |dn>
            #             = rho_dn * <dn| c_up(tau') |updn><updn| cdag_up(tau) |dn><dn| n_dn(0) |dn>
            return - rho['dn'] * exp(E['updn']*(tau - taup) + E['dn']*(taup - tau))
        elif block == ('dn','up'):
            # chi3ph_updn = rho_up * <up| c_dn(tau') cdag_dn(tau) n_up(0) |up>
            #             = rho_up * <up| c_dn(tau') |updn><updn| cdag_dn(tau) |up><up| n_up(0) |up>
            return - rho['up'] * exp(E['updn']*(tau - taup) + E['up']*(taup - tau))
        elif block == ('up','up'):
            # chi3ph_upup = 0, There is no state for which <m| c_up(tau') cdag_up(tau) n_up(0) |m> is finite
            return 0.
        elif block == ('dn','dn'):
            # chi3ph_dndn = 0, There is no state for which <m| c_dn(tau') cdag_dn(tau) n_dn(0) |m> is finite
            return 0.
        else:
            print "Invalid Block"
            raise


# -------- Initialize chi3ph_tau Container ---------

tau_mesh = MeshImTime(beta, 'Fermion', 2*n_iw+1)
chi3ph_tau = Block2Gf_from_struct(mesh=MeshProduct(tau_mesh, tau_mesh), struct=gf_struct)

for tau, taup in chi3ph_tau[('up','dn')].mesh:
    for block in [('up','up'), ('up','dn'), ('dn','dn'), ('dn','up')]:
        chi3ph_tau[block][tau, taup] = chi3ph(block, tau.value, taup.value)

# -------- Save in archive ---------
res["chi3ph_tau"] = chi3ph_tau


# ================ chi2pp_tau ==================

# Exact expressions for Equal-time two-particle Green function
# chi2_pp(tau) =   <T cdag_up(tau^{+}) c_up(0^{+}) cdag_dn(tau) c_dn(0)>
#              = - <T cdag_up(tau^{+}) cdag_dn(tau) c_up(0^{+}) c_dn(0)>

# -------- Initialize chi2ph_tau Container ---------

chi2pp_tau = Block2Gf_from_struct(MeshImTime(beta, 'Boson', 2*n_iw+1), struct=gf_struct)

for tau in chi2pp_tau[('up','dn')].mesh:
    for block in [('up','up'), ('up','dn'), ('dn','dn'), ('dn','up')]:
        chi2pp_tau[block][tau] = chi3pp(block, tau.value, tau.value)

# -------- Save in archive ---------
res["chi2pp_tau"] = chi2pp_tau


# ================ chi2ph_tau ==================

# Exact expressions for Equal-time two-particle Green function
# chi2_ph(tau) =   <T cdag_up(tau^{+}) c_up(tau) cdag_dn(0^{+}) c_dn(0)>

# -------- Initialize chi2ph_tau Container ---------

chi2ph_tau = Block2Gf_from_struct(MeshImTime(beta, 'Boson', 2*n_iw+1), struct=gf_struct)

for tau in chi2ph_tau[('up','dn')].mesh:
    for block in [('up','up'), ('up','dn'), ('dn','dn'), ('dn','up')]:
        chi2ph_tau[block][tau] = chi3ph(block, tau.value, tau.value)

# -------- Save in archive ---------
res["chi2ph_tau"] = chi2ph_tau


# ================ chi3pp_iw ==================

# Exact expressions for Equal-time two-particle Green function
# chi3_pp(iw,iw') = - <cdag(iw) cdag(iw') c(0^{+}) c(0)>
# Compare https://hal.archives-ouvertes.fr/tel-01247625v1 Page 217

def f(i, j, k):
    if rho[k] == rho[i]:
        return lambda iw, iwp: 1. / (iwp + E[j] - E[k]) * ( rho[i] + rho[j] ) / (iw + E[i] - E[j]) + beta * kronecker(iw, -iwp) * rho[i] / (iwp + E[j] - E[i])
    else:
        return lambda iw, iwp: 1. / (iwp + E[j] - E[k]) * (  ( rho[k] - rho[i] ) / (iw + iwp + E[i] - E[k]) + ( rho[i] + rho[j] ) / (iw + E[i] - E[j]) )


def chi3pp(block, iw, iwp):

    if block == ('up','dn'):
        # chi3pp = - <i|cdag_up|j> <j|cdag_dn|k> <k|c_up c_dn|i> + <i|cdag_dn|j> <j|cdag_up|k> <k|c_up c_dn|i>
        # Take care of state ordering inside |updn>, i.e. c_dn |updn> = -|up>
        return f("updn","dn","0")(-iw,-iwp) + f("updn","up","0")(-iwp,-iw)
    elif block == ('dn','up'):
        # chi3pp = - <i|cdag_dn|j> <j|cdag_up|k> <k|c_dn c_up|i> + <i|cdag_up|j> <j|cdag_dn|k> <k|c_dn c_up|i>
        # Take care of state ordering inside |updn>, i.e. c_dn |updn> = -|up>
        return f("updn","up","0")(-iw,-iwp) + f("updn","dn","0")(-iwp,-iw)
    elif block == ('up','up'):
        # chi3pp_upup = 0, There is no state for which <m| cdag_up(tau^{+}) cdag_up(tau') c_up(0^{+}) c_up(0) |m> or
        # or  <m| cdag_up(tau') cdag_up(tau) c_up(0^{+}) c_up(0) |m> is finite
        return 0.
    elif block == ('dn','dn'):
        # chi3pp_dndn = 0, There is no state for which <m| cdag_dn(tau^{+}) cdag_dn(tau') c_dn(0^{+}) c_dn(0) |m> or
        # <m| cdag_dn(tau') cdag_dn(tau) c_dn(0^{+}) c_dn(0) |m> is finite
        return 0.
    else:
        print "Invalid Block"
        raise

def chi3pp_conn(block, iw, iwp):
    bl1, bl2 = block
    dens = G_iw[bl2].density()
    return chi3pp(block, iw, iwp) - G[bl1](iw) * G[bl2](iwp) + kronecker(bl1, bl2) * G[bl1](iw) * G[bl2](iwp)

# -------- Initialize chi3pp_iw Container ---------

iW_mesh = MeshImFreq(beta, 'Boson', n_iw)
chi3pp_iw = Block2Gf_from_struct(mesh=MeshProduct(iw_mesh, iW_mesh), struct=gf_struct)
chi3pp_conn_iw = chi3pp_iw.copy()

for iw, iW in chi3pp_iw[('up','dn')].mesh:
    for block in [('up','up'), ('up','dn'), ('dn','dn'), ('dn','up')]:
        iw1 = iw.value
        iw3 = iW.value - iw.value
        chi3pp_iw[block][iw, iW] = chi3pp(block, iw1, iw3)
        chi3pp_conn_iw[block][iw, iW] = chi3pp_conn(block, iw1, iw3)

# -------- Save in archive ---------
res["chi3pp_iw"] = chi3pp_iw
#res["chi3pp_conn_iw"] = chi3pp_conn_iw


# ================ chi3ph_iw ==================

# Exact expressions for Equal-time two-particle Green function
# chi3_ph(iw,iw') = <cdag(iw) c(iw') n(0)>
# Compare https://hal.archives-ouvertes.fr/tel-01247625v1 Page 217

def chi3ph(block, iw, iwp):

    if block == ('up','dn'):
        # chi3ph_updn = <i|cdag_up|j> <j|c_up|k> <k|n_dn|i> - <i|c_up|j> <j|cdag_up|k> <k|n_dn|i>
        return f("updn", "dn", "updn")(-iw, iwp) - f("dn", "updn", "dn")(iw, -iwp)
    elif block == ('dn','up'):
        # chi3ph_dnup = <i|cdag_dn|j> <j|c_dn|k> <k|n_up|i> - <i|c_dn|j> <j|cdag_dn|k> <k|n_up|i>
        return f("updn", "up", "updn")(-iw, iwp) - f("up", "updn", "up")(iw, -iwp)
    elif block == ('up','up'):
        # chi3ph_dnup = <i|cdag_up|j> <j|c_up|k> <k|n_up|i>
        return f("up", "0" ,"up")(-iw, iwp) + f("updn", "dn" ,"updn")(-iw, iwp)
    elif block == ('dn','dn'):
        # chi3ph_dnup = <i|cdag_dn|j> <j|c_dn|k> <k|n_dn|i>
        return f("dn", "0" ,"dn")(-iw, iwp) + f("updn", "up" ,"updn")(-iw, iwp)
    else:
        print "Invalid Block"
        raise

def chi3ph_conn(block, iw, iwp):
    bl1, bl2 = block
    dens = G_iw[bl2].density()
    return chi3ph(block, iw, iwp) - beta * kronecker(iw,iwp) * G[bl1](iw) * dens + kronecker(bl1, bl2) * G[bl1](iw) * G[bl2](iwp)

# -------- Initialize chi3ph_iw Container ---------

chi3ph_iw = Block2Gf_from_struct(mesh=MeshProduct(iw_mesh, iW_mesh), struct=gf_struct)
chi3ph_conn_iw = chi3ph_iw.copy()

for iw, iW in chi3ph_iw[('up','dn')].mesh:
    for block in [('up','up'), ('up','dn'), ('dn','dn'), ('dn','up')]:
        iw1 = iw.value
        iw2 = iW.value + iw.value
        chi3ph_iw[block][iw, iW] = chi3ph(block, iw1, iw2)
        chi3ph_conn_iw[block][iw, iW] = chi3ph_conn(block, iw1, iw2)

# -------- Save in archive ---------
res["chi3ph_iw"] = chi3ph_iw
#res["chi3ph_conn_iw"] = chi3ph_conn_iw


# ================ chi2pp_iw ==================

# Exact expressions for Equal-time two-particle Green function
# chi2_pp(iW) = - < (cdag cdag)(iW) c(0^{+}) c(0)>
# Compare https://hal.archives-ouvertes.fr/tel-01247625v1 Page 217

def f(i, j):
    if E[i] == E[j]:
        return lambda iW: beta * kronecker(iW, complex(0.,0.))
    else:
        return lambda iW: ( rho[j] - rho[i] ) / (iW + E[i] - E[j])

def chi2pp(block, iW):

    if block == ('up','dn'):
        # chi2pp = - <updn|cdag_up cdag_dn|0> <0|c_up c_dn|updn>
        # Take care of state ordering inside |updn>, i.e. c_dn|updn> = -|up>
        return f("updn","0")(iW)
    elif block == ('dn','up'):
        # chi2pp = - <i|cdag_dn cdag_up|j> <j|c_dn c_up|i>
        # Take care of state ordering inside |updn>, i.e. c_dn |updn> = -|up>
        return f("updn","0")(iW)
    elif block == ('up','up'):
        # chi3pp_upup = 0, There is no state for which <m| cdag_up(tau^{+}) cdag_up(tau) c_up(0^{+}) c_up(0) |m> or
        # or  <m| cdag_up(tau') cdag_up(tau) c_up(0^{+}) c_up(0) |m> is finite
        return 0.
    elif block == ('dn','dn'):
        # chi3pp_dndn = 0, There is no state for which <m| cdag_dn(tau^{+}) cdag_dn(tau) c_dn(0^{+}) c_dn(0) |m> or
        # <m| cdag_dn(tau') cdag_dn(tau) c_dn(0^{+}) c_dn(0) |m> is finite
        return 0.
    else:
        print "Invalid Block"
        raise

# -------- Initialize chi3pp_iw Container ---------

iW_mesh = MeshImFreq(beta, 'Boson', n_iw)
chi2pp_iw = Block2Gf_from_struct(iW_mesh, struct=gf_struct)

for iW in chi2pp_iw[('up','dn')].mesh:
    for block in [('up','up'), ('up','dn'), ('dn','dn'), ('dn','up')]:
        chi2pp_iw[block][iW] = chi2pp(block, iW)

# -------- Save in archive ---------
res["chi2pp_iw"] = chi2pp_iw


# ================ chi2ph_iw ==================

# Exact expressions for Equal-time two-particle Green function
# chi2_ph(iW) = <(cdag c)(iW) n(0)>
# Compare https://hal.archives-ouvertes.fr/tel-01247625v1 Page 217

def chi2ph(block, iW):

    if block == ('up','dn'):
        # chi2ph_updn = <i|n_up|j> <j|n_dn|i>
        return f("updn", "updn")(iW)
    elif block == ('dn','up'):
        # chi2ph_dnup = <i|n_dn|j> <j|n_up|i>
        return f("updn", "updn")(iW)
    elif block == ('up','up'):
        # chi2ph_dnup = <i|n_up|j> <j|n_up|i>
        return f("up", "up")(iW) + f("updn", "updn")(iW)
    elif block == ('dn','dn'):
        # chi2ph_dnup = <i|n_dn|j> <j|n_dn|i>
        return f("dn", "dn")(iW) + f("updn", "updn")(iW)
    else:
        print "Invalid Block"
        raise

# def chi2ph_conn(block, iW):
    # return 0.0

# -------- Initialize chi2ph_iw Container ---------

chi2ph_iw = Block2Gf_from_struct(mesh=iW_mesh, struct=gf_struct)

for iW in chi2ph_iw[('up','dn')].mesh:
    for block in [('up','up'), ('up','dn'), ('dn','dn'), ('dn','up')]:
        chi2ph_iw[block][iW] = chi2ph(block, iW)

# -------- Save in archive ---------
res["chi2ph_iw"] = chi2ph_iw


# ================ Two-particle Green function G2_iw ==================

# Exact expressions for two-particle Green function according to
# Eq. (B.11) and (B.12) of S. Pairault, et al., EPJ B 16, 85-105 (2000)

# Function to return G2c_iw('up','up') and G2c_iw('up','dn')
# in frequency notation according to Paper: <T c_1 c_2 cdag_3 cdag_4>
def G2c_iw_EPJ(iw_1, iw_2, iw_3):

    iw_4 = iw_1 + iw_2 - iw_3

    x_1_up = iw_1 + mu + h
    x_1_dn = iw_1 + mu - h
    x_2_up = iw_2 + mu + h
    x_2_dn = iw_2 + mu - h
    x_3_up = iw_3 + mu + h
    x_3_dn = iw_3 + mu - h
    x_4_up = iw_4 + mu + h
    x_4_dn = iw_4 + mu - h

    x_1_bar_up = iw_1 + mu + h - U
    x_1_bar_dn = iw_1 + mu - h - U
    x_2_bar_up = iw_2 + mu + h - U
    x_2_bar_dn = iw_2 + mu - h - U
    x_3_bar_up = iw_3 + mu + h - U
    x_3_bar_dn = iw_3 + mu - h - U
    x_4_bar_up = iw_4 + mu + h - U
    x_4_bar_dn = iw_4 + mu - h - U

    G2c_iw_upup = beta * U**2 * n['dn'] * (1. - n['dn']) * (kronecker(iw_2, iw_3) - kronecker(iw_1, iw_3)) / (x_1_up * x_1_bar_up * x_2_up * x_2_bar_up)

    if abs(U - 2. * mu) < 1E-14:
        G2c_iw_updn = kronecker(iw_1, -iw_2) * (3. * beta + beta * exp((mu - h) * beta) + beta * exp((mu - h) * beta)) \
        / (2. + exp((mu + h) * beta) + exp((mu - h) * beta))**2 * (1. / x_1_bar_dn + 1. / x_2_bar_up) * (1. / x_3_bar_up + 1. / x_4_bar_dn)
    else:
        G2c_iw_updn = (n['up'] + n['dn'] - 1.) / (iw_1 + iw_2 + 2. * mu - U) * (1. / x_1_bar_dn + 1. / x_2_bar_up) * (1. / x_3_bar_up + 1. / x_4_bar_dn)

    if abs(h) < 1E-14:
        G2c_iw_updn -= kronecker(iw_1, iw_3) * beta * exp(mu * beta) / (1. + 2. * exp(mu * beta) + exp((2 * mu - U) * beta)) \
          * (1. / x_1_dn - 1. / x_3_bar_up) * (1. / x_4_up - 1. / x_2_bar_up)
    else:
        G2c_iw_updn += (n['up'] - n['dn']) / (iw_1 - iw_3 - 2. * h) * (1. / x_1_dn - 1. / x_3_bar_up) * (1. / x_4_up - 1. / x_2_bar_up)

    G2c_iw_updn += beta * U**2 * kronecker(iw_2, iw_3) * (exp((2. * mu - U) * beta) - exp(2. * mu * beta)) / Z**2 * 1. / (x_1_dn * x_1_bar_dn * x_2_up * x_2_bar_up)

    G2c_iw_updn += (n['up'] - 1.) / (x_1_dn * x_3_bar_up * x_4_dn) + (1. - n['up']) / (x_1_dn * x_2_bar_up * x_3_bar_up) \
       + (1. - n['dn']) / (x_1_bar_dn * x_2_up * x_3_bar_up) + (n['dn'] - 1.) / (x_2_up * x_3_bar_up * x_4_dn)

    G2c_iw_updn += (1. - n['dn']) / (x_1_dn * x_2_up * x_4_dn) + (1. - n['dn']) / (x_1_dn * x_2_up * x_3_up)

    G2c_iw_updn += (1. - n['dn']) / (x_1_bar_dn * x_3_up * x_4_bar_dn) + (n['dn'] - 1.) / (x_1_bar_dn * x_2_up * x_3_up) \
       + (1 - n['dn']) / (x_2_bar_up * x_3_up * x_4_bar_dn) + (n['dn'] - 1.) / (x_1_dn * x_2_bar_up * x_3_up)

    G2c_iw_updn += -n['up'] / (x_1_bar_dn * x_2_bar_up * x_4_bar_dn) - n['up'] / (x_1_bar_dn * x_2_bar_up * x_3_bar_up)

    return G2c_iw_upup, G2c_iw_updn

# The full vertex
def F_iw_EPJ(iw_1, iw_2, iw_3):

    iw_4 = iw_1 + iw_2 - iw_3

    F_upup, F_updn = G2c_iw_EPJ(iw_1, iw_2, iw_3)

    F_upup *= 1. / (G['up'](iw_1) * G['up'](iw_2) * G['up'](iw_3) * G['up'](iw_4))
    F_updn *= 1. / (G['up'](iw_1) * G['dn'](iw_2) * G['up'](iw_3) * G['dn'](iw_4))

    return F_upup, F_updn

# The two-particle Green function
def G2_iw_EPSJ(iw_1, iw_2, iw_3):

    # Calculate connected part
    G2_upup, G2_updn = G2c_iw_EPJ(iw_1, iw_2, iw_3)

    # Add disconnected part
    iw_4 = iw_1 + iw_2 - iw_3
    G2_upup += beta * ( kronecker(iw_2, iw_3) - kronecker(iw_1, iw_3) ) * G['up'](iw_1) * G['up'](iw_2)
    G2_updn += beta * kronecker(iw_2, iw_3) * G['up'](iw_1) * G['up'](iw_2)

    return G2_upup, G2_updn
