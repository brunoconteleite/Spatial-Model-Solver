
# This is a test for solving a system of non-linear equations for "my" spatial
# model of economy and geography. The model is based upon Allen, Arkolakis (2014).
# I start with a version of my model without agglomeration/congestion effects on
# amenity utility. This is a "economy on a line" version of the model!

from scipy import optimize
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

# Setting working directory - start from user (~/) from Mac:
os.chdir('Dropbox/Research/jm project/')

# Setting rendering of plots to allow for Latex
rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

# --------------
# Inputs:

    # * Geography G(S):
    # N: number of locations;
    # L0 = {L_i^0}: vector of initial labor allocation;
    # {A_i}: vector of productivities;
    # T = {tau_ij}: matrix of bilateral trade costs;
    # M = {mu_ij}: matrix of bilateral migration frictions.

    # * Parameters (for the moment):
    # theta: dispersion of Frechet;
    # s: CES.



# --------------
# 1. Defining functions:

# TAU matrix:
def tau(t,N):
    "Calculates the bilateral trade costs' matrix"
    "Parametrization: {tau_ij} = {e^(tau*|i-j|)}"
    T = np.nan * np.ones((N,N))
    for i in range(N):
        for j in range(N):
            T[i,j] = np.exp(t*np.absolute(i-j))
    return T

# P_i:
def ds_price(w,A,T):
    "This calculates the Dixit-Stiglitz price index P_i,"
    "already vectorized for all locations!; thus the outcome"
    "is a Nx1 array"
    p = np.divide(w,A)
    p = np.power(p,(1-s))
    p = np.matmul(p,
        np.power(T,(1-s)))
    p = np.power(p,(1/(1-s)))
    return p

# Pi shares -- I calculate it as a NxN symmetric matrix in which every element
# is a PI_ij.
def Pi(w,A,M,T):
    "Calculates the shares of migration as a NxN transition matrix"
    "The formula in matrix algebra is"
    "Pi = phi x Phi x T"
    p = ds_price(w,A,T)
    phi = np.power(np.divide(w,p),theta)    # this is the "partial" nominator
    Phi = np.power(M,(-theta))              # this is the denominator
    Phi = np.matmul(phi,Phi)
    Phi = np.power(Phi,(-1))
    Pi = np.outer(phi,Phi)                  # this is the matrix of bilateral shares
    Pi = np.multiply(Pi,np.power(M,(-theta)))
    return Pi

def Labor(w,A,M,T,L0):
    "Calculates the labor allocation in equilibrium"
    P = Pi(w,A,M,T)
    return np.matmul(P,L0)

def eq_condition(w,A,M,T,L0):
    "Evaluates equilibrium condition at a certain level of {w} given G(S)"
    "output: vector of Nx1 equations that must equal zero in equilibrium"
    L = Labor(w,A,M,T,L0)
    P = ds_price(w,A,T)
    omega = np.multiply(L,np.power(A,(1-s)))
    omega = np.power(omega,(-1))
    psi = np.multiply(w,L)
    psi = np.divide(psi,np.power(P,(1-s)))
    psi = np.matmul(np.power(T,1-s),psi)
    return w - np.multiply(omega,psi)


# --------------
# 3. All seems set to look for a solution! I will use a GMM-like procedure
# to try to minimize the square of the objective function.

def f(w,A,M,T,L0):
    "Objective function to be minimized"
    return np.matmul(eq_condition(w,A,M,T,L0),eq_condition(w,A,M,T,L0))

# setting up parameters:
N = 31
L0 = np.ones(N)
A = np.ones(N)
# A2 = np.linspace(1,2,N)

# s = 3 # by Sotelo
s = 6.4 # by Ruggieri
theta = 7 # by me
# theta = 7 # by Morten, Oliveira

# To create trade and migration frictions, I follow standard parametrization by
# setting them as
t = 0.01 # tau
T = tau(t,N)       # bilateral frictions' matrix
# M = T            # bilateral migration frictions' matrix
M = np.ones((N,N))

# Benchmark optimization:
w0 = np.ones(N)
w = optimize.fmin(f,w0, xtol=1e-8, args=(A,M,T,L0,), maxfun=1e10)
# w2 = optimize.fmin(f,w, xtol=1e-8, args=(A2,M,T,L0,), maxfun=1e10)
# w2 = optimize.fmin(f,w0, xtol=1e-8, args=(A,M,T2,L0,), maxfun=1e10)
# w3 = optimize.fmin(f,w0, xtol=1e-8, args=(A,M,T3,L0,), maxfun=1e10)

# Benchmark set of equilibrium allocations:
w1 = w
L1 = Labor(w1,A,M,T,L0)
P1 = ds_price(w1,A,T)

# --------------
# 4. Some experiments with different parameters - use previous solution as
# initial values:

# Increasing trade costs:
T2 = tau(t+0.02,N)

w2 = optimize.fmin(f,w, xtol=1e-8, args=(A,M,T2,L0,), maxfun=1e10)
L2 = Labor(w2,A,M,T2,L0)

# Uneven productivities (initial trade costs):
A2 = np.linspace(1,1.05,N)
w3 = optimize.fmin(f,w, xtol=1e-8, args=(A2,M,T,L0,), maxfun=1e10)
L3 = Labor(w3,A2,M,T,L0)

# Different dispersion of taste shocks:
theta=1.1
w4 = optimize.fmin(f,w, xtol=1e-8, args=(A,M,T2,L0,), maxfun=1e10)
L4 = Labor(w4,A,M,T2,L0)

# --------------
# 5. Plotting results:

# baseline plot:
plt.plot(np.linspace(0,1,N), L1, label=r'$L_i$')
plt.plot(np.linspace(0,1,N), P1, '--', label=r'$P_i$')
plt.plot(np.linspace(0,1,N), w1, ':', label=r'$w_i$')
plt.legend()
plt.ylabel(r'$L_i$, $P_i$, $w_i$')
plt.xlabel(r'Location $i$')
# plt.show()
plt.savefig('results/figures/model_eq_labor.png', dpi=300, bbox_inches="tight")
plt.clf()

# different trade costs:
plt.plot(np.linspace(0,1,N), L1, label=r'$\tau$=0.01')
plt.plot(np.linspace(0,1,N), L0, ':', label=r'$\tau$=0')
plt.plot(np.linspace(0,1,N), L2, '--', label=r'$\tau$=0.03')
plt.legend()
plt.ylabel(r'$L_i$')
plt.xlabel(r'Location $i$')
# plt.show()
plt.savefig('results/figures/model_eq_labor_diff_tradecosts.png', dpi=300, bbox_inches="tight")
plt.clf()

# uneven productivities:
plt.plot(np.linspace(0,1,N), L1, label=r'Symmetric $A_i$')
plt.plot(np.linspace(0,1,N), L3, '--', label=r'Uneven $A_i$')
plt.plot(np.linspace(0,1,N), A2, ':', label=r'$A_i$')
plt.legend()
plt.ylabel(r'$L_i$')
plt.xlabel(r'Location $i$')
# plt.show()
plt.savefig('results/figures/model_eq_labor_uneven_Ai.png', dpi=300, bbox_inches="tight")
plt.clf()

# Different theta:
plt.plot(np.linspace(0,1,N), L2, label=r'$\theta=7$')
plt.plot(np.linspace(0,1,N), L4, '--', label=r'$\theta=2$')
plt.legend()
plt.ylabel(r'$L_i$')
plt.xlabel(r'Location $i$')
# plt.show()
plt.savefig('results/figures/model_eq_labor_diffe_theta.png', dpi=300, bbox_inches="tight")
plt.clf()





# # Mg Utility:
# def mg_utility(i):
#     p = ds_price(i);
#     p2 = w[i];
#     return p2/p
#
# # \Pi_ij:
# def Pi(i,j):
#     "This calculates the share of workers that migrate from i to j"
#     return
# # Mg Utility:
# def mg_utility(i):
#     p = ds_price(i);
#     p2 = w[i];
#     return p2/p
#
# # \Pi_ij:
# def Pi(i,j):
#     "This calculates the share of workers that migrate from i to j"
#     return
# # Mg Utility:
# def mg_utility(i):
#     p = ds_price(i);
#     p2 = w[i];
#     return p2/p
#
# # \Pi_ij:
# def Pi(i,j):
#     "This calculates the share of workers that migrate from i to j"
#     return
# # Mg Utility:
# def mg_utility(i):
#     p = ds_price(i);
#     p2 = w[i];
#     return p2/p
#
# # \Pi_ij:
# def Pi(i,j):
#     "This calculates the share of workers that migrate from i to j"
#     return
