import numpy as np
import os
import timeit
import model_multisector_CD as model

def optimal_scale(a,gdp,w,L,A0,T,s,theta,mu):
    '''
    Function to be minimized so to find the optimal scale that links productivity
    data in the units I collected with GDP data measured in US$. Need as inputs
    the agricultural productivities (A) and the agricultural gdp (gdp). The rest
    is standard.
    '''
    A = np.multiply(a,A0)
    P = model.ds_price(w,A,T,s)
    nc = np.int(np.size(A)/len(A)) # number of sectors (columns of A)
    Phi = np.repeat(np.multiply(w,L),nc).reshape(A.shape)
    Phi = np.divide(Phi, np.power(P,1-s))
    Phi = np.matmul(np.power(T,(1-s)), Phi)
    Psi = np.repeat(w,nc).reshape(A.shape)
    # Psi = np.array(Psi,dtype=np.complex)    # add this to avoid problems in the power
    Psi = np.power(np.divide(Psi,A),(1-s))
    # Psi = np.real(Psi) # back to real numbers
    Psi = np.multiply(Psi,mu)
    f = gdp-np.sum(np.multiply(Psi,Phi))
    return f

def demand_diff(a,w,L,A0,T,s,theta,mu):
    '''
    Evaluates equilibrium condition at a certain level of {w} given G(S)
    output: vector of Nx1 equations that must equal zero in equilibrium.
    Takes labour as given (not calculated). A0 is the vector of known A_i,
    "a" is the vector of A_i^K (to be calculated).
    '''
    A = np.column_stack((A0,a))
    P = model.ds_price(w,A,T,s)
    nc = np.int(np.size(A)/len(A)) # number of sectors (columns of A)
    Phi = np.repeat(np.multiply(w,L),nc).reshape(A.shape)
    Phi = np.divide(Phi, np.power(P,1-s))
    Phi = np.matmul(np.power(T,(1-s)), Phi)
    Psi = np.repeat(w,nc).reshape(A.shape)
    # Psi = np.array(Psi,dtype=np.complex)    # add this to avoid problems in the power
    Psi = np.power(np.divide(Psi,A),(1-s))
    # Psi = np.real(Psi) # back to real numbers
    Psi = np.multiply(Psi,mu)
    if nc>1:
        f = np.sum(np.multiply(Psi,Phi),axis=1)
    else:
        # in the case of one sector, the shapes get a bit messy: (N,1) rather than (N,)
        f = np.multiply(Psi.reshape(w.shape),Phi.reshape(w.shape))
    f = np.multiply(w,L) - f
    return f

def match_A(w,L,a,A0,T,s,theta,mu,p,ww,maxiter):
    '''
    Iterates over initial A guess to find model implied one. The iteration rule
    works as follows:
    - I calculate difference in the observed and calculated demands using the
      model and the data with demand_diff();
    - If the maximum of the difference is higher than the tolerance (1e-p),
      iterate:
      * A1 = A0 * demand_diff*(weight), where weight gets smaller the closest we
      get to the tolerance.

     Inputs:
     a,A0: initial guess of unobserved productivity, oberved product;
     p: precision of the tolerance/stopping criterium (1e-p);
     ww: weight to the dimishing rule of the weight;
     maxiter: maximum iterations to avoid endless loops.
    '''
    # Starting...
    ttt = timeit.default_timer()
    t = 10**(-p)
    cc = 0 # counter

    # output vectors:
    zz=np.repeat(np.nan,(np.int(maxiter+1))*w.shape[0]).reshape(w.shape[0],np.int(maxiter+1)) # store the differences
    xx=np.repeat(np.nan,maxiter+1) # store the tolerances
    aa=np.repeat(np.nan,(np.int(maxiter+1))*w.shape[0]).reshape(w.shape[0],np.int(maxiter+1)) # store the productivities

    # check and store demand diff.
    c = demand_diff(a,w,L,A0,T,s,theta,mu)
    zz[:,cc] = c
    aa[:,cc] = a

    # calculate and store tolerance:
    # tol = np.amax(np.absolute(c)) # all differences super close to zero
    # tol = np.sum(np.absolute(c)) # a bit meaningless
    tol = np.abs(np.sum(c)) # excess demand close to zero (MUST BE SUUUUPER PRECISE!)
    xx[cc] = tol

    while tol>t:
        # creating the weight on the iterative updating:
        # we = np.geomspace(1,10**(-p),num=p+1)
        # we = we[tol>we][0]*(ww)
        we = np.power(ww,2)

        # # adjustment: two possible types:
        # # - initial one with adding a weigthed demand diff:
        # # adj = np.add(a,np.multiply(c,we))
        # adj = np.add(a,c)
        # a = adj

        # - second: using exponentials and weigthing at each step:
        # adj = np.multiply(np.divide(c,np.abs(c)),we*np.power(np.abs(c),we)) # super small steps
        adj = np.multiply(np.divide(c,np.abs(c)),ww*np.power(np.abs(c),we)) # smal steps
        # adj = np.multiply(np.divide(c,np.abs(c)),np.power(np.abs(c),we)) # larger steps
        adj = np.add(a,adj)
        # adding weights:
        a = np.multiply(adj,ww)+np.multiply(a,(1-ww))

        # need to restrict the productivities to be > 0!
        a = a.clip(0.00000001)

        # re-calc. demand difference:
        c = demand_diff(a,w,L,A0,T,s,theta,mu)

        # tolerance update:
        # tol = np.amax(np.absolute(c))
        # tol = np.sum(np.absolute(c))
        tol = np.abs(np.sum(c)) # excess demand close to zero (MUST BE SUUUUPER PRECISE!)
        cc = cc+1

        # storing output:
        xx[cc] = tol
        zz[:,cc] = c
        aa[:,cc] = a
        if cc==maxiter:
            break
        # if cc in np.arange(0,maxiter,1000):
        #     print("Iterating (" + str(cc) + ")...")
    if tol<t:
        print("Iteration terminated successfully")
        print("    Current tolerance value: " + str(tol))
        print("    Iterations: " + str(cc))
        print("    Running time (secs): " + str(np.round(timeit.default_timer()-ttt, decimals = 4)))
    else:
        print("Iteration terminated UNsuccessfully")
        print("    Current tolerance value: " + str(tol))
        print("    Iterations: " + str(cc))
        print("    Running time (secs): " + str(np.round(timeit.default_timer()-ttt, decimals = 4)))
    return a, xx, zz, aa

def epsilon(w,L,A,T,s,theta,mu):
    '''
    explain here
    '''
    P = model.ds_price(w,A,T,s)
    nc = np.int(np.size(A)/len(A)) # number of sectors (columns of A)
    Phi = np.repeat(np.multiply(w,L),nc).reshape(A.shape)
    Phi = np.divide(Phi, np.power(P,1-s))
    Phi = np.matmul(np.power(T,(1-s)), Phi)
    Psi = np.repeat(w,nc).reshape(A.shape)
    # Psi = np.array(Psi,dtype=np.complex)    # add this to avoid problems in the power
    Psi = np.power(np.divide(Psi,A),(1-s))
    # Psi = np.real(Psi) # back to real numbers
    Psi = np.multiply(Psi,mu[:-1])
    if nc>1:
        f = np.sum(np.multiply(Psi,Phi),axis=1)
    else:
        # in the case of one sector, the shapes get a bit messy: (N,1) rather than (N,)
        f = np.multiply(Psi.reshape(w.shape),Phi.reshape(w.shape))
    f = np.multiply(w,L) - f
    return f

def special_price(w,a,T,s):
    '''
    This calculates the Dixit-Stiglitz price index P_i, already vectorized for all
    locations-sectors! The diff. from the other function is that it considers
    a = A^{s-1}, to avoid the inversion.
    '''
    nc = np.int(np.size(a)/len(a)) # number of sectors (columns of A)
    p = np.repeat(w,nc).reshape(a.shape)
    p = np.multiply(np.power(p,(1-s)),a)
    p = np.matmul(np.power(T,(1-s)),p)
    p = np.power(p,(1/(1-s)))
    # p = np.real(p)                      # back to real number again (actually complex thing should not be a problem)
    return p

def optimal_a(w,L,a,A,T,s,theta,mu):
    '''
    Calculate model implied optimal A_i^K for looking for a fixed point. Inputs
    are A (Nx(K-1) matrix, all agricultural productivities) and a (Nx1, the prods
    of the non-agric. sector)
    '''
    eps = epsilon(w,L,A,T,s,theta,mu) # this is the nominator, must be a Nx1 array
    # now the denominator (dependent on A_i^K, also Nx1):
    P = model.ds_price(w,a,T,s)
    # P = special_price(w,a,T,s)
    nc = np.int(np.size(a)/len(a)) # number of sectors (columns of a; must be one!)
    Phi = np.repeat(np.multiply(w,L),nc).reshape(a.shape)
    Phi = np.divide(Phi, np.power(P,1-s))
    Phi = np.matmul(np.power(T,(1-s)), Phi)
    Psi = np.repeat(w,nc).reshape(a.shape)
    # Psi = np.array(Psi,dtype=np.complex)    # add this to avoid problems in the power
    Psi = np.power(Psi,(1-s))
    # Psi = np.real(Psi) # back to real numbers
    Psi = np.multiply(Psi,mu[-1:])
    # in the case of one sector, the shapes get a bit messy: (N,1) rather than (N,)
    aa = np.multiply(Psi.reshape(w.shape),Phi.reshape(w.shape))
    aa = np.divide(eps,aa).clip(0.00000001)
    aa = np.power(aa,(1/(s-1)))
    return aa

def fixed_point_A(w,L,a,A,T,s,theta,mu,p,ww,maxiter):
    '''
    explain here
    '''

    # Starting...
    ttt = timeit.default_timer()
    t = 10**(-p)
    cc = 0 # counter
    er=False # error to be updated if we get into an increasing zone.

    # output vectors:
    xx=np.repeat(np.nan,maxiter+1) # store the tolerances
    zz=np.repeat(np.nan,(np.int(maxiter+1))*w.shape[0]).reshape(w.shape[0],np.int(maxiter+1)) # store the productivities

    # iterate to fixed point:
    tol = 1
    while tol > t:
        # calculating optimal A given initial guess:
        aa = optimal_a(w,L,a,A,T,s,theta,mu)

        # checking and storing tolerance:
        tol = np.max(np.abs(np.subtract(aa,a)))
        xx[cc] = tol

        # updating it using weights:
        a = ww*aa + (1-ww)*a
        a = a.clip(0.00000001)
        zz[:,cc] = a
        cc = cc+1

        if cc>100:
            er=(np.sum(xx[(cc-100):(cc)]-xx[(cc-100-1):(cc-1)]>0)==100)
            if er==True:
                break # stop if my iterating rule is increasing the tolerance for more than 100 iterations.

        if cc==maxiter:
            break

    if tol<t:
        print("Iteration terminated successfully")
        print("    Current tolerance value: " + str(tol))
        print("    Iterations: " + str(cc))
        print("    Running time (secs): " + str(np.round(timeit.default_timer()-ttt, decimals = 4)))
    else:
        print("Iteration terminated UNsuccessfully")
        print("    Current tolerance value: " + str(tol))
        print("    Iterations: " + str(cc))
        print("    Running time (secs): " + str(np.round(timeit.default_timer()-ttt, decimals = 4)))
    return a, xx, zz, er

def optimal_w(w,L0,A,u,T,s,theta,mu):
    '''
    Calculate model-implied optimal wages as a function of Nx1 wages. To be used
    in a fixed-point algorithm to find equilibrium income.
    '''
    P = model.ds_price(w,A,T,s)
    L = model.Labor(model.Pi(w,P,u,theta,mu),L0)
    nc = np.int(np.size(A)/len(A)) # number of sectors (columns of A)
    Phi = np.repeat(np.multiply(w,L),nc).reshape(A.shape)
    Phi = np.divide(Phi, np.power(P,1-s))
    Phi = np.matmul(np.power(T,(1-s)), Phi)
    Psi = np.repeat(np.power(L,(-1)),nc).reshape(A.shape)
    Psi = np.divide(Psi,np.power(A,(1-s)))
    Psi = np.multiply(Psi,mu)
    if nc>1:
        f = np.sum(np.multiply(Psi,Phi),axis=1)
    else:
        # in the case of one sector, the shapes get a bit messy: (N,1) rather than (N,)
        f = np.multiply(Psi.reshape(w.shape),Phi.reshape(w.shape))
    f = np.power(f,(1/s))
    return f

def fixed_point_w(w,L0,A,u,T,s,theta,mu,p,ww,maxiter):
    '''
    explain here
    '''

    # Starting...
    ttt = timeit.default_timer()
    t = 10**(-p)
    cc = 0 # counter
    er=False # error to be updated if we get into an increasing zone.

    # output vectors:
    xx=np.repeat(np.nan,maxiter+1) # store the tolerances
    # zz=np.repeat(np.nan,(np.int(maxiter+1))*w.shape[0]).reshape(w.shape[0],np.int(maxiter+1)) # store the productivities

    # iterate to fixed point:
    tol = 1
    while tol > t:
        # calculating optimal A given initial guess:
        aa = optimal_w(w,L0,A,u,T,s,theta,mu)

        # checking and storing tolerance:
        tol = np.max(np.abs(np.subtract(aa,w)))
        xx[cc] = tol

        # updating it using weights:
        w = ww*aa + (1-ww)*w
        # a = a.clip(0.00000001)
        # zz[:,cc] = a
        cc = cc+1

        if cc>100:
            er=(np.sum(xx[(cc-100):(cc)]-xx[(cc-100-1):(cc-1)]>0)==100)
            if er==True:
                break # stop if my iterating rule is increasing the tolerance for more than 100 iterations.

        if cc==maxiter:
            break

    if tol<t:
        print("Iteration terminated successfully")
        print("    Current tolerance value: " + str(tol))
        print("    Iterations: " + str(cc))
        print("    Running time (secs): " + str(np.round(timeit.default_timer()-ttt, decimals = 4)))
    else:
        print("Iteration terminated UNsuccessfully")
        print("    Current tolerance value: " + str(tol))
        print("    Iterations: " + str(cc))
        print("    Running time (secs): " + str(np.round(timeit.default_timer()-ttt, decimals = 4)))
    # return a, xx, zz, er
    return w, xx, er

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
