

import numpy as np

def mu_likelihood(n, mu):#Okay
    return np.sum(n*np.log(mu)-mu)

def mu_score(n,mu): #Okay
    return n/mu-1

def mu(c,r): 
    return np.einsum("ia,abt,jb->ijt",c,r,c) 

def likelihood(c,r,n):
    return mu_likelihood(n,mu(c,r))

def r_score(c,r,n=None,back=None):
    if back is None and n is None:
        raise TypeError("Either n or the backprop should be specified")
    if back is None:
        back = mu_score(n,mu(c,r))
    return np.einsum("ijt,ia,jb->abt",back,c,c)

def c_score(c,r,n=None,back=None):
    if back is None and n is None:
        raise TypeError("Either n or the backprop should be specified")
    if back is None:
        back = mu_score(n,mu(c,r))
    return np.einsum("kjt,gbt,jb->kg",back,r,c) + \
           np.einsum("ikt,agt,ia->kg",back,r,c)