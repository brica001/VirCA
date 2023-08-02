# Import libraries
import numpy as np

# Define alternative equations
#def pend(y, t, b, c):
#    theta, omega = y
#    dydt = [omega, -b*omega - c*np.sin(theta)]
#    return dydt

def capsid_q_2(c0,t,k1,a2):
    c1, c2 = c0
    dc1dt = -2*k1*c1*c1 + a2*c2
    dc2dt = +k1*c1*c1 - a2*c2 
    dcdt = [dc1dt, dc2dt]

    return dcdt

def test(c0,t,v):
    c1, c2 = c0
    dc1dt = -2*v[0]*c1*c1 + v[1]*c2
    dc2dt = +v[0]*c1*c1 - v[1]*c2
    dcdt = [dc1dt, dc2dt]

    return dcdt

def capsid_q(c0,t,q,k,a):
    # q: number of subunits
    # k: vector with the association rate constant of a monomer to each partial capsid
    # a: vector with dissociation rate constant of each partical capsid into a monmer and a partial capsid
    
    # Initial conditions
    c = c0

    # Declare vector containing rates
    dcdt = np.zeros(q)
    
    # Compact system of differential equations
    dcdt[0] = -2*k[0]*c[0]*c[0] + 2*a[1]*c[1] - c[0]*np.sum(np.multiply(k[1:q-1],c[1:q-1])) + np.sum(np.multiply(a[2:],c[2:]))
    dcdt[1:q-1] = c[0]*(np.multiply(k[0:q-2],c[0:q-2])) - c[0]*(np.multiply(k[1:q-1],c[1:q-1])) - np.multiply(a[1:q-1],c[1:q-1]) + np.multiply(a[2:],c[2:])
    dcdt[q-1] = c[0]*k[q-2]*c[q-2] - a[q-1]*c[q-1]

    return dcdt
