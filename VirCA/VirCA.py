# Main VirCA's python script
# Date (start): 2022-03-31
# Date (end): xxxx-xx-xx

# Standard libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Physical constants
import constants as ct

# Physical inputs
import inputs as inp

# Reference dimensions
import dimensions as dim
L = dim.L ### Length
E = dim.E ### Energy
t = dim.t ### Time (use lower case to avoid confusion with Temperature)
C = dim.C ### Concentration
# Capsid geometry


# Model parameters
# (Defined based on q, Deltag, T, subunit characteristics, capsid geometry...) 
import parameters as para
## Number of subunits 
q = inp.q

## Binding rates
#k1 = 0.25
k=para.k
## Decay rates


a = para.a
## Combine
#param = (k1,a2)
#param = [k1,a2]

# Set initial concentrations
c_t0 = np.zeros(q)
c_t0[0] = inp.c0
# c_t0[0] = 10000
#c_t0 = [inp.c0, 0.0]

# Define sampling time
t = np.linspace(0, 100, 10001)

# Obtain dynamics
#from mastereqs import capsid_q_2 as syseq
#from mastereqs import test as syseq
from mastereqs import capsid_q as syseq
#from mastereqs import test
#sol = odeint(pend, y0, t, args=(b, c))
#sol = odeint(capsid_q_2, y0, t, args=param)
sol = odeint(syseq, c_t0, t, args=(q,k,a,))
#sol = odeint(test, y0, t, args=(param,)
# print(sol,np.ones(len(sol)))
# Total subunits
nvec = np.arange(1,q+1)
ctot = np.zeros((len(sol)))
for m in range(0,len(sol)):
    ctot[m] = sum(np.multiply(nvec,sol[m,:]))
print(ctot)


# Plot solution
plt.plot(t, sol[:, 0], 'b', label='c1(t)')
plt.plot(t,sol[:,1],'g', label='c2(t)')
plt.plot(t, sol[:, -1], 'y', label='cq(t)')
# for n in range(0,q):
#     plt.plot(t, sol[:, n], label=n+1)

# plt.plot(t,ctot, label='ctot')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()



quit()

