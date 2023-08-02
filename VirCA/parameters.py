# Parameter values for the model
## The parameters are determined from physical constants and physical inputs from the system.
import inputs as inp
import constants as ct
import matplotlib.pyplot as plt
import numpy as np
q = inp.q
Drel= inp.D_rel


dL=inp.Rh*2###interaction distance defined as 2 times hydrodynamic radius of subunits 
###along interaction axis (units: m)

L=np.exp(-inp.Dg_hr)*4*ct.pi*inp.Rc/q###interaction length perpindicular to interaction axis (units:m)

A=[L*dL*(np.sqrt((n+1)*(q-n+1))) for n in range(q)] ###total interaction area between subunits for each partial capsid
## Binding rates (1/s)
### Model based on effective Smoluchowski rate
k=[(Drel[n]/A[n]) for n in range(q)] #Relative to partial capsid, how many subunits pass through interaction area A

# plt.plot(k,label='k')
## Decay rates (1/s)
a=[k[n-1]*np.exp(inp.DelG[n]-inp.DelG[n-1]) for n in range(1,q)]
a.append(0)# completed capsids are stable and don't unbind or bind
a[0]=0 #free subunits do not decay, unecessary as index not called in master eqs, but true as a value in the array of unbdinding kernels.
plt.plot(a,label='a')
plt.title('binding and unbinding kernels')
plt.ylabel('1/sec')
plt.xlabel('n')
plt.ylim(-10,10)
plt.grid()
plt.legend()
plt.show()
