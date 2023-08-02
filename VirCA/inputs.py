# List of inputs associated with the physical system to be modeled
import constants as ct
from scipy import constants as cts
import numpy as np

# Fundamental inputs in "human readable (hr)" (intuitive units) for capsid assembly.da 
q_hr = 90 ### Number of subunits forming the capsid. Units: Dimensionless
Dg_hr = -12### Mean binding energy of a subunit in the full capsid. Units: kB*T0 where T0 is normal temperature (293.15 K).
T_hr = 25 ### Temperature of the system. Units: Celsius.
Rs_hr = 3 ## Subunit effective radius. Units: nm.
c0_hr = 1000 ### Initial concentration of subunits. Units: nM (nano molar = nano mol/L)
Rh_hr=2.5 ### hydrodynamic radius estimate for folded globular protein (HBsAg,mw=24kDa)
Rc_hr=18 ###capsid radius (HBV) Units: nM
### Note: Rational for working in nM: 1 molecule in a bacteria (1 um^3 in volume) would be 1.66 nM (so ~1 nM is an insightful unit)

### Warning!!! Water viscosity given at T=25C. Module needed to include temperature dependence.
nu_hr = 0.890 ### Water viscosity. Units: cP (centiPoise).

# Converstion of inputs to SI units
# These are the values used for consistency in operations.
import constants as ct
from scipy.constants import convert_temperature
import matplotlib.pyplot as plt

q = q_hr
Dg = Dg_hr*ct.kB*ct.T0 ## See Dg_hr. Units: Joules (J)
T = convert_temperature(T_hr, 'Celsius', 'Kelvin') ## See T_hr. Units: Kelvin (K)
Rs = Rs_hr*1E-9 ## See Rs_hr. Units: meters (m)
nu = nu_hr*1E-3 ## See nu_hr. Units: Pascal x seconds (Pa*s)
c0 = c0_hr*1E-6#*cts.N_A ## See c0. Units: 1/meter^3 (1/m^3)
Rh = Rh_hr*1E-9 ## See Rh_hr. Units: meters (m)
Rc= Rc_hr*1E-9 ## See Rc_hr. Units: meters (m)

D1 = ct.kB*T/(6*ct.pi*nu*Rs) ### Subunit diffusion. Based on Stokes-Einstein equation. Units: m^2/sec
D=[D1/(i+1) for i in range(q)] ### Approximates diffusion for partial capsids decreasing by 1/n
D_rel=[D[n]+D[0]for n in range(q)] ### Approximates relative diffusion as D1+Dn

DelG = [n*Dg_hr-n*np.log(c0)+ct.pi/2*np.abs(Dg_hr)/np.sqrt(q)*np.sqrt(n*(q-n)) for n in range(q)] ### Total free energy of partial capsids
###Ignoring -ln(cs) and kT, they vanish since we are only interested in the unitless differences in energies 