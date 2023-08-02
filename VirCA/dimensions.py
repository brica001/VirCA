# Reference units for the fundamental dimensions in the model (SI based)

# Use reference units based on relevant inputs
import inputs as inp
import constants as ct

L = inp.Rs ### Length: subunit radius
E = ct.kB*inp.T ### Energy: kT for input temperature
t = 1/(4*ct.pi*(inp.Rs+inp.Rs)*(inp.D1+inp.D1)*inp.c0)  ### Time: inverse of subunit encounter rate.
C = inp.c0 ### Concentration: initial subunit concentration
