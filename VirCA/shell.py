# Module calculating the properties of (partial) shells (or capsids) for different models.

# Libraries
import numpy as np

class ReturnValue:
    def __init__(self,polymer,length,surface,surfgeom,surffrac):
        self.polymer = polymer
        self.length = length
        self.surface = surface
        self.surfgeom = surfgeom
        self.surffrac = surffrac

def sphere(q,R,subar):
    # Object
    ## The function is defined as a class.
    ## The return values can be accessed accordingly.
    ## If the inputs are arrays, the outpus will be arrays accordingly

    # Model
    ## This model assumes that partial capsids are spherical caps.

    # Inputs
    ## q: number of subunits in the fully formed capsid.
    ## R: capsid radius.
    ## subar: area of the subunit exposed on the external side of the capsid

    #  Return values
    ## polymer: number of subunits in the (partial) shell.
    ## length: length of the boundary.
    ## surface: surface area of the (partial) shell based on the area of the subunits
    ## surfrac: relative surface with respect the geometrical surface of the shell.

    # Code
    ## Polymer
    n = np.arange(1,q+1)
    n.astype(int)
    polymer = n

    # Length
    length = 4*np.pi*R/q*np.sqrt(n*(q-n))
    ## Formula from Luque Santolaria, Thesis, 2011 (page 230)

    # Surface based on units
    surface = subar*n

    # Surface based on geometry
    costh = 1 - 2*n/q
    h = R*(1-costh) ## height
    surfgeom = 2*np.pi*R*h

    # Relative surface between the subunit model and geometrical model
    surffrac = surface/surfgeom

    return ReturnValue(polymer,length,surface,surfgeom,surffrac)

# Examples
## capsid = sphere(10,1)
## print(capsid.polymer)
## print(capsid.length)
## print(capsid.surface)
