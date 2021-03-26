# AnySim - Framework for solving arbitrary linear systems

## Class hierarchy
~~~
AnySim : root class implementing the Modified Born series iteration
 │
 └─ GridSim : base class for all simulations on a grid
     │
     ├─ DiffuseSim : solves the diffusion equation
     └─ HelmholtzSim : solves the Helmholtz equation

Medium : root class implementing operator 1-G for grid-based potentials
 │       also implements computation of centering & scaling matrices Tl Tr and V0
 ├─ TensorMedium: associates a matrix to each grid point, operator V=1-G is implemented as a matrix-vector multiplication
 ├─ DiagonalMedium: associates a diagonal matrix to each grid point, functionally equivalent to TensorMedium but more efficient
 └─ ScalarMedium: associates a scalar to each grid point
   
Classes that implement the Transform interface:
   FourierTransform

Helper classes:
   SimGrid
   DisplayCallback
   TerminationCondition
   Source
   State
~~~

## Example: Diffusion equation

## Implementing a solver
The general structure of the modified Born series algorithm is implemented already. To implement a solver for a specific linear problem, one needs to implement a simulation object inheriting from AnySim or one of its derived classes (such as GridSim). See DiffuseSim for an example. The following methods should be implemented:
  * The constructor
  * AnalyzeDimensions
  * Start
 
### The constructor
The constructor takes all information needed to describe a specific linear system (e. g. a refractive index map) and a set of options. The returned object fully describes the linear system and all details of the simulation. 

The constructor should check the validity of all inputs and fill in defaults for missing options. Importantly, it should set the properties `medium`, `propagator`, and `transform`. 

Medium: This object corresponds to operator `G:= 1-V = 1-Tl (V_raw-V0) Tr`, and typically is
  implemented as a multiplication with a scattering potential in the 
  spatio-temporal domain.
  For the diffusion equation, for example, a DiffusionMedium object is used,
  which performs a multiplication with the absorption-(inverse)diffusion tensor.
  In preparing the medium, the scattering potential is first shifted by ~V0~ to minimze ‖V‖, and then scaled to have ‖V‖<1. The matrices responsible for this scaling (Tl and Tr) are stored in the Medium object.

Propagator: This object corresponds to scaled operator `Tl (L+1)^(-1) Tr`, and typically
  is implemented as a multiplication with a fixed function in the spatio-temporal frequency domain.
  For the diffusion equation, the DiffusePropagator object corresponds
  to a spatio-temporal differential operator, applied in the frequency
  domain.

Transform: This operator transforms between the domains of V and L
  Typically, this is just a Fourier transform.


## Algorithm implementation
### Preconditioning
    The initial linear equation (L+V) u=s is first converted to
    Tl^(-1) (L'+V') Tr^(-1) u = s

    with: L' = Tl (L+V0) Tr       and   V' = Tl (V-V0) Tr
    introducing: u' = Tr^(-1) u   and   s' = Tl s
    we get the equation in the form:
    (L'+V') u' = s'

    with Tl and Tr diagonal invertible pre-conditioning matrices. V0
    is an operator. Tl, Tr, and V0 are chosen such that
        % 1. L', V' and are dimensionless
        % 2. The operator norm ||V-V0|| is minimized
        % 3. ||V'|| < 1 (but as close to 1 as conveniently possible)
        %
    
    Implementation:
    Tl, V0 and Tr are computed by the Medium object and stored in the
    'scaling' structure property.
    The Medium operator stores V'
    The Source operator stores s'
    The Propagator operator stores G' = (L' + 1)^(-1)
    After the iterations finishes, the algorithm returns u = Tr u'

### Iteration
In the manuscript, the following iteration is derived:
    u -> (G Li G + 1 - G) u + G Li s

which is implemented as:
    t1 = G u + s            Medium.mix_source
    t1 -> Li t1             Propagator.propagate
    u -> u + G (t1 - u)     Medium.mix_field

requires:
1 temporary storage (t1)
1 field storage (u)
1 potential storage (G)
1 propagator storage (Li, may be computed on the fly in some implementations!)


### Bookkeeping
While running the algorithm, a State object is passed to all operator
function calls. The State object is of handle type (a reference) and
contains bookkeeping information (such as the iteration number) and can be
used to store diagnostics and debugging information.


## GridSim
### Data Storage
GridSim objects work with simulation data that can be represented on a
regular grid. The data may be scalar, vector, or matrix-valued.

Internally, all data is stored as a matrix field in an N-dimensional array 'u'.
The first two dimensions of 'u' correspond to the size of a single value.
For scalar simulations, these dimensions are [1,1]. For vector-valued
data, the dimensions are [N_components,1], and for matrix-valued data
they are [N, M]. Matrix-valued operators (G and V) follow the same data layout.

When data is passed to the user (when returning from exec()), spurious
dimensions are removed. So, a Nx x Ny scalar simulation will return a Nx x Ny
array.

Each dimension may have a different pixel pitch and unit, this metadata
is stored in a SimGrid object.
