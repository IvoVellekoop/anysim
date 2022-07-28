# AnySim - Framework for solving arbitrary linear systems

## Class hierarchy
~~~
AnySim : root class implementing the Modified Born series iteration
 │
 └─ GridSim : base class for all simulations on a grid
     │
     ├─ DiffuseSim : solves the diffusion equation
     ├─ Pantograph : solves the pantograph equation
     └─ HelmholtzSim : solves the Helmholtz equation

Helper classes:
   Grid
   DisplayCallback
   TerminationCondition
   State
   *Options classes : contain a list of options (and code for validating them)
                      for constructing an object of the corresponding type
~~~

## Implementing a solver
The general structure of the modified Born series algorithm is implemented
already. To implement a solver for a specific linear problem, one needs to
implement a simulation object inheriting from AnySim or one of its derived
classes (such as GridSim). See DiffuseSim for an example. The following
methods should be implemented:
  * The constructor. Processes options and sets 'medium', 'transform' and 
    'propagator' properties to the correct operators.
 
### The constructor
The constructor takes all information needed to describe a specific linear
system (e. g. a refractive index map) and a set of options. The returned
object fully describes the linear system and all details of the simulation. 

The constructor should check the validity of all inputs and fill in
defaults for missing options. Importantly, it should set the properties
`medium` and `propagator`. Note that the medium_adj and propagator_adj 
properties are experimental at the moment. 

`.medium`: Returns a handle to a function implementing x -> (1-V) x.
  Typically, this is implemented as a multiplication with a scattering potential in the 
  spatio-temporal domain.
  For the diffusion equation, for example, the function
  performs a multiplication with the absorption-(inverse)diffusion tensor.
  In preparing the medium, the scattering potential is first
  shifted by ~V0~ to minimze ‖V‖, and then scaled to have ‖V‖<1. The matrices
  responsible for this scaling (Tl and Tr) are stored in the Medium object.

`.propagator`: Returns a handle to a function implementing x -> (L+1)^{-1} x
  Typically, this is implemented as a multiplication with a fixed function in the spatio-temporal frequency domain.
 

## Algorithm implementation
### Preconditioning
    The initial linear equation (L_raw+V_raw) u_raw=s_raw
    is first converted to (L+V) u = s

    with: L = Tl (L_raw+V_0) Tr     and   V = Tl (V_raw-V0) Tr
    introducing: u = Tr^(-1) u_raw   and  s = Tl s_raw

    with Tl and Tr diagonal invertible pre-conditioning matrices. V0
    is an operator. Tl, Tr, and V0 are chosen such that
        % 1. L, V and are dimensionless
        % 2. The operator norm ||V-V0|| is minimized
        % 3. ||V|| < 1 (but as close to 1 as conveniently possible)
    
    Implementation:
    Tl, V0 and Tr are computed using the center_scale function. Currently,
    for vector potentials Tl = Tr are diagonal matrices, and for scalar
    potentials Tl = 1 and Tr is a scalar.

### Iteration
In the manuscript, the following iteration is derived:
    u -> (B Li B + 1 - B) u + B Li s

which is implemented as:
    t1 = B u + s
    t1 -> Li t1
    u -> u + alpha B (t1 - u)

requires:
1 temporary storage (t1)
1 field storage (u)
1 potential storage (B)
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

When storing scalar fields, 'u' is an N-dimensional array. For vector fields
or tensor fields, dimensions are added at the start of the array. For example
to store a 3-element vector in a 3-D simulation: size(u) = [3, Nx, Ny, Nz]

Each dimension may have a different pixel pitch and unit, this metadata
is stored in a SimGrid object.
