# AnySim - Framework for solving arbitrary linear systems

## Class hierarchy
~~~
AnySim : root class implementing the Modified Born series iteration
 |
 |- GridSim : base class for all simulations on a grid
     |
     |- DiffuseSim : solves the diffusion equation

TensorMedium: object implementing the 'Medium' interface for tensor fields
 |
 |- DiffuseMedium : constructs a medium for DiffuseSim using absorption and
                    diffusion coefficients
~~~

## Implementing a solver
The general structure of the modified Born series algorithm is implemented already. To implement a solver for a specific linear problem, one needs to implement a simulation object inheriting from AnySim or one of its derived classes (such as GridSim). See DiffuseSim for an example. The following methods should be implemented:
  * The constructor
  * AnalyzeDimensions
  * Start
  
to four different operators:
- Medium: This object corresponds to operator V, and typically is
  implemented as a multiplication with a scattering potential in the 
  spatio-temporal domain.
  For the diffusion equation, for example, a DiffusionMedium object is used,
  which performs a multiplication with the absorption-(inverse)diffusion tensor.

- Propagator: This object corresponds to operator (L+1)^(-1), and typically
  corresponds to a multiplication with a fixed function in the
  spatio-temporal frequency domain.
  For the diffusion equation, the DiffusePropagator object corresponds
  to a spatio-temporal differential operator, applied in the frequency
  domain.

- Transform: This operator transforms between the domains of V and L
  Typically, this is just a Fourier transform.

- Source: This operator applies the source term

== Algorithm implementation ==
= Preconditioning = 
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

= Iteration =
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


= Bookkeeping = 
While running the algorithm, a State object is passed to all operator
function calls. The State object is of handle type (a reference) and
contains bookkeeping information (such as the iteration number) and can be
used to store diagnostics and debugging information.


== GridSim==
= Data Storage =
Grid-based simulations store data as multi-dimensional arrays, where the 
dimensions correspond to spatial and temporal dimensions of the simulated
problem. 

Each dimension may have a different pixel pitch and unit, this metadata
is stored in a SimGrid object.

When working with vector or tensor fields fields, the vector index/tensor
index corresponds to the _first_ dimension in the array (first two
dimensions for tensor array).
For example, an NxM field of 3-element vectors is stored as a 3 x N x M
multidimensional array. 
This storage scheme has the advantage of better memory locality for matrix
multiplications and it is slightly easier to program in MATLAB (e.g.
one can use pagefun for gpu arrays). 
The disadvantages are 1) if the vector length is not a power of 2, reads
and writes are not aligned and not coalesced. (although the performance hit
is not too bad on modern hardware?)
2) the fft may be less efficient for strided data.

== ==

See test_diffusion for an example 


ISSUES: L+V
