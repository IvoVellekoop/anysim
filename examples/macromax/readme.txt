# Macromax implementation with SplitRichardson

Three examples are included that should run as-is:
 - `exampleScalar`: scalar wave equation
 - `exampleVectorial`: vectorial wave equation
 - `exampleAnistropic`: vectorial wave equation with birefringent material
once the code is loaded on the Matlab path. Convenience functions `setPath` and `resetPath` are provided for that.

The Macroscopic Maxwell solver algorithm is contained in solveMacroscopicMaxwell.m. Its core iteration is compared to the more general splitrichardson solver.
Type `helpwin spitricharson` for specifications and usage instructions.