% run all tests for symmetric systems
test_diffusion
sym_results = table; clear table;

test_diffusion_diagonal
sym_results(end+1) = table(2); clear table;

test_diffusion_tensor
sym_results(end+1) = table(2); clear table;

test_matrix
sym_results(end+1) = table(2); clear table;

test_helmholtz
asym_results = table; clear table;

test_helmholtz_inhomogeneous
asym_results(end+1) = table(2); clear table;

test_helmholtz_logo
asym_results(end+1) = table(2); clear table;

test_pantograph
asym_results(end+1) = table(2); clear table;

test_pantograph_analytical
asym_results(end+1) = table(2); clear table;

disp(sym_results)
disp(asym_results)
