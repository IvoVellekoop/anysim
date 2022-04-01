clear all; close all;

% run all tests for symmetric systems
test_diffusion_diagonal
asym_results = table; clear table;
clearvars -except sym_results asym_results;

test_diffusion_tensor
asym_results(end+1) = table(2); clear table;
clearvars -except sym_results asym_results;

test_matrix
sym_results = table; clear table;
clearvars -except sym_results asym_results;

test_helmholtz_1D
asym_results(end+1) = table(2); clear table;
clearvars -except sym_results asym_results;

test_helmholtz_2D
asym_results(end+1) = table(2); clear table;
clearvars -except sym_results asym_results;

test_pantograph_complex
asym_results(end+1) = table(2); clear table;
clearvars -except sym_results asym_results;

%test_pantograph_analytical
%asym_results(end+1) = table(2); clear table;
%clearvars -except sym_results asym_results;

disp(sym_results)
disp(asym_results)

f = fopen("simulation_results.tex", "w+");
for line = sym_results
    fprintf(f, "%s\n", line);
end
fprintf(f, "\n\n");
for line = asym_results
    fprintf(f, "%s\n", line);
end
fclose(f);
