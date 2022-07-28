%% Runs the examples and produces the benchmark data and figures used in the manuscript

clear all; close all;
chdir(fileparts(mfilename('fullpath')));
addpath("..");

test_diffusion
all_results = table; clear table;
clearvars -except all_results;

test_diffusion_tensor
all_results(end+1) = table(2); clear table;
clearvars -except all_results;

test_helmholtz_1D
all_results(end+1) = table(2); clear table;
clearvars -except all_results;

test_helmholtz_2D
all_results(end+1) = table(2);
all_results(end+1) = table(3); clear table;
clearvars -except all_results;

test_pantograph_complex
all_results(end+1) = table(2); clear table;
clearvars -except all_results;

test_pantograph_converge_diverge
%%
disp(all_results)

f = fopen("simulation_results.tex", "w+");
for line = all_results
    fprintf(f, "%s\n", line);
end
fclose(f);
