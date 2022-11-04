%% Runs the examples and produces the benchmark data and figures used in the manuscript

clear all; close all;
chdir(fileparts(mfilename('fullpath')));
addpath("..");

default_options = struct(...
    callback = TextCallback(), ...
    forward_operator = true, ...
    precision = 'single', ...
    gpu_enabled = true,...
    termination_condition = TerminationCondition(tolerance = 1E-3, iteration_count = 25E4)...
);

legacy_options = default_options;
legacy_options.legacy_mode = true;

pantograph_options = default_options;
pantograph_options.precision = 'double'; % PantographF is using sparse matrices, which can only be of double precision at the moment

tests = {...
    "diffusion isotropic (homogeneous slab)", "test_diffusion", default_options;...
    "diffusion anisotropic (Fig.~\ref{fig:diffusion})", "test_diffusion_tensor", default_options;...
    "Helmholtz 1-D (glass plate, $n=1.5$)", "test_helmholtz_1D", default_options; ...
    "Helmholtz 2-D $\bbR$ bias (Fig.~\ref{fig:helmholtz})", "test_helmholtz_2D", legacy_options;...
    "Helmholtz 2-D $\bbC$ bias (Fig.~\ref{fig:helmholtz})", "test_helmholtz_2D", default_options;...
    "Helmholtz 2-D $\bbR$ bias, low-contrast", "test_helmholtz_2D_low_contrast", legacy_options;...
    "Helmholtz 2-D $\bbC$ bias, low-contrast", "test_helmholtz_2D_low_contrast", default_options;...
    "pantograph non-symmetrised (Fig.~\ref{fig:pantograph})", "test_pantograph", pantograph_options; ...
};

preconditioners = {...
    struct(preconditioner = "circulant"),...
    struct(preconditioner = "none"),...
    struct(preconditioner = "moborn"),...
    struct(preconditioner = "shift", preconditioner_shift = 0.5),...
    struct(preconditioner = "shift", preconditioner_shift = 1),...
    struct(preconditioner = "shift", preconditioner_shift = 2),...
};


all_results = cell(numel(preconditioners), numel(tests));

for p_i = 1:numel(preconditioners)
    pre = preconditioners{p_i};
    disp(pre);

    for t_i = 1:size(tests, 1)
        test = tests(t_i, :);
        results = run_test(test, pre);
        all_results{p_i, t_i} = results;
    end
end

%% Store in file
f = fopen("simulation_results.tex", "w+");
for p_i = 1:numel(preconditioners)
    pre = preconditioners{p_i};
    fprintf(f, "===\n%s\n===\n\n", jsonencode(pre));
    
    results = all_results{p_i, 1};
    for r_i = 1:numel(results)
        r = results(r_i);
        fprintf(f, "& %s", r.name);
    end
    fprintf(f, "\\\\\n\\hline\n");

    for t_i = 1:size(tests, 1)
        % write table data
        fprintf(f, "%s ", tests{t_i, 1});
        results = all_results{p_i, t_i};
        for r_i = 1:numel(results)
            r = results(r_i);
            if r.flag == 0
                fprintf(f, " & %d", r.iter);
            else
                fprintf(f, " & -"); % stagnated or diverged
            end
        end
        fprintf(f, "\\\\\n");
    end
    fprintf(f, "\n");
end
fclose(f);
edit simulation_results.tex

%     if opt.measure_time
%         data = data + sprintf('\\\\\n%s, exec. time [s]', strrep(filename, '_', ' '));
%         for r = rstore
%             if r.iter < Nit && r.flag == 0
%                 data = data + sprintf(" & %.1f", r.time);
%             elseif r.flag == 3
%                 data = data + sprintf(" & s"); % stagnates, may be due too low machine precision
%             else
%                 data = data + sprintf(" & -");
%             end
%         end
%     end


function r = run_test(test, pre)
    global option_override; %#ok<GVMIS> 
    option_override = test{3};
    fields = fieldnames(pre);
    for k = 1:numel(fields)
        option_override.(fields{k}) = pre.(fields{k});
    end
    disp(test{2});
    disp(pre);
    eval(test{2});
    r = results;
    r(end) = []; % remove the anysim test.
    r = struct2table(r);
    r = removevars(r, {'value'});
    disp(r);
    r = table2struct(r);
    close all;
end



