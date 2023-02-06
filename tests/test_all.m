%% Runs the examples and produces the benchmark data and figures used in the manuscript

clear all; close all;

%%
chdir(fileparts(mfilename('fullpath')));
addpath("..");

default_options = struct(...
    callback = TextCallback(), ...
    forward_operator = true, ...
    precision = 'single', ...
    gpu_enabled = true,...
    termination_condition = TerminationCondition(tolerance = 1E-3, iteration_count = 30000)...%25E3)...
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
    "Helmholtz 2-D $\bbR$ bias, dielectric", "test_helmholtz_2D_low_contrast", legacy_options;...
    "Helmholtz 2-D $\bbC$ bias, dielectric", "test_helmholtz_2D_low_contrast", default_options;...
    "pantograph non-symmetrised (Fig.~\ref{fig:pantograph})", "test_pantograph", pantograph_options; ...
};

preconditioners = {...
    struct(preconditioner = "none"),...
    struct(preconditioner = "moborn"),...
    struct(preconditioner = "shift", preconditioner_shift = 0.5),...
    struct(preconditioner = "circulant"),...
};
%%

all_results = cell(numel(preconditioners), size(tests, 1));

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
save test_all_results.mat
%%
f = fopen("simulation_results.tex", "w+");
for p_i = 1:size(all_results, 1)
    tabulate_results(f, p_i, 0, preconditioners, all_results, tests);
end
for p_i = 1:size(all_results, 1)
    if preconditioners{p_i}.preconditioner == "shift"
        tabulate_results(f, p_i, 1, preconditioners, all_results, tests);
    end
end
fclose(f);
edit simulation_results.tex

%%

function tabulate_results(f, p_i, inner, preconditioners, all_results, tests)
    % construct a table of the results for preconditioner with index p_i,
    % and write that to file f. When inner = 1, shows the inner
    % iterations in the table (for shift preconditioner only)
    % marks the best result (fewest iterations) in bold.
    
    pre = preconditioners{p_i};
    results = all_results{p_i, 1};

    % temporary: remove bicgstab(l) column:
    results(4) = [];
    if inner
        name = pre.preconditioner + " (inner iterations)";
    else
        name = pre.preconditioner;
    end
    fprintf(f, "\n\\textbf{%s}\\\\\n", name);
    fprintf(f, "\\begin{tabular}{l|c|c|c|c|c|c|c}");

    for r_i = 1:numel(results)
        r = results(r_i);
        name = strrep(r.name, "Rich", "FP"); % temporary!
        fprintf(f, "& %s", name);
    end
    fprintf(f, "\\\\\n\\hline\n");

    for t_i = 1:size(all_results, 2)
        % write table data
        fprintf(f, "%s ", tests{t_i, 1});
        results = all_results{p_i, t_i};
        % temporary: remove bicgstab(l) column:
        results(4) = [];
        
        % find best (fastest and lowest number of iterations)
        best_rit = inf;
        best_t = inf;
        for r_i = 1:numel(results)
            r = results(r_i);
            if inner
                rit = r.iter_intern;
            else
                rit = r.iter;
            end
            if r.flag == 0
                best_rit = min(rit, best_rit);
                best_t = min(r.time, best_t);
            end
        end
        
        for r_i = 1:numel(results)
            fprintf(f, " & "); 
            r = results(r_i);
            if r.flag == 0
                if inner
                    rit = r.iter_intern;
                else
                    rit = r.iter;
                end
                str = format_number(rit);
                if rit == best_rit
                    str = "\minit{" + str + "}";
                end
                if r.time == best_t
                    str = "\fast{" + str + "}";
                end
                fprintf(f, "%s", str);
            else
                switch r.flag
                    case 1 
                        fl = "m"; % reached maximum number of iterations did not converge to tolerance
                    case 2
                        fl = "c"; % ill conditioned
                    case 3
                        fl = "s"; % stagnated
                    case 4
                        fl = "d"; % quantity too large or too small
                    case 5
                        fl = "i"; % inner iteration did not converge
                end    
                fprintf(f, fl); % stagnated or diverged
            end
        end
        fprintf(f, "\\\\\n");
    end
    fprintf(f, "\\end{tabular}\n");
    fprintf(f, "\n");
end


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

function str = format_number(num)
    if num < 1000
        str = sprintf("%d", num);
    elseif num < 100000
        str = sprintf("%6.1f k", num/1000);
    else
        str = sprintf("%6.1f M", num/1000000);
    end
end
