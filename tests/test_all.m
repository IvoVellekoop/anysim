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
    struct(preconditioner = "shift", preconditioner_shift = 0.5),...
    struct(preconditioner = "circulant"),...
    struct(preconditioner = "moborn"),...
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
    pre = preconditioners{p_i};
    fprintf(f, "\n\\textbf{%s}\\\\\n", pre.preconditioner);
    fprintf(f, "\\begin{tabular}{l|c|c|c|c|c|c|c}");
    results = all_results{p_i, 1};
    % temporary: remove bicgstab(l) column:
    results(4) = [];

    for r_i = 1:numel(results)
        r = results(r_i);
        fprintf(f, "& %s", r.name);
    end
    fprintf(f, "\\\\\n\\hline\n");

    for t_i = 1:size(all_results, 2)
        % write table data
        fprintf(f, "%s ", tests{t_i, 1});
        results = all_results{p_i, t_i};
        % temporary: remove bicgstab(l) column:
        results(4) = [];
        
        % find best
        best_r_i = 0;
        best_rit = inf;
        for r_i = 1:numel(results)
            r = results(r_i);
            rit = max(r.iter_intern, r.iter);
            if r.flag == 0 && rit < best_rit
                best_rit = rit;
                best_r_i = r_i;
            end
        end

        for r_i = 1:numel(results)
            r = results(r_i);
            if r_i == best_r_i
                fprintf(f, " & \\textbf{");
            else
                fprintf(f, " & ");
            end
            if r.flag == 0
                if (r.iter_intern == 0)
                    fprintf(f, "%s", format_number(r.iter));
                else
                    fprintf(f, "%s (%s)", format_number(r.iter_intern), format_number(r.iter));
                end
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
            if r_i == best_r_i
                fprintf(f, "}");
            end
        end
        fprintf(f, "\\\\\n");
    end
    fprintf(f, "\\end{tabular}\n");
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

function str = format_number(num)
    if num < 1000
        str = sprintf("%d", num);
    elseif num < 100000
        str = sprintf("%6.1f k", num/1000);
    else
        str = sprintf("%6.1f M", num/1000000);
    end
end
