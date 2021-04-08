function opt = set_defaults(defaults, opt)
%SET_DEFAULTS Takes structure opt and fills in missing fields based on
%'defaults' structure.
    if nargin == 1
        opt = defaults;
    else
        f = fieldnames(defaults);
        for i = 1:length(f)
            if ~isfield(opt, f{i})
                opt.(f{i}) = defaults.(f{i});
            elseif isstruct(opt.(f{i}))
                % recursively call set_defaults on structures
                opt.(f{i}) = set_defaults(defaults.(f{i}), opt.(f{i}));
            end
        end
    end
end

