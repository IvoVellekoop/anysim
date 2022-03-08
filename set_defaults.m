function opt = set_defaults(defaults, opt)
    arguments
        defaults struct
        opt struct
    end
    
    %SET_DEFAULTS Fills in default values for missing fields
    %   SET_DEFAULTS(DEFAULTS, OPT) Takes structure OPT and adds all fields that
    %   are present in DEFAULTS and missing in OPT.
    %
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

