function opt = override_options(opt)
%OVERRIDE_OPTIONS Replaces options by global overrides, if set
% When the global variable 'option_override' exists, replace the
% options set in 'opt' by the overrides given in 'option_override'.
    global option_override;
    if isstruct(option_override)
        fields = fieldnames(option_override);
        for k = 1:numel(fields)
            opt.(fields{k}) = option_override.(fields{k});
        end
    end
end

