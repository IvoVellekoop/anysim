function V = extend(V, defaults)
    %EXTEND Extends vector V with default values
    % The returned vector has the same length as DEFAULTS,
    % where any missing values are copied from DEFAULTS
    % gives an error if length(V) > length(defaults) 
    V = [V defaults(length(V)+1:length(defaults))];
end

