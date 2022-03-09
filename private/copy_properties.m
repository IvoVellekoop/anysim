function target = copy_properties(target, source)
    %COPY_PROPERTIES Copies properties from source to target
    %   If a property is missing in source or target, it is skipped
    props = intersect(fields(source), properties(target));
    for p = 1:length(props)
        target.(props{p}) = source.(props{p});
    end
end

