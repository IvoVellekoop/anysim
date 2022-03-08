function sz = compatible_size(sz1, sz2)
    arguments
        sz1 double {mustBeVector} 
        sz2 double {mustBeVector} 
    end
    %COMPATIBLE_SIZE   Returns the size of the array that would result when 
    %                  adding arrays of SZ1 and SZ2.
    %
    %   SZ = COMPATIBLE_SIZE(SZ1, SZ2) is equivalent to size(zeros(sz1) + zeros(sz2))
    %   SZ1 and SZ2 are compatible if all elements are either equal, or if
    %   one of the elements equals 1 (corresponding to singleton expansion).
    %   If the sizes are compatible, returns the maximum size along each
    %   dimension. If they are not compatible, throw an error.
    %   Note: the size vectors are padded with ones to have the same length,
    %   the resulting vector sz may still have trailing ones.
    %   
    
    % pad size vectors to have the same length
    ldif = length(sz1) - length(sz2);
    sz1 = [sz1(:); ones(-ldif, 1)];
    sz2 = [sz2(:); ones(ldif, 1)];
    compatible = sz1 == sz2 | sz1 == 1 | sz2 == 1;
    if ~all(compatible) 
        disp(sz1);
        disp(sz2);
        error('Incompatible sizes');
    end
    sz = max(sz1, sz2);
end

