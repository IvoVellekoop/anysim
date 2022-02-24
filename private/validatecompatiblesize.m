function validatecompatiblesize(sz1, sz2)
%VALIDATECOMPATIBLESIZE(SZ1, SZ2) Checks if SZ1 and SZ2 are sizes of arrays that can be added together
%   Sizes are compatible if for all elements: the values are equal, or one of the elements is 1
%   The size vectors are padded with ones to have the same length

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
end

