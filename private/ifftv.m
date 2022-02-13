function x = ifftv(x)
%FFTV Vector inverse Fourier transform
%   Performs an inverse Fourier transform over all dimensions, except for the first.
%   This is useful for vector fields, where the first dimension corresponds
%   to the vector elements.
    for d=2:ndims(x)
        if size(x, d) > 1
            x = ifft(x, [], d);
        end
    end
end

