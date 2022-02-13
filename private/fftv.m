function x = fftv(x)
%FFTV Vector Fourier transform
%   Performs a Fourier transform over all dimensions, except for the first.
%   This is useful for vector fields, where the first dimension corresponds
%   to the vector elements.
    for d=2:ndims(x)
        if size(x, d) > 1
            x = fft(x, [], d);
        end
    end
end

