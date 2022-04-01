function x = fftm(x)
%FFTV Matrix field Fourier transform
%   Performs an inverse Fourier transform over all dimensions, except for the first two.
%   This is useful for matrix fields, where the two dimensions correspond
%   to the elements of the matrix.
    for d=3:ndims(x)
        if size(x, d) > 1
            x = fft(x, [], d);
        end
    end
end

