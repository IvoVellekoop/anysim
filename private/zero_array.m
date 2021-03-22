function D = zero_array(N, opt)
%ZERO_ARRAY Constructs a zero array of size N with proper formatting (double/single, gpu or not)
%   OPT.precision   = 'single' or 'double'
%   OPT.gpu_enabled = 'true' or 'false'
%
    Dscalar = data_array(0, opt); % avoid code duplication with data_array
    D = zeros(N, 'like', Dscalar);
end

