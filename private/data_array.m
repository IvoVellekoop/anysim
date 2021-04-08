function D = data_array(data, opt)
%DATA_ARRAY Constructs array with proper formatting (double/single, gpu or not)
%   data_array(X, OPT)      converts data in X to formatting specified in OPT
%   OPT.precision   = 'single' or 'double'
%   OPT.gpu_enabled = 'true' or 'false'
%
    switch opt.precision
    case 'single'
        D = single(data);
    case 'double'
        D = double(data);
    otherwise
        error('Precision should be single or double. %s is not supported', opt.precision);
    end
    if opt.gpu_enabled
        D = gpuArray(D);
    else
        D = gather(D);
    end
end

