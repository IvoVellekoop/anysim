classdef FourierTransform
    %FOURIERTRANSFORM FFT transform engine for AnySim
    %   Transforms from r to k space using a fast Fourier transform
    % 
    %   We need to transform over all dimensions except for the first one
    %   Unfortunately, MATLAB does not allow this, so we transform
    %   over all dimensions separately at the moment.
    %
    properties
        real_signal     % true to indicate that the signal (in r space) is real.
                        % imaginary part is discarded.
    end
    
    methods
        function obj = FourierTransform(opt)
            opt = set_defaults(struct('real_signal', false), opt);
            obj.real_signal = opt.real_signal;
        end
        function u = k2r(obj, u, state)  %#ok<INUSD>
            for d=3:ndims(u)
                if size(u,d) > 1
                    u = ifft(u, [], d);
                end
            end
            if obj.real_signal
                u = real(u);
            end
        end
        function u = r2k(obj, u, state)  %#ok<INUSD,INUSL>
            for d=3:ndims(u)
                if size(u,d) > 1
                    u = fft(u, [], d);
                end
            end
        end
        function u = r2r(obj, u, state)  %#ok<INUSD,INUSL>
        end
    end    
end

