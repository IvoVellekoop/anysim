classdef AnySimOptions
    %ANYSIMOPTIONS Options shared across all ANYSIM simulations
    %   (c) 2022. Ivo Vellekoop
    %
    properties
        % flag to determine if simulation are run
        % on the GPU (default: run on GPU if we have one)
        gpu_enabled (1,1) logical = gpuDeviceCount > 0;
        gpu_device (1,1) {mustBeInteger} = 1;

        % flag to determine if single precision or 
        % double precision calculations are used.
        % Note that on a typical GPU, double
        % precision calculations are about 10
        % times as slow as single precision.
        precision string {mustBeTextScalar, mustBeMember(precision, ["single", "double"])} = "single";

        % When set to 'true', the obj.operator
        % property will hold the scaled  forward operator
        % without preconditioning: H = L+V.
        % Defaults to 'false' because this operator
        % is not needed for regular use and may be expensive
        % to generate.
        forward_operator (1,1) logical = false;

        % value of ‖V‖ after scaling. 0.95 is a compromise, the best
        % value depends on the size of the simulation.
        V_max (1,1) double = 0.95
        
        % value of α used in the Richardson iteration in the basic
        % AnySim algorithm. Optimum is 1 for positive definite systems
        % and between 0.5 and 1 for general accretive systems.
        alpha (1,1) double = 0.75

        % termination condition and callback
        % to specify a different callback, either set the
        % opt.callback.call function directly (as below), or
        % set the opt.callback.handle function to a constructor
        % for a callback object (see DisplayCallBack)
        termination_condition = TerminationCondition();
        termination_condition_interval = 16;
        callback = TextCallback();
        callback_interval = 25;
    end
    methods
        function obj = validate(obj)
            if (obj.gpu_enabled)
                if (gpuDeviceCount > 0)
                    gpu = gpuDevice(obj.gpu_device); % select GPU device
                    disp(['GPU found. Performing simulations on: ', gpu.Name]);
                else
                    warning('No GPU found');
                    obj.gpu_enabled = false;
                end
            end
        end
    end
end
