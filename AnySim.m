classdef (Abstract) AnySim < handle
    %ANYSIM Simulation system for linear operator equations
    %   (c) 2019. Ivo Vellekoop
    %
    % This is the base class for all simulation objects. It implements
    % the split-Richardson algorithm in an abstract way.
    % It is the responsibility of the constructor of the simulation
    % object to implement the operators for the medium,
    % propagator, source, transform and scaling
    %
    properties
        % operators (see readme.md) for a detailed explanation
        medium;     % operator V+1. This object implements two functions:
                    % mix(phi, prop_phi) = (V+1)^2 prop_phi + (V+1) phi
                    % mix_final(phi, prop_phi) = (V+1) prop_phi + phi
                    % also includes fields for preconditioning matrices Tl, V0 and Tr
        propagator; % operator (L-1)^-1
        transform;  % operators (V+1) and (L-1)^-1 are each implemented
                    % in their own domain (typically real-space and 
                    % Fourier-transformed space). The transform object
                    % transforms between both domains.
        opt;        % simulation options
    end
        
    methods
        function obj = AnySim(opt)
            % This function sets the default options for the simulation.
            % Subclasses can override this function to set defaults for
            % options they define. In that case, subclasses should
            % call defaults() on the superclass here
            % (e.g. opt = AnySim.defaults).
            % All classes should call set_defaults(opt, ClassName.defaults)
            % at the start of the constructor.
            
            % flag to determine if simulation are run
            % on the GPU (default: run on GPU if we have one)
            defaults.gpu_enabled = gpuDeviceCount > 0;

            % flag to determine if single precision or 
            % double precision calculations are used.
            % Note that on a typical GPU, double
            % precision calculations are about 10
            % times as slow as single precision.
            defaults.precision = 'single';
            
            % default termination condition (see tc_relative_error)
            defaults.termination_condition.handle = @TerminationCondition;
            defaults.termination_condition.interval = 16;
            defaults.callback.handle = @(u, state, opt) disp(state.iteration);
            defaults.callback.interval = 16;
            obj.opt = set_defaults(defaults, opt);
        end
        
        function [u, state] = exec(obj, source)
            % EXEC executes the simulation for the given source
            % See README.md for a detailed explanation of the algorithm

            %% Initialize state and source
            % The iteration is implemented as:
   % t1 = G u + s            Medium.mix_source
   % t1 -> Li t1             Propagator.propagate
   % u -> u + G (t1 - u)     Medium.mix_field
            [u, state] = obj.start();
            
            while state.running()
                % Gu + s
                t1 = obj.medium.mix_source(u, source, state); 
                
                % Li t1
                t1 = obj.transform.r2k(t1, state);
                t1 = obj.propagator.apply(t1, state);
                t1 = obj.transform.k2r(t1, state);
                
                % in case r domain is non-stationary: transform u_r
                % to the proper domain (does nothing yet)
                u = obj.transform.r2r(u, state);
                
                % u + G (t1-u)
                u = obj.medium.mix_field(u, t1, state);
                state.next(u);
            end
            
            % u -> Tr u (convert to non-scaled solution)
            % + optional final processing (in grid-based simulations the solution
            % is cropped to the roi)
            u = obj.finalize(u, state);
        end
        
        function c = coordinates(obj, d)
            % SIM.COORDINATES(D) returns the coordinates in dimension
            % D, cropped to the region of interest. Coordinates are 0 for
            % the start of the ROI, and increasing.
            c = obj.grid.crop(obj.grid.coordinates(d));
        end
    end
    methods (Abstract, Access=protected)
        [u, state] = start(obj)
        u = finalize(obj, u, state)
    end
end


