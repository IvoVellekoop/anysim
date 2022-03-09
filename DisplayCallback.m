classdef DisplayCallback
    %DISPLAYCALLBACK Callback for displaying the a cross section of
    %   the current solution for a GridSim-based simulation
    %   
    %   Don't create an instance of this class directly, instead
    %   specify it as a callback in the options passed to the Sim object.
    %   For example: 
    %     opt.callback = @DisplayCallback
    %     opt.cross_section = @(u) u(3, :, end/2)
    %     sim = DiffuseSim(D, a, opt)
    %
    %   Options:
    %     cross_section    selects the data to display. It holds a cell
    %                      array, with either a relative position where to 
    %                      make the cross section (between 1/end and 1), or []
    %                      to indicate taking all data along that dimension.
    %                      For vector-data, the first element is an integer
    %                      indicating which componetn to display.
    %                      For example, to display the 3th component of a
    %                      4-element vector field, over a cross section at
    %                      x=end/2, and all y and z.
    %                      {4, 1/2, ':', ':'}
    %                       
    %     show_boundaries  when false (default), only shows the region of
    %                      interest. When true, includes the boundaries.
    %     show_convergence when true (default), shows a plot of ‖Δ𝜓‖^2 to
    %                      monitor convergence
    %
    properties (SetAccess = private)
        % Grid object with scaling and cropping information
        grid Grid
        
        % scaling operator/matrix: we want to plot the solution in non-scaled form  
        Tr

        %component % index of the vector component we are displaying
        imageplot % true for 2-dimensional cross sections, false otherwise
        coord1 % coordinates for 'x' axis
        coord2 % coordinates for 'y' axis
        label1 % label for 'x' axis
        label2 % label for 'y' axis
        selection % subsref structure to access cross section to display
    end

    properties
        % for tensor data: 2-element index of component to display.
        % for vector data: 1-element index of component to display.
        % for scalar data: ignored
        % missing values: default to 1
        component { mustBePositive, mustBeInteger } = []

        % select data region to display. For each dimension, either indicate
        % ':' (display all data along that dimension) or a fraction between
        % 0 (first 'row') or 1 (last 'row).
        % missing values: use 0.5 for all other dimensions.
        cross_section (:,1) = {':', ':'}

        % true to draw a graph of magnitude of each update step
        show_convergence logical = true
    end
    
    methods
        function obj = DisplayCallback(opt)
            arguments
                opt.?DisplayCallback
            end
            obj = copy_properties(obj, opt);
        end

        function obj = prepare(obj, sim)
            arguments
                obj
                sim GridSim  % This feedback function can only be used with grid-based simulations!
            end
            obj.grid = sim.grid; 
            obj.Tr = sim.Tr;
            
            % convert component and cross_section to a structure that can
            % be used in subsref
            indices = cell(length(obj.grid.N_u), 1);

            % component selection
            valuedim = length(obj.grid.N_components);
            obj.component = extend(obj.component, ones(valuedim, 1)); % missing values: 1
            for i = 1:valuedim
                indices{i} = obj.component(i);
            end
            
            % roi cross section
            dims = [];
            for i = 1:obj.grid.N_dim
                % missing values default to 0.5
                if (i > length(obj.cross_section))
                    obj.cross_section{i} = 0.5;
                end

                pos = obj.cross_section{i};
                if pos == ':'
                   indices{i  + valuedim} = ':';
                   dims = [dims i]; %#ok<AGROW>
                else
                   N = numel(obj.grid.coordinates(i));
                   indices{i + valuedim} = min(max(round(pos * (N-1) + 1), 1), N);
                end
            end
            obj.selection.type = '()';
            obj.selection.subs = indices;
            
            %% find out in what dimensions the cross section was taken,
            if length(dims) > 2 || length(dims) < 1
                error('DisplayCallback: Cross section must be 1 or 2 dimensional');
            end
        
            %% Prepare labels and coordinates
            if length(dims) == 2 && all(sim.grid.N(dims) > 1)
                obj.imageplot = true;
                obj.label2 = sprintf("y [%s]", obj.grid.pixel_unit(dims(2)));
                obj.coord2 = sim.grid.coordinates(dims(2));
            else
                obj.imageplot = false;
            end
            obj.label1 = sprintf("x [%s]", obj.grid.pixel_unit(dims(1)));
            obj.coord1 = sim.grid.coordinates(dims(1));
        end
        
        function call(obj, u, ~, state)
            if obj.show_convergence && ~isempty(state.residuals)
                subplot(2, 1, 1);
                semilogy(state.residual_its, state.residuals / max(state.residuals));
                xlabel('Iteration');
                ylabel('‖Δ𝜓‖^2 (normalized)'); 
                subplot(2, 1, 2);
            end
            u = obj.grid.crop(u);
            u = real(fieldmultiply(obj.Tr, u));
            u = squeeze(subsref(u, obj.selection));

            if obj.imageplot
                imagesc(obj.coord1, obj.coord2, u.');
                xlabel(obj.label1);
                ylabel(obj.label2);
                title(state.iteration);
                colorbar;
                axis image;
            else
                plot(obj.coord1, u);
                xlabel(obj.label1);
                title(state.iteration);
            end
            drawnow();
        end
    end
end

