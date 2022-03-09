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
        % SimGrid object with scaling and cropping information
        grid SimGrid
        
        % scaling operator/matrix: we want to plot the solution in non-scaled form  
        Tr

        %component % index of the vector component we are displaying
        imageplot % true for 2-dimensional cross sections, false otherwise
        coord1 % coordinates for 'x' axis
        coord2 % coordinates for 'y' axis
        label1 % label for 'x' axis
        label2 % label for 'y' axis
        valuedim % 0 for scalar fields, 1 for vector fields, 2 for matrix fields, etc.
        selection % subsref structure to access cross section to display
    end

    properties
        % select data to display
        cross_section = {1}

        % when true, shows all simulation data. when false, removes the boundaries first
        show_boundaries logical = false

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
            
            % Check if 'u' is a vector/matrix field. In that case,
            % by default display the first element of the vector/matrix,
            obj.valuedim = length(obj.grid.N_components);
            
            % convert cross_section to a structure that can be used in
            % subsref
            dims = [];
            indices = cell(length(obj.grid.N_u), 1);

            for i = 1:length(obj.grid.N_u)
                % depending on the data type, the first 0-2 elements of 
                % cross_section are treated as absolute indices
                % (for instance {1} selects the x component in a vector
                % field)
                if (i <= obj.valuedim)
                    if ~isempty(obj.cross_section)
                        component = obj.cross_section{1};
                    else
                        component = 1;
                    end
                    validateattributes(component, "numeric", "scalar");
                    indices{i} = obj.cross_section{1};
                    continue;
                end
                if (i > length(obj.cross_section)) % nothing specified: assume ':' for first two non-singleton dimensions
                    if length(dims) < 2 && obj.grid.N_u(i) > 1
                        indices{i} = ':';
                        dims = [dims i - obj.valuedim]; %#ok<AGROW>
                        continue;
                    else
                        pos = 0.5; % assume center position for all followin dimensions
                    end
                else
                    % specified position
                    pos = obj.cross_section{i};
                    if pos == ':'
                        indices{i} = ':';
                        dims = [dims i]; %#ok<AGROW>
                        continue;
                    end
                end
                if ~obj.show_boundaries
                    bw = obj.grid.boundaries_width(i - obj.valuedim);
                    N = obj.grid.N_roi(i - obj.valuedim);
                    indices{i} = max(round(pos * N), 1) + floor(bw);
                else
                    indices{i} = max(round(pos * obj.grid.N_u(i)), 1);
                end
            end
            obj.selection.type = '()';
            obj.selection.subs = indices;
            
            %% find out in what dimensions the cross section was taken,
            if length(dims) > 2 || length(dims) < 1
                error('DisplayCallback: Cross section must be 1 or 2 dimensional');
            end
            if length(dims) == 2
                obj.imageplot = true;
            else
                obj.imageplot = false;
                dims(2) = dims(1);
            end
            
            %% Prepare labels and coordinates
            obj.label1 = sprintf("x [%s]", obj.grid.pixel_unit(dims(1)));
            obj.label2 = sprintf("y [%s]", obj.grid.pixel_unit(dims(2)));
            obj.coord1 = sim.grid.coordinates(dims(1));
            obj.coord2 = sim.grid.coordinates(dims(2));
            if ~obj.show_boundaries
                obj.coord1 = obj.grid.crop(obj.coord1);
                obj.coord2 = obj.grid.crop(obj.coord2);
            end
        end
        
        function call(obj, u, ~, state)
            if obj.show_convergence && ~isempty(state.diffs)
                subplot(2, 1, 1);
                semilogy(state.diff_its, state.diffs / max(state.diffs));
                xlabel('Iteration');
                ylabel('‖Δ𝜓‖^2 (normalized)'); 
                subplot(2, 1, 2);
            end
            if ~obj.show_boundaries
                u = obj.grid.crop(u, obj.valuedim);
            end
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

