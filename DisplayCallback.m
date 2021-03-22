classdef DisplayCallback
    %DISPLAYCALLBACK Callback for displaying the current value
    %   of the solution. Options:
    %   	CROSS_SECTION: vector with a number
    %
    %   Example:
    %   % Select this object as callback
    %   opt.callback.handle = @DisplayCallback;
    %
    %   % Specify the cross section to show
    %   % may return 1-D or 2-D data (resulting in line plot or image)
    %   opt.callback.cross_section = [4, 0, 0.5]
    %
    
    properties
        grid % simulation grid
        Tr % scaling matrix: we want to plot the solution in non-scaled form  
        cross_section % handle of function to select data to display
        show_boundaries % when true, shows all simulation data. when false, removes the boundaries first
        component % index of the vector component we are displaying
        imageplot % true for 2-dimensional cross sections, false otherwise
        coord1 % coordinates for 'x' axis
        coord2 % coordinates for 'y' axis
        label1 % label for 'x' axis
        label2 % label for 'y' axis
    end
    
    methods
        function obj = DisplayCallback(sim, opt)
            %DISPLAYCALLBACK Don't call this function directly
            % see example in class documentation.
            default.cross_section = @(u) u;
            default.show_boundaries = false;
            opt = set_defaults(default, opt);
            
            if ~isa(sim, 'GridSim')
                error("CB_IMAGE Feedback function can only be used with grid-based simulations");
            end
            obj.grid = sim.grid;
            obj.cross_section = opt.cross_section;
            obj.show_boundaries = opt.show_boundaries;
            
            % find out which component the cross-section function is returning
            % and along which dimensions the data is cropped.
            % This is a bit of a hack!
            hack = zeros(sim.N) + (1:sim.N(1)).';
            hack = obj.cross_section(hack);
            obj.component = hack(1);
    
            % find out in what dimensions the cross section was taken,
            % and locate appropriate coordinates for plotting
            dims = 0:length(sim.N);
            dims = dims(size(hack) > 1);
            if length(dims) > 2 || length(dims) < 1
                error('DisplayCallback: Cross section must be 1 or 2 dimensional');
            end
            if length(dims) == 2
                obj.imageplot = true;
            else
                obj.imageplot = false;
                dims(2) = dims(1);
            end
            obj.label1 = sprintf("x [%s]", obj.grid.unit(dims(1)));
            obj.label2 = sprintf("y [%s]", obj.grid.unit(dims(2)));
            obj.coord1 = obj.grid.coordinates(dims(1));
            obj.coord2 = obj.grid.coordinates(dims(2));
            if ~obj.show_boundaries
                obj.coord1 = obj.grid.crop(obj.coord1);
                obj.coord2 = obj.grid.crop(obj.coord2);
            end
            
            obj.Tr = sim.medium.Tr(obj.component, obj.component);
            if ~isdiag(sim.medium.Tr)
                error('only works for diagonal matrices!');
            end
        end
        
        function call(obj, u, state)
            if ~obj.show_boundaries
                u = obj.grid.crop(u, 1);
            end
            u = obj.Tr * squeeze(real(obj.cross_section(u)));
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

