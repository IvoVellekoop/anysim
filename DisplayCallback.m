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
    %     cross_section    handle of function that crops 'u' to select the 
    %                      data to display
    %                      note: at the moment, cross_section can only crop
    %                      'u' and should not perform any other processing
    %                      such as scaling.
    %     show_boundaries  when false (default), only shows the region of
    %                      interest. When true, includes the boundaries.
    %     show_convergence when true (default), shows a plot of ‖Δ𝜓‖^2 to
    %                      monitor convergence
    %
    
    properties
        sim % simulation object (used for cropping etc.)
        Tr % scaling matrix: we want to plot the solution in non-scaled form  
        cross_section % handle of function to select data to display
        show_boundaries % when true, shows all simulation data. when false, removes the boundaries first
        component % index of the vector component we are displaying
        imageplot % true for 2-dimensional cross sections, false otherwise
        coord1 % coordinates for 'x' axis
        coord2 % coordinates for 'y' axis
        label1 % label for 'x' axis
        label2 % label for 'y' axis
        show_convergence % true to draw a graph of magnitude of each update step
    end
    
    methods
        function obj = DisplayCallback(sim, opt)
            %DISPLAYCALLBACK Don't call this function directly
            % see example in class documentation.
            if ~isa(sim, 'GridSim')
                error("CB_IMAGE Feedback function can only be used with grid-based simulations");
            end
            
            default.cross_section = @(u) u;
            default.show_boundaries = false;
            default.show_convergence = true;
            opt = set_defaults(default, opt);
            
            obj.sim = sim;
            obj.cross_section = opt.cross_section;
            obj.show_boundaries = opt.show_boundaries;
            obj.show_convergence = opt.show_convergence;
            
            %% find out which component the cross-section function is returning
            % and along which dimensions the data is cropped.
            % This is a bit of a hack!
            hack = sim.to_external(zeros(sim.N) + (1:sim.N(1)).');
            hack = obj.cross_section(hack);
            obj.component = hack(1);
    
            %% find out in what dimensions the cross section was taken,
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
            
            %% Prepare labels and coordinates
            obj.label1 = sprintf("x [%s]", sim.grid.unit(dims(1)));
            obj.label2 = sprintf("y [%s]", sim.grid.unit(dims(2)));
            obj.coord1 = sim.grid.coordinates(dims(1));
            obj.coord2 = sim.grid.coordinates(dims(2));
            if ~obj.show_boundaries
                obj.coord1 = sim.grid.crop(obj.coord1);
                obj.coord2 = sim.grid.crop(obj.coord2);
            end
            
            obj.Tr = sim.medium.Tr(obj.component, obj.component);
            if ~isdiag(sim.medium.Tr)
                error('only works for diagonal matrices!');
            end
        end
        
        function call(obj, u, state)
            if obj.show_convergence && ~isempty(state.diffs)
                subplot(2, 1, 1);
                semilogy(state.diffs / max(state.diffs(:)));
                xlabel('Interation');
                ylabel('‖Δ𝜓‖^2 (normalized)'); 
                subplot(2, 1, 2);
            end
            if ~obj.show_boundaries
                u = obj.sim.grid.crop(u, 2);
            end
            u = obj.Tr * squeeze(real(obj.cross_section(obj.sim.to_external(u))));
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

