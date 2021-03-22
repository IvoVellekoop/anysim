classdef TensorMedium
    %TENSORMEDIUM Helper class to implement a tensor field potential
    %   Detailed explanation goes here
    %
    %   (c) 2019. Ivo Vellekoop
    properties
        Vp % potential array = V + 1 = Tl*(Vraw-V0)*Tr + 1
        V0 % background potential
        Tl % left scaling matrix
        Tr % right scaling matrix
    end
    methods
        function obj = TensorMedium(V_raw, V_raw_min, V_max)
            error("Not implemented completely yet")
            Nc = size(V_raw, 1);
            if size(V_raw, 2) ~= Nc
                error('Tensors in scattering potential must be square')
            end

            % Calculate scaling and centering matrices to minimize norm of V
        
            %% 1. Solve smallest circle problem for eV0h element separately
            centers = zeros(Nc, Nc);
            radii = zeros(Nc, Nc);
            for r=1:Nc
                for c=1:Nc
                    [centers(r, c), radii(r, c)] = smallest_circle(V_raw(r, c, :));
                end
            end
            
            %% 2. Equilibrate matrix so that all elements have magnitude <= 1
            % note: this is not necessarily optimal, but it is a simple
            % way to (more or less) evenly scale the variance in all dimensions
            radii = max(radii, V_min);
            [P,R,C] = equilibrate(radii);
            obj.Tl = P'*R*P;
            obj.Tr = C;
            obj.V0 = centers;
            
            %% 3. Scale to have operator norm 1
            % Apply centering and equilibration scaling
            V = pagemtimes(pagemtimes(obj.Tl, V_raw - obj.V0), obj.Tr);

            % finally scaling step to have ||V|| < 0.95
            Vnorm = max(pagefun(@norm, V), [], 'all');
            V = V * 0.95 / Vnorm; %todo: make 'slack' factor of 0.95 customizable
            obj.Tl = obj.Tl * sqrt(Vnorm);
            obj.Tr = obj.Tr * sqrt(Vnorm);
            obj.Vp = data_array([], opt, V + eye(Nc)); % V+1
        end
                
        function phi = mix(phi, prop_phi)
            % phi => (V+1)^2 prop_phi + (V+1) phi
            phi = phi + pagemtimes(obj.Vp, prop_phi);
            phi = pagemtimes(obj.Vp, phi);
        end
        function phi = mix_final(phi, prop_phi)
            % phi => (V+1) prop_phi + phi
            phi = phi + pagemtimes(obj.Vp, prop_phi);
        end
    end
end

