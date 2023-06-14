%how to use. create obj, object = NLS_stochastic_solver,
% then object.result = solve(object)
% then plotting(object)

classdef NLS_stochastic_solver
    properties
        h = 0.00001;
        no_timesteps = 1000000;
        delta = 0.0001;
        result = [];
        random = false;
    end
    methods
        function obj = NLS_stochastic_solver(h, no_timesteps, delta, varargin)
            obj.h = h;
            obj.no_timesteps = no_timesteps;
            obj.delta = delta;
            if length(varargin) == 2
                obj.random = [varargin{1} , varargin{2}];
            end
        end
        function y = solve2(obj)
            y = zeros(1,1,obj.no_timesteps);
        end
        function y = solve(obj)
            h = obj.h;
            sqrt_h = h^0.5;
            no_timesteps = obj.no_timesteps;
            no_timesteps_in_a_period = ceil(2*pi/h);
            My_approximation_for_2pi = no_timesteps_in_a_period*h;
            delta = obj.delta;
            no_timesteps_to_sim_Bloop = ceil(no_timesteps/no_timesteps_in_a_period)*no_timesteps_in_a_period;



            %here is the standard code for a sample path
            epsilon_array = zeros(1,1,no_timesteps);
            eta_array = zeros(1,1,no_timesteps);
            kappa_array = zeros(1,1,no_timesteps);
            f_1_array = zeros(1,1,no_timesteps);
            f_2_array = zeros(1,1,no_timesteps);
            f_3_array = zeros(1,1,no_timesteps);

            if obj.random
                rng(obj.random(1))
            end 
            epsilon = normrnd(0, 1, [1, no_timesteps_to_sim_Bloop]);
            if obj.random
                rng(obj.random(2))
            end 
            eta = normrnd(0, 1, [1, no_timesteps_to_sim_Bloop]);

            %create the sequence of random matrices Sn
            sequence = zeros(3, 3, no_timesteps);

            %create Sn

            %dX = V_1 X kappa dt + V_2 X dtau (can be extended)

            %generate P and Q, both Brownian loop.
            W_1 = sqrt_h*cumsum(epsilon);
            W_2 = sqrt_h*cumsum(eta);
            time = h*cumsum(ones(1,no_timesteps));


            %TODO make this work, need to choose different
            for i = 1:ceil(no_timesteps/no_timesteps_in_a_period);
                lower = (i-1)*no_timesteps_in_a_period + 1;
                upper = min((i*no_timesteps_in_a_period), no_timesteps)
                P(1,((i-1)*no_timesteps_in_a_period + 1):min((i*no_timesteps_in_a_period), no_timesteps)) = W_1(1,((i-1)*no_timesteps_in_a_period + 1):min((i*no_timesteps_in_a_period), no_timesteps)) + (-1)* W_1(no_timesteps_in_a_period*i)*time(1, lower:upper)/My_approximation_for_2pi;
                Q(1, lower:upper) = W_2(1, lower:upper) + (-1)* W_2(no_timesteps_in_a_period)*time(1, lower:upper)/My_approximation_for_2pi;
            end
            %kappa (we dont need the epsilon)
            %kappa = (P.^2 + Q.^2 + delta^2).^(0.5);clear a
            kappa = (P.^2 + Q.^2).^(0.5);

            % tau = f_1 dP + f_2 dQ + f_3 dt
            f_1 = (delta^2 - P.^2).*Q ./((delta^2 + P.^2).^2 + P.^2 .* Q.^2);
            f_2 = P.*(delta^2 + P.^2) ./ ((delta^2 + P.^2).^2 + P.^2 .* Q.^2);

            f_3 = W_1(no_timesteps)*time/My_approximation_for_2pi.*f_1 + W_2(no_timesteps)*time/My_approximation_for_2pi.*f_2 + ((-2)*(P.^3.*Q.*(delta^2 + P.^2) + P.*Q.*((delta^2 + P.^2).^2 + P.^2 .* Q.^2)) - ( delta^2 - P.^2).*Q.*(P.*Q.^2 + 2*P.*(delta^2 + P.^2))) ./ ((delta^2 + P.^2).^2 + P.^2 .* Q.^2).^2;


            %V_1
            V_1 = [ 0 1 0; -1 0 0; 0 0 0];

            %V_2
            V_2 = [ 0 0 0; 0 0 1; 0 -1 0];




            epsilon_array(1,1,:) = epsilon(1, 1:no_timesteps);
            eta_array(1,1,:) = eta(1, 1:no_timesteps);


            kappa_array(1,1,:) = kappa;

            %f_1 = permute(f_1, [1 3 2]);
            %f_2 = permute(f_2, [1 3 2]);
            %f_3 = permute(f_3, [1 3 2]);


            f_1_array(1,1,:) = f_1;
            f_2_array(1,1,:) = f_2;
            f_3_array(1,1,:) = f_3;

            identity_array = zeros(3,3,no_timesteps) + [1 0 0; 0 1 0; 0 0 1];



            %Omega, can see details in Overleaf doc.
            Omega = (V_1.*kappa_array + V_2.*f_3_array)*h + V_2*sqrt_h.*(epsilon_array.*f_1_array + eta_array.*f_2_array) ;


            %now use rodriguez formula to approx exp
            expOmega = zeros(3,3,no_timesteps);

            theta = sqrt(Omega(1,2,:).^2 + Omega(3,3,:).^2 + Omega(2,3,:).^2);
            expOmega = identity_array + sinc( theta / pi ) .* Omega + 0.5* sinc(theta / (2*pi)).^2 .* pagemtimes(Omega , Omega);




            %intitialise X_0 as the identity matrix
            X = identity_array;

            %calculate solution, Nth term is exp(Omega_n)*exp(Omega_{n-1})*...etc
            X(:,:,1) = expOmega(:, :, 1);
            for n = 2:no_timesteps
                X(:, :, n) = expOmega(:,:, n)*X(:,:, n-1);
            end



            %with initial condition y_0 = [0 0 1]^T create a solution path for each
            % X.

            y = zeros(3,1, no_timesteps);
            y(3,1,:) = ones(no_timesteps, 1);
            y = pagemtimes(X, y);
        end
        function a_plot = plotting(obj)
            y = obj.result;
            xcoord = permute(y(1,1,:), [3 2 1]);
            ycoord = permute(y(2,1,:), [3 2 1]);
            zcoord = permute(y(3,1,:), [3 2 1]);


            a_plot = scatter3(xcoord, ycoord, zcoord,10, '*')
            xlabel("x")
            ylabel("y")
            zlabel("z")
            xlim([-1, 1])
        end
    end
end
