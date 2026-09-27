classdef rmm_estimator
    properties
        x
        sigma
        Q
        Q_min
        R
        R_min
        H
        iner
        alpha_q
        alpha_r
        controller
        dt
        
    end
    
    methods
        
        function obj = rmm_estimator(x_init, sigma_init, Q_init, R_init, H, Q_min, R_min, alpha_q, alpha_r, iner, controller, dt)
            obj.x = x_init;
            obj.sigma = sigma_init;
            obj.Q = Q_init;
            obj.R = R_init;
            obj.iner = iner;
            obj.alpha_q = alpha_q;
            obj.alpha_r = alpha_r;
            obj.Q_min = Q_min;
            obj.R_min = R_min;
            obj.H = H;
            obj.dt = dt;
            obj.controller = controller;
            
        end

        function obj = set_init(obj, x_init)
            obj.x = x_init;
        end
        
        function [obj, phikm1] = update(obj, mag_mom, earth_field_b, z)
            
            w = obj.x(1:3);
            rmm = obj.x(4:6);

            xdot_hat = [inv(obj.iner)*(cross(mag_mom, earth_field_b) + cross(rmm, earth_field_b) - cross(w, obj.iner*w));
                        0;
                        0;
                        0];
            
            [~, k_v, eps] = get_parameters(obj.controller);
            normb2 = norm(earth_field_b)^2;
            w_skew = Skew(w);
            Jw = obj.iner*w;
            Jw_skew = Skew(Jw);
            B_skew = Skew(earth_field_b);

            A = [obj.iner\(Jw_skew - w_skew*obj.iner - 1/normb2 *B_skew * eps*k_v*obj.iner* B_skew)  -inv(obj.iner)*B_skew; %A
                    zeros(3) zeros(3)];

            phikm1 = expm(A*obj.dt);
            %Agregar terminos de mayor orden
            %phikm1 = eye(6) + A*dt + 0.5*A^2*dt;

            %Al eliminar la dependencia de u, B podría ser una identidad
            B = eye(6);
            %B = [-inv(iner)*B_skew zeros(3);
            %    zeros(3) eye(3)];
            
            gammkm1 = obj.dt*phikm1*B;
            %gammkm1 = B*dt;
            
            %Pasito de prediccion
            x_hat = [w; rmm] + xdot_hat*obj.dt;
            sigma_hat = phikm1*obj.sigma*phikm1' + gammkm1*obj.Q*gammkm1';

            %Cambiado respecto del paper, se usa la sensibilidad (robotica) en vez
            %de solo la R
            K = sigma_hat*obj.H'/(obj.H*sigma_hat*obj.H' + obj.R);
            dy = z - obj.H*x_hat;

            obj.x = x_hat + K*dy;
            obj.sigma = (eye(6) - K*obj.H)*sigma_hat;
            
            
            Q_avg = K*(dy*dy')*K';
            obj.Q = (1-obj.alpha_q)*obj.Q + (obj.alpha_q)*(Q_avg + obj.Q_min);

            res = z - obj.H*obj.x;
            obj.R = (1-obj.alpha_r)*obj.R + (obj.alpha_r)*(res*res' + obj.H*sigma_hat*obj.H' + obj.R_min);



            %Cosas que no andan:
            %Runge para orden 4 para la integracion del sistema no lineal
            %Hacer moving average para Q_avg
            %Ver los exponentes de Floquet, aunque da que hablar. 


            %Investigar (aunque capaz la rta es un no):
            %Encontrar un paso de tiempo adecuado para lograr un control del ángulo
            %alrededor del nadir

            %Probar:
            %Realizacion de montecarlo
        end
        
        function [x, sigma, Q, R] = get_estimates(obj)
            x = obj.x;
            sigma = obj.sigma;
            Q = obj.Q;
            R = obj.R;
        end
        
    end
    
end