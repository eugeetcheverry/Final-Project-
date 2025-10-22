classdef controller
    properties
        k_p
        k_v
        eps
        iner
    end
    methods
        function obj = controller(kp_ini, kv_ini, eps_ini, iner_ini)
            obj.k_p = kp_ini;
            obj.k_v = kv_ini;
            obj.eps = eps_ini;
            obj.iner = iner_ini;
        end
        
        function [k_p, k_v, eps] = get_parameters(obj)
            k_p = obj.k_p;
            k_v = obj.k_v;
            eps = obj.eps;
        end
        
        function obj = update(kp_new, kv_new, eps_new, obj)
            obj.k_p = kp_new;
            obj.k_v = kv_new;
            obj.eps = eps_new;
        end
        
        function u = get_control_action(dq, pointing, signq4, dw, obj)
            u = - obj.eps*obj.k_v*obj.iner*dw;
            if pointing
                u = u - obj.eps*obj.eps*obj.k_p*inv(obj.iner)*dq*signq4;
            end   
        end
    end
end