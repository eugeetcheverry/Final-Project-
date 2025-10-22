classdef controller
    properties
        k_p
        k_v
        eps
        iner
    end
    methods
        function obj = controller(k_p, k_v, eps, iner)
            obj.k_p = k_p;
            obj.k_v = k_v;
            obj.eps = eps;
            obj.iner = iner; 
        end
        
        function [k_p, k_v, eps] = get_parameters(obj)
            k_p = obj.k_p;
            k_v = obj.k_v;
            eps = obj.eps;
        end
        
        function obj = k_v_update(obj, k_v)
            obj.k_v = k_v;
        end
        function u = get_control_action(obj, dq, signq4, dw, pointing)
            
            u = - obj.eps*obj.k_v*obj.iner*dw;
            if pointing
                u = u - obj.eps*obj.eps*obj.k_p*inv(obj.iner)*dq*signq4;
            end
        end
    end
end