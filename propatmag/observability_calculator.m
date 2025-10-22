classdef observability_calculator
    properties
        H
        Phiks
        N
    end
    methods
        function obj = observability_calculator(H, N)
            obj.H = H;
            obj.N = N;
            obj.Phiks = {};
        end
        
        function obj = update(obj, Phik)
            if size(obj.Phiks,2) < (obj.N - 2)
                obj.Phiks{end+1} = Phik;
            else
                obj.Phiks = obj.Phiks(2:end);
                obj.Phiks{end+1} = Phik;
            end
        end
        
        function O = get_O(obj)
            if size(obj.Phiks,2) == (obj.N-2)
                O = obj.H;
                newmat = obj.H;
                for i = 1:(obj.N-2)
                    newmat = newmat*obj.Phiks{i};
                    O = [O ; newmat];
                end
            else
                O = zeros(size(obj.H, 1));
            end
        end
        
        function M = get_M(obj)
            O = get_O(obj);
            M = O'*O;
        end
    end
end