function [gra_torq] = gravity_torque (r_sat, c_b_i, inert)
    % [gra_torq] = gravity_torque(r_sat, c_b_i, inert)
    %   Function to evaluate the gravity gradient torque from satellite
    %   attitude.
    % inputs:
    %   r_sat
    %       Satellite orbit position in inertial frame (m) (kepel_statvec)
    %   c_b_i
    %      Satellite attitude matrix, with respect to the inertial frame.
    %      v_body = c_b_i * v_inertial
    %   inert
    %      Inertia matrix in kg m^2, in satellite coordinates
    %
    % outputs:
    %   gra_torq 
    %      Gravity gradient torque (spacecraft reference coordinates)
    %      in Nm
    %

    % Valdemir Carrara, Dec 2016

    r_nor       = norm(r_sat);
    vert        = c_b_i*(r_sat/r_nor);
    ggamp       = 3*3.9860064e14/r_nor/r_nor/r_nor;
    gra_torq    = gg_torque (ggamp, inert, vert);

return