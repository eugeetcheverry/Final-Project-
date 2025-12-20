function [Tc] = fTc(mean_elem_chief)

    Tc      = eye(6);
    Tc(1,1) = 1/mean_elem_chief(1);
    Tc(2,6) = cos(mean_elem_chief(3));
    Tc(6,6) = sin(mean_elem_chief(3));
