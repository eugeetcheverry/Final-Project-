function M_solar = solar_pressure(v_sun, T, S, r_cp, C_r)
    %Solar pressure:
    %  I  = sigma*T^4*Sup*Coef_reflective
    %     --------------
    %     4*pi*d_sol
    %sigma -> constante de boltzmann
    %r_cp: centro de masa -> centro de presion
    %T: temperatura Sol
    %S: Superficie
    I = ((5.67e-8)*T^4) * 4*pi*(695.7e12)*S*C_r / (4*pi*(1.5e11)^2);
    v_norm = norm(v_sun);
    v_s = v_sun/v_norm;
    f_sun = -(v_s*I);
    M_solar = cross(r_cp, f_sun);
end