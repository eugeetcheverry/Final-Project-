function [orbit_va] = orbit_prop(t, orbit_rv, M, F_control, q)%, d)
% Orbit_Prop propagates the satellite's orbit
%   [orbit_va] = orbit_prop(t, orbit_rv, F_control, M) takes in time (unused but
%   needed for ode45), ECI position and speed in [m] and [m/s] as a vector, mass
%   in [kg], control force in body coordinates in [N], and attitude quaternion and
%   returns ECI velocity and acceleration in [m/s] and [m/s2]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = orbit_rv(1:3);
v = orbit_rv(4:6);

F_g = gravity(r,M);
% [gx,gy,gz] = gravitysphericalharmonic(r', 'EGM2008', 20);
% F_g = M*[gx;gy;gz];
% F_d = -v * 0.5*2.2*1.56E-13*norm(v)*0.56 * d;

F_total = F_g + quatrmx(q)*F_control;% + F_d;
a = F_total/M;

orbit_va = [v; a];

% agregar presion solar