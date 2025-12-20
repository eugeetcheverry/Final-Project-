function out = kepler2rv(a, ecc, inc, RAAN, w, nu)

%---------------------------- Constants ----------------------------------------------%
mu_earth = 3.986004418e14;

%
%--------------------------------------------------------------------------------------
p = a*(1-ecc ^2);
r_0 = p / (1 + ecc * cos(nu));
%
%%--------------- Coordinates in the perifocal reference system Oxyz -----------------%
%
% position vector coordinates
x = r_0 * cos(nu);
y = r_0 * sin(nu);
%
%
% velocity vector coordinates
Vx_ = -(mu_earth/p)^(1/2) * sin(nu);
Vy_ = (mu_earth/p)^(1/2) * (ecc + cos(nu));
%
%
%%-------------- the geocentric-equatorial reference system OXYZ ---------------------%
%
% position vector components X, Y, and Z
X = (cos(RAAN) * cos(w) - sin(RAAN) * sin(w) * cos(inc)) * x + (-cos(RAAN) * sin(w) - sin(RAAN) * cos(w) * cos(inc)) * y;
Y = (sin(RAAN) * cos(w) + cos(RAAN) * sin(w) * cos(inc)) * x + (-sin(RAAN) * sin(w) + cos(RAAN) * cos(w) * cos(inc)) * y;
Z = (sin(w) * sin(inc)) * x + (cos(w) * sin(inc)) * y;
% velocity vector components X', Y', and Z'
Vx = (cos(RAAN) * cos(w) - sin(RAAN) * sin(w) * cos(inc)) * Vx_ + (-cos(RAAN) * sin(w) - sin(RAAN) * cos(w) * cos(inc)) * Vy_;
Vy = (sin(RAAN) * cos(w) + cos(RAAN) * sin(w) * cos(inc)) * Vx_ + (-sin(RAAN) * sin(w) + cos(RAAN) * cos(w) * cos(inc)) * Vy_;
Vz = (sin(w) * sin(inc)) * Vx_ + (cos(w) * sin(inc)) * Vy_;

out = [X, Y, Z, Vx, Vy, Vz]';