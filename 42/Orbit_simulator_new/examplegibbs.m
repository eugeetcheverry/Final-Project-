% ˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜˜
% Example_5_01
% ˜˜˜˜˜˜˜˜˜˜˜˜
%
% This program uses Algorithm 5.1 (Gibbs’ method) and
% Algorithm 4.1 to obtain the orbital elements from the data
% provided in Example 5.1.
%
% deg - factor for converting between degrees and
% radians
% pi - 3.1415926...
% mu - gravitational parameter (kmˆ3/sˆ2)
% r1, r2, r3 - three coplanar geocentric position vectors (km)
% ierr - 0 if r1, r2, r3 are found to be coplanar
% 1 otherwise
% v2 - the velocity corresponding to r2 (km/s)
% coe - orbital elements [h e RA incl w TA a]
% where h = angular momentum (kmˆ2/s)
% e = eccentricity
% RA = right ascension of the ascending
% node (rad)
% incl = orbit inclination (rad)
% w = argument of perigee (rad)
% TA = true anomaly (rad)
% a = semimajor axis (km)
% T - period of elliptic orbit (s)
%
% User M-functions required: gibbs, coe_from_sv
% ------------------------------------------------------------
clear
deg = pi/180;
global mu
%...Input data for Example 5.1:
mu = 398600;
r1 = [-294.32 4265.1 5986.7];
r2 = [-1365.4 3637.6 6346.8];
r3 = [-2940.3 2473.7 6555.8];

r1 = [1.543462818353342  -0.924262600026847   6.845292502564548]*1000;
r2 = [   0.960915110271401  -0.944221349993601   6.948230554228246]*1000;
r3 = [   0.356701631919502  -0.957833373058413   7.003385502390747]*1000;

r1 = [1.543462981553771  -0.924263347785903   6.845294158976667]*1000;
r2 = [0.960917089041437  -0.944221786758880   6.948244023493811]*1000;
r3 = [0.356702416138607  -0.957834632657853   7.003389911221596]*1000;

%...
%...Echo the input data to the command window:
fprintf('---------------------------------------------------')
fprintf('\n Example 5.1: Gibbs Method\n')
fprintf('\n\n Input data:\n')
fprintf('\n Gravitational parameter (kmˆ3/sˆ2) = %g\n', mu)
fprintf('\n r1 (km) = [%g %g %g]', r1(1), r1(2), r1(3))
fprintf('\n r2 (km) = [%g %g %g]', r2(1), r2(2), r2(3))
fprintf('\n r3 (km) = [%g %g %g]', r3(1), r3(2), r3(3))
fprintf('\n\n');
%...Algorithm 5.1:
[v2, ierr] = gibbs(r1, r2, r3);
%...If the vectors r1, r2, r3, are not coplanar, abort:
if ierr == 1
fprintf('\n These vectors are not coplanar.\n\n')
return
end
%...Algorithm 4.1
coe = coe_from_sv(r2,v2);
h = coe(1);
e = coe(2);
RA = coe(3);
incl = coe(4);
w = coe(5);
TA = coe(6);
a = coe(7);
%...Output the results to the command window:
fprintf(' Solution:')
fprintf('\n');
fprintf('\n v2 (km/s) = [%g %g %g]', v2(1), v2(2), v2(3))
fprintf('\n\n Orbital elements:');
fprintf('\n Angular momentum (kmˆ2/s) = %g', h)
fprintf('\n Eccentricity = %g', e)
fprintf('\n Inclination (deg) = %g', incl/deg)
fprintf('\n RA of ascending node (deg) = %g', RA/deg)
fprintf('\n Argument of perigee (deg) = %g', w/deg)
fprintf('\n True anomaly (deg) = %g', TA/deg)

kepler = [a e incl/deg RA/deg w/deg TA/deg]
