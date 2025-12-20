function [statvec] = kepel_statvec (kepel)
%  [statvec] = kepel_statvec (kepel)
%	The function kepel_statvec transforms the keplerian
%	elements kepel into the corresponding state vector in
%	the same reference system.
%
% Input:
%	kepel
% 		vector containing the keplerian elements:
%		(1) - semimajor axis of the orbit in meters.
%		(2) - eccentricity.
%		(3) - inclination in radians.
%		(4) - right ascension of ascending node in radians.
%		(5) - argument of perigee in radians.
%		(6) - true anomaly in radians.
%
% Output:
%	statvec
%  		state vector in meters and meters/second.
%

global MU_EARTH

a    = kepel(1);	% semi-major axis
e    = kepel(2);	% eccentricity
i    = kepel(3);
RAAN = kepel(4);
w    = kepel(5);
nu   = kepel(6);
p    = a*(1-e^2);

cosnu = cos(nu); sinnu = sin(nu);
rpqw = [p*cosnu/(1+e*cosnu), p*sinnu/(1+e*cosnu), 0]';
vpqw = [-sqrt(MU_EARTH/p)*sinnu, sqrt(MU_EARTH/p)*(e+cosnu), 0]';
cosraan = cos(RAAN); sinraan = sin(RAAN);
cosw = cos(w); sinw = sin(w);
cosi = cos(i); sini = sin(i);
MAT = [cosraan*cosw-sinraan*sinw*cosi, -cosraan*sinw-sinraan*cosw*cosi, sinraan*sini; sinraan*cosw+cosraan*sinw*cosi, -sinraan*sinw+cosraan*cosw*cosi, -cosraan*sini; sinw*sini, cosw*sini, cosi];
r = MAT*rpqw;
v = MAT*vpqw;

statvec = [r; v]';

