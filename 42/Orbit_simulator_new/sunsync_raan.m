function raan = sunsync_raan (eq_cross_time, gst0)
%
%   This subroutine obtains the orbital right ascension of the ascending 
%   node, given the ascending equator crossing hour and the Greenwich 
%   sidereal time.
% 
% Inputs:
%   eq_cross_time
%       Equator crossing hour of the ascending node, in radians.
%   gst0
%       Greenwich sidereal time at 0:00:00 h UTC
%
% Outputs:
%   raan
%       Right ascension of the ascending node, in rad.
% 
%   Valdemir Carrara         June 2015      V 1.0
% 
      raan      = eq_cross_time + gst0;

