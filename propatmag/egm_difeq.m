function dxdt = egm_difeq (time, x, flag, mjd, dsec, ext_acc)
% dxdt = egm_difeq (time, x, flag, ext_forces)
%
%   The egm_difeq function returns with the time derivative of the 
%   satellite orbit position and velocity. The gravity model is internally 
%   computed by the egm functions. However, the egm model has to be
%   initialized by calling the egm_read_dat function. See the egm_read_data
%   function for description.
%
% inputs:
%   time
%       Propagation time (s)
%   x
%       Orbit state vector: position (m) (1:3) and velocity (m/s) (4:6)
%   flag
%       See ODEFILE
%   mjd
%       Modified Julian date referred to 1950.
%   dsec
%       Day time UTC in seconds of initial orbit ephemeris. Current time
%       is (dsec + time)
%   ext_acc
%       External specific forces (acceleration) acting on the spacecraft 
%       (m/s2) including control (thrusters, orbit correction and 
%       disturbances), except gravity.
% outputs:
%   dxdt 
%       Time derivative (velocity and acceleration) of the state vector
% 
% author:
%   Valdemir Carrara,   July, 2017.

xip     = x(4:6);
gwst    = gst(mjd, dsec + time);
se      = inertial_to_terrestrial(gwst, x');
xe      = se(1:3);
ae      = [egm_acc(xe); 0; 0; 0];
ai      = terrestrial_to_inertial(gwst, ae');
vip     = ext_acc + ai(1:3)';

dxdt = [xip; vip];
