function [quaternions] = ezyxquat (euler_angles)
% [quaternions] = ezyxquat (euler_angles)
%   To transform Euler angles from a Z-Y-X (yaw-pitch-roll) rotation
%   into quaternions.
% inputs:
%   euler_angles
%       Euler angles from a Z-Y-X rotation, in radians (3)
% outputs:
%   quaternions
%       quaternions corresponding to the Euler angles 
%

rot_mat = ezyxrmx(euler_angles);
      
quaternions = rmxquat(rot_mat);
