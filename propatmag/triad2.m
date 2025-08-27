function c_triad = triad2(vb, wb, vr, wr)
% [c_triad] = triad(vb, wb, vr, wr)
% purpose:
%   To calculate the rotation matrix (attitude) from two pairs of given 
%   vectors, known in two different reference frames, using the TRIAD 
%   algorithm.
% inputs:
%   vb, wb
%       Body vectors v and w (3x1)
%   vr, wr
%       Reference vectors v and w (3x1)
%
% output
%   c_triad
%       Transformation matrix: u_b = c_triad * u_r
%       c_triad = c_b_r
%
% Valdemir Carrara, Oct, 2012
%
    c_triad = triad(vb, wb)'*triad(vr, wr);
