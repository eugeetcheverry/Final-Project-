function [] = egm_read_data(egm_data_file, nmax)
% [] = egm_read_data(egm_data_file, nmax)
%
% Function to read and store in memory an Earth Gravity Model (egm) file.
% An error message will be shown if the data file coudn't be
% opened, or the file couldn't be found or data is readen with error
%
% inputs:
%   egm_data_file
%       data file name (string). Current version has:
%           egm_96.dat      96 model of 360th order
%           egm_10.dat      10th model of egm of 30th order
% outputs:
%   [none]
%
% The egm data file shall contain (ASCII, blank separator):
% first register:
%   egm_order
%       Maximum order number of the egm model (integer)
%   egm_lenght
%       Lenght in records of the egm model starting in the second register
%       = (egm_order+1)*(egm_order+2)/2 - 3 (integer)
%   egm_scale_factor
%       Convertion factor (scale factor) of the coefficients, including
%       Earth radius correction and unit convertion (float). It should be 
%       computed from:
%       egm_scale_factor = convertion_factor*egm_earth_radius/6378136.3,
%       where convertion_factor is the coefficients scale, and
%       egm_earth_radius is the standard Earth radius of the egm model
% following registers (each one):
%   n
%       first (zonal) indice, starting in 2 up to egm_order (integer)
%   m   
%       second (tesseral) indice, starting in 0 up to egm_order (integer)
%   cc
%       cosine coeficients of the tesseral harmonics (Cnm) (float)
%   sc
%       sine coeficients of the tesseral harmonics (Snm) (float)
%   ce
%       Cnm variances (not used) (float)
%   se
%       Snm variances (not used) (float)
%
% example (egm_96.dat):
%     360          65338     1.000000000000e+00 
%   2   0 -4.84165371736e-04  0.00000000000e+00  3.5610635e-11  0.0000000e+00 
%   2   1 -1.86987635955e-10  1.19528012031e-09  1.0000000e-30  1.0000000e-30 
%   2   2  2.43914352398e-06 -1.40016683654e-06  5.3739154e-11  5.4353269e-11 
%   3   0  9.57254173792e-07  0.00000000000e+00  1.8094237e-11  0.0000000e+00 
%   3   1  2.02998882184e-06  2.48513158716e-07  1.3965165e-10  1.3645882e-10 
%   3   2  9.04627768605e-07 -6.19025944205e-07  1.0962329e-10  1.1182866e-10 
%   3   3  7.21072657057e-07  1.41435626958e-06  9.5156281e-11  9.3285090e-11
%   ....
% 360 359  1.83971631467e-11 -3.10123632209e-11  5.0033977e-11  5.0033977e-11 
% 360 360 -4.47516389678e-25 -8.30224945525e-11  5.0033977e-11  5.0033977e-11 
%
% author:
%   Valdemir Carrara, July 2017
%

global egm_order egm_length egm_conv_f egm_cc egm_sc
global egm_ae egm_gm egm_pn egm_qn egm_ip egm_nmax

[funit, message] = fopen(egm_data_file, 'r');
if funit < 0
    disp(message)
    stop
end

egm_order = fscanf(funit, '%d', 1);
egm_length  = fscanf(funit, '%d', 1);
egm_length  = (egm_order+2)*(egm_order+1)/2 - 3;
egm_conv_f  = fscanf(funit, '%e', 1);

cf = fscanf (funit, '%g', [6, egm_length]);

egm_cc   = [0 0 0 cf(3, :)];
egm_sc   = [0 0 0 cf(4, :)];

fclose (funit);

egm_ae     = 6378136.3;	% m
egm_gm     = 3986004.415e+8; % m**3/s**2
egm_pn     = zeros(1, egm_order + 1);
egm_qn     = egm_pn;
egm_ip     = egm_pn;
egm_ip(1)  = 0;

for n = 2: egm_order + 1
    egm_ip(n)  = egm_ip(n-1) + n - 1;
end

if (nargin > 1)
    if (nmax < egm_order)
        egm_nmax = nmax;
    else
        egm_nmax = egm_order;
    end
else
    egm_nmax = egm_order;
end

return;
