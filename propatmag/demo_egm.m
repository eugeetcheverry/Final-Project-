% *************************************************************************
%                            egm_demo
%
%   Program to demonstrate how to use the Earth Gravity Model functions
%   and the numeric orbit propagator
%
%   Valdemir Carrara, july 2017
%
%**************************************************************************

% angle convertion
deg2rad     = pi/180;

% Inputs:
% Orbit keplerian elements: (meters, radians)
kepel = [7000000, 0.01, 98*deg2rad, 0, 35*deg2rad, 0];   

% Orbit state vector:
stat = kepel_statvec(kepel);
    
% Compute the variations in keplerian elements due to the Earth oblateness
delk = delkep(kepel); % just to compare with the numeric orbit propagator

year = 2017;
mjd = djm(17, 7, year);     % Format (day, month, year)

% Ephemerides time:
dfra = time_to_dayf (23, 0, 0);    % UTC time in (hour, minute, sec)

% Propagation time in seconds:
tstart  = 0;        % initial time (sec)
tstep   = 1;        % step time (sec)
tend    = 6000;     % end time (~1 orbit)
n       = fix(tend/tstep);  % data size

% Select the egm model and the polynomial order
%egm_read_data('egm_10.dat', 9) % use this to select both
egm_read_data('egm_10.dat') % use this for maximum polynomial order

% ODE solver precision:
options = odeset('abstol', 1e-4, 'reltol', 1e-4);

% data storage
z1      = zeros([1, n]);
z3      = [z1; z1; z1];
r_time      = z1;
r_xo        = z3;
r_vo        = z3;
r_sma       = z1;
r_ecc       = z1;
r_inc       = z1;
r_raan      = z1;
r_par       = z1;
r_ma        = z1;
r_dist      = z1;
r_rx        = z1;
r_ry        = z1;

% initial input
dist_acc = [0; 0; 0];
cont_acc = [0; 0; 0];

% data index
ic          = 1;

% Orbit propagation
for t = tstart:tstep:tend

    % Analytical orbit propagation
    kp_an   = kepel + delk*t;
    
    % Convert from keplerian elements to state vector
    sv_an   = kepel_statvec(kp_an);
    xi_an   = sv_an(1:3)';
    vi_an   = sv_an(4:6)';

    % Orbit reference frame rotation matrix
    c_i_o   = orbital_to_inertial_matrix(kp_an);
    
    % ODE Solver parameters
    tspan   = [t, t+tstep/2, t+tstep];

    %external force:
    ext_acc   = dist_acc + cont_acc;
    
    % Numeric integration (ODE45)
    [T, Y] = ode45('egm_difeq', tspan, stat, options, mjd, dfra, ext_acc);

    sv_nm   = Y(3, :)';         % propagated state vector
    xi_nm   = sv_nm(1:3);       % propagated inertial posititon vector
	vi_nm   = sv_nm(4:6);       % propagated inertial velocity vector
    stat    = sv_nm;            % state vector update
    
    % numerically propagated keplerian elements
    kp_nm   = statvec_kepel(sv_nm');
    
    % eccentric anomaly
    ea_nm   = kepler(kp_nm(6), kp_nm(2));
    
    % geocentric distance
    dist    = kp_nm(1)*(1 - kp_nm(2)*cos(ea_nm));
        
    % orbit control acceleration (if any)
 	cont_acc = [0; 0; 0];
 
    % disturbance specific forces (if any)
 	dist_acc = [0; 0; 0];

    % Store data to be plotted
    r_time(ic)  = t;
    r_xo(:, ic) = c_i_o'*(xi_nm - xi_an)/1000;
    r_vo(:, ic) = c_i_o'*(vi_nm - vi_an);
    r_dist(ic)  = dist/1000;
    r_rx(ic)    = kp_nm(1)*(cos(ea_nm) - kp_nm(2))/1000;
    r_ry(ic)    = kp_nm(1)*sin(ea_nm)*sqrt(1 - kp_nm(2)^2)/1000;
    r_sma(ic)   = kp_nm(1) - kp_an(1)/1000;
    r_ecc(ic)   = kp_nm(2) - kp_an(2);
    r_inc(ic)   = kp_nm(3) - kp_an(3);
    r_raan(ic)  = kp_nm(4) - kp_an(4);
    r_par(ic)   = kp_nm(5) - kp_an(5);
    r_ma(ic)    = proximus(kp_nm(6), kp_an(6)) - kp_an(6);
    ic          = ic + 1;

end

close all

% Output visualization
plot(r_time, r_xo(1,:), 'r', r_time, r_xo(2,:), 'g', r_time, r_xo(3,:), 'b');
xlabel('Time (s)')
ylabel('Satellite position (km)')
title('Satellite position in orbit frame')

figure
plot(r_time, r_vo(1,:), 'r', r_time, r_vo(2,:), 'g', r_time, r_vo(3,:), 'b');
xlabel('Time (s)')
ylabel('Satellite velocity (m/s)')
title('Satellite velocity in orbit frame')

figure
plot(r_time, r_dist);
xlabel('Time (s)')
ylabel('Distance (km)')
title('Geocentric distance')

figure
plot(r_rx, r_ry);
xlabel('Orbit plane - x (km)')
ylabel('Orbit plane - y (km)')
title('Orbit')

figure
plot(r_xo(2,:), r_xo(1,:));
xlabel('Along track position (km)')
ylabel('Zenith position (km)')
title('Satellite position in orbit plane')

figure
plot(r_xo(3,:), r_xo(1,:));
xlabel('Cross track position (km)')
ylabel('Zenith position (km)')
title('Cross track satellite position')

figure
plot(r_time, r_sma);
xlabel('Time (s)')
ylabel('Relative semi major axis (km)')
title('Semi major axis variation')

figure
plot(r_time, r_ecc);
xlabel('Time (s)')
ylabel('Relative eccentricity')
title('Eccentricity variation')

figure
plot(r_time, r_inc/deg2rad);
xlabel('Time (s)')
ylabel('Relative inclination (deg)')
title('Orbit inclination variation')

figure
plot(r_time, r_raan/deg2rad);
xlabel('Time (s)')
ylabel('Relative right ascention (deg)')
title('Right ascention of ascending node variation')

figure
plot(r_time, r_par/deg2rad);
xlabel('Time (s)')
ylabel('Relative perigee argument (deg)')
title('Perigee argument variation')

figure
plot(r_time, r_ma/deg2rad);
xlabel('Time (s)')
ylabel('Relative mean anomaly (deg)')
title('Mean anomaly variation')



    

