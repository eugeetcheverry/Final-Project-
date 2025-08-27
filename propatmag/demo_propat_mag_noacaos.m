% demo_propat
%
%       Program to show how to use PROPAT
%       Type demo_propat from the Matlab prompt

% Inputs:
% Orbit keplerian elements:
kepel = [7000000, 0.01, 98, 0, 0, 0];   % see function delkep
% Orbit state vector:
stat = kepel_statvec(kepel);

% Velocidad Orbital
GM = 3.986004418e14;
orb_period = sqrt(4*pi^2*kepel(1)^3/GM);
omeg_0 = 2*pi/orb_period;

dlyap = 0.0001;

% Attitude elements in Euler angles of a 3-1-3 (z-x-z) rotation
eulzxz = [30, 50, 20]'*pi/180;   % converted from degrees to radians

% Attitude in quaternions
quat = ezxzquat(eulzxz) + sqrt(2)*dlyap*[1;-1;1;0]/sqrt(3);        % converted from Euler angles
quat = quat/norm(quat);

% Angular velocity vector in body frame:
w_ang = 0.05*[0.1, 0, 0.5]' + sqrt(2)*dlyap*[1;-1;1]/sqrt(3);           % in radians/sec
    
% Compute the variations in keplerian elements due to the Earth oblateness
delk = delkep(kepel);

% Ephemerides date in Modified Julian date:
year = 2009;
mjd = djm(13, 4, year);     % Format (day, month, year)
mjdo = djm(1, 1, year);     % modified julian date of 1/1/year
mjd1 = djm(1, 1, year+1);   % modified julian date of 1/1/(year+1)
year_frac = year + (mjd - mjdo)/(mjd1 - mjdo);  % year and fraction

% Ephemerides time:
dfra = time_to_dayf (10, 20, 0);    % UTC time in (hour, minute, sec)

% Propagation time in seconds:
tstart = 0;              % initial time (sec)
tstep = 5;               % step time (sec)
tend = 10*orb_period;    % end time (10 minutes)

% Inertia matrix of axis-symetric rigid body:
iner = [8 0 0; 0 10 0; 0 0 12];         % in kg*m*m
%iner = [27 0 0; 0 17 0; 0 0 25];       % in kg*m*m

% Inverse inertia matrix:
invin = inv(iner);

% Initial control torque:
contq = [0 0 0]';

% Magnetic moment torque flag and moment:
flag_mag = 1;   % 1=compute magnetic moment / 0=discard magnetic moment
mag_mom = [0; 0; 0.1];      % in A.m

% ODE solver precision:
options = odeset('abstol', 1e-4, 'reltol', 1e-4);

% Initial vectors
time = tstart;          % to store time
euler = eulzxz*180/pi;  % Euler angles
omeg = w_ang;           % Angular velocity
orbit = stat';          % Orbit elements (state vector)
keorb = kepel';         % Orbit elements (keplerian)

gamavg = 0*eye(3);  % Gamma avg inicial
maxmagmom = 20;      % Tope momento magnético en Am2
vdq = [0;0;0];
vdw = [0;0;0];
u = [0;0;0];
vu = u;
signq4 = 0;
vext_torq = u;
vmag_mom = u;

% Attitude and orbit propagation
for t = tstart:tstep:tend

    % Orbit propagation
    kep2 = kepel + delk*t;
    
    % To convert from keplerian elements to state vector (if needed)
    stat = kepel_statvec(kep2);

    % Perturbation torques:
    ambt = [0 0 0]';
    A    = quatrmx(quat);
    exb  = A(:,1);
    ambt = 3*omeg_0^2*cross(exb,iner*exb);
    
    % External torques (perturbation + control)
    ext_torq = ambt + contq;
    
    % Initial attitude vector:
    att_vec = [quat; w_ang]';         % Transposed

    % ODE Solver parameters
    tspan = [t, t+tstep/2, t+tstep];
    
    % Numeric integration (ODE45)
    if flag_mag == 0
        [T, Y] = ode45('rigbody', tspan, att_vec, options, ext_torq, iner, invin);
    else
        % To convert from inertial state vector to terrestrial vector
        geoc = inertial_to_terrestrial(gst(mjd, dfra+t), stat);

        % Earth's magnetic field
        sphe = rectangular_to_spherical(geoc);
        alt = sphe(3)/1000;
        elong = sphe(1);
        colat = pi/2 - sphe(2); 
        earth_field = 1.e-9*igrf_field(year_frac, alt, colat, elong)';

        % Magnetic Control
        k_p = 0.01;            % ganancia proporcional
        k_v = 0.04;           % ganancia derivativa
        eps = 0.05;            % epsilon
        earth_field_b = quatrmx(quat)*earth_field;

        dq = quat(1:3);
        dw = w_ang;
        
        if (norm(dw)<0.005 && quat(4)*quat(4)>0.1 && signq4==0) 
          signq4=sign(quat(4));
        end

        u = -eps*eps*k_p*dq*signq4 - eps*k_v*iner*dw;
        %gamavg = gamavg + Skew(earth_field_b)*Skew(-earth_field_b);
        %gamavg = gamavg / norm(gamavg);
        %u = pinv(gamavg)*u;
        
        normb2 = norm(earth_field)^2;

        mag_mom = 1/normb2 * cross(earth_field_b, u);
        
        if max(abs(mag_mom)) > maxmagmom             % torqrod saturation with bisection
          mag_mom = mag_mom*maxmagmom/max(abs(mag_mom));
        end
        
        ext_torq = ext_torq + cross(mag_mom, earth_field_b);

        [T, Y] = ode45('rigbody', tspan, att_vec, options, ext_torq, iner, invin, ...
            mag_mom, earth_field);
    end
    
    t

    att_vec = Y(3, :)';         % propagated attitude vector
    quat = att_vec(1:4);        % propagated quaternion
	w_ang = att_vec(5:7);       % propagated angular velocity

    eulzxz = quatezxz(quat);    % euler angles

    % attitude control torque logic (if any)
 	cont_torq = [0; 0; 0];
 
    % Store data to be plotted
    time = cat(2, time, t);
    euler = cat(2, euler, eulzxz*180/pi);
    omeg = cat(2, omeg, w_ang);
    orbit = cat(2, orbit, stat');
    keorb = cat(2, keorb, kep2');
    vdq = [vdq dq];
    vdw = [vdw dw]; 
    vmag_mom = [vmag_mom mag_mom];
    vext_torq = [vext_torq ext_torq];

end

close all

% Output visualization
figure(1)
plot(time/orb_period, vdq);
xlabel('Time (orbits)')
ylabel('Quaternion 1-3')
title('Attitude in quaternion vector')

figure(2)
plot(time/orb_period, omeg);
hold on;
xlabel('Time (orbits)')
ylabel('Angular velocity (rad/s)')
title('Attitude angular velocity')

figure(3)
subplot(3, 1, 1);
plot(time/orb_period, keorb(4,:)*180/pi);
xlabel('Time (orbits)')
ylabel('Ascencion node (deg)')
title('Right ascencion of the ascending node')

subplot(3, 1, 2);
plot(time/orb_period, keorb(5,:)*180/pi);
xlabel('Time (orbits)')
ylabel('Perigee argument (deg)')
title('Perigee argument')

subplot(3, 1, 3);
plot(time/orb_period, keorb(6,:)*180/pi);
xlabel('Time (orbits)')
ylabel('Mean anomaly (deg)')
title('Mean anomaly')

figure(4)
subplot(2, 1, 1);
plot(time/orb_period, orbit(1:3,:)/1000);
xlabel('Time (orbits)')
ylabel('Position (km)')
title('Satellite inertial position ')

subplot(2, 1, 2);
plot(time/orb_period, orbit(4:6,:)/1000);
xlabel('Time (orbits)')
ylabel('Velocity (km/s)')
title('Satellite velocity')

tini=ceil(length(time)*1.195/4);
figure(4);clf;
subplot(1,3,1);
title('Phase plane w_x/q_x');
plot(vdw(1,tini:end),vdq(1,tini:end),'r');hold on;grid on;xlabel("w_x [rad/sec]");ylabel("q_x");
subplot(1,3,2);
title('Phase plane w_y/q_y');
plot(vdw(2,tini:end),vdq(2,tini:end),'g');hold on;xlabel("w_y [rad/sec]");ylabel("q_y");grid on;
subplot(1,3,3);
title('Phase plane w_z/q_z');
plot(vdw(3,tini:end),vdq(3,tini:end),'b');hold on;xlabel("w_z [rad/sec]");ylabel("q_z");grid on;
grid on;

figure(22);clf;
plot3(1000*vext_torq(1,tini:end),1000*vext_torq(2,tini:end),1000*vext_torq(3,tini:end));grid on;
title('Magnetorquers Control Magnetic Torque [mNm]')
grid on

% Exponente de Lyapunov: requiere dos corridas, una se graba y compara
% luego contra la segunda. 
if dlyap==0
    
save('mwq1.mat','vdw','vdq','vmag_mom', 'vext_torq');

else

newomeg = vdw;
newvdq = vdq;
newext_torq = vext_torq;
newmag_mom = vmag_mom;

load('mwq1.mat','vdw','vdq','vmag_mom', 'vext_torq');

logdif = log( (1/dlyap)*sqrt( (newomeg(1,tini:end)-vdw(1,tini:end)).^2 + (newomeg(2,tini:end)-vdw(2,tini:end)).^2 + (newomeg(3,tini:end)-vdw(3,tini:end)).^2 + (newvdq(1,tini:end)-vdq(1,tini:end)).^2 + (newvdq(2,tini:end)-vdq(2,tini:end)).^2 + (newvdq(3,tini:end)-vdq(3,tini:end)).^2 ));

figure(23);clf;
title('Logarithm of attitude state differences norm');hold on;
plot(time(tini:end),logdif,'k');hold on;grid on;

end
