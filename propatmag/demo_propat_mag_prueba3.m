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

% Attitude elements in Euler angles of a 3-1-3 (z-x-z) rotation
eulzxz = [30, 50, 20]'*pi/180;   % converted from degrees to radians

% Attitude in quaternions
quat = ezxzquat(eulzxz);        % converted from Euler angles
quat = quat/norm(quat);

% Angular velocity vector in body frame:
w_ang = [10, 10, 10]'*pi/180;           % in radians/sec
    
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
tstep = 10;               % step time (sec)
tend = 20*orb_period;    % end time (10 minutes)

% Inertia matrix of axis-symetric rigid body:
iner = [0.005 0 0; 0 0.005 0; 0 0 0.005];         % in kg*m*m
%iner = [0.0022 0 0; 0 0.0020 0; 0 0 0.0022];         % in kg*m*m
%iner = [27 0 0; 0 17 0; 0 0 25];       % in kg*m*m
factorI = 5/2.2;

% Inverse inertia matrix:
invin = inv(iner);

% Initial control torque:
contq = [0 0 0]';

% Magnetic moment torque flag and moment:
flag_mag = 1;   % 1=compute magnetic moment / 0=discard magnetic moment
mag_mom = [0; 0; 0];      % in A.m

% ODE solver precision:
options = odeset('abstol', 1e-4, 'reltol', 1e-4);

% Initial vectors
time = tstart;          % to store time
euler = eulzxz*180/pi;  % Euler angles
omeg = w_ang;           % Angular velocity
orbit = stat';          % Orbit elements (state vector)
keorb = kepel';         % Orbit elements (keplerian)

gamavg = 0*eye(3);      % Gamma avg inicial
maxmagmom = 0.02*5;       % Tope momento magnético en Am2
vdq = [0;0;0];
vdqs = [0;0;0];
vdw = [0;0;0];
u = [0;0;0];
vu = u;
signq4 = 0;
vext_torq = u;
vmag_mom = u;

earthradius = 6371000;

% Attitude and orbit propagation
for t = tstart:tstep:tend

    % Orbit propagation
    kep2 = kepel + delk*t;

    % To convert from keplerian elements to state vector (if needed)
    stat = kepel_statvec(kep2);

    % Perturbation torques:
    A    = quatrmx(quat);
    exb  = A(:,1);
    ambt = 3*omeg_0^2*cross(exb,iner*exb);

    % Peor caso
    %ambt = 1e-10*[1 1 1]';
    
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

        earth_field_ned = 1.e-9*igrf_field(year_frac, alt, colat, elong)';
        earth_field_ecef = ned2ecef(earth_field_ned',[sphe(2);sphe(1)]);
        earth_field_eci = terrestrial_to_inertial(gst(mjd, dfra+t), [earth_field_ecef 0 0 0])';
        earth_field_eci = earth_field_eci(1:3);
        earth_field_b = quatrmx(quat)*earth_field_eci;

        % Magnetic Control
        k_p = 0.00001*factorI;           % ganancia proporcional
        k_v = 0.75/factorI;                  % ganancia derivativa
        eps = 0.01;                         % epsilon
        %earth_field_b = quatrmx(quat)*earth_field;

        dq = quat(1:3);
        dw = w_ang;
        
        dqs = cross(quatrmx(quat)*sun_dir(mjd, dfra+t),[0;0;1]);
        angles = asin(norm(dqs));
        versors = dqs/norm(dqs);
        dqs = sin(angles/2)*versors;
        dqs4 = cos(angles/2);
                
        %if (norm(dw)<0.001 && quat(4)*quat(4)>0.1 && signq4==0) 
        %  signq4=sign(quat(4));
        %end
        if (norm(dw)<0.002 && dqs4*dqs4>0.1 && signq4==0) 
          signq4=sign(dqs4);
        end
        
        % Control Derivativo
        u = - eps*k_v*iner*dw;
        
        % Eclipse
        eclipse = 0;
        if (orbit(1:3,end)'*sun_dir(mjd, dfra+t)<0)
            perpsun = (eye(3) - sun_dir(mjd, dfra+t)*sun_dir(mjd, dfra+t)')*orbit(1:3,end);
            if (norm(perpsun)<earthradius)
                eclipse = 1;
            end
        end
        %eclipse=0;
        
        % Sun pointing: término proporcional
        if (eclipse == 0 && signq4*signq4>0)
            u = u - eps*eps*k_p*inv(iner)*dqs*signq4;
        end
        
        if (eclipse == 1 && signq4*signq4>0)
            u = 20*u;
        end
        
        earth_field_bx = earth_field_b;
        
        normb2 = norm(earth_field)^2;

        mag_mom = 1/normb2 * cross(earth_field_b, u);
        
        if max(abs(mag_mom)) > maxmagmom             % torqrod saturation with bisection
          mag_mom = mag_mom*maxmagmom/max(abs(mag_mom));
        end
        
        mag_mom = mag_mom + maxmagmom*[0.003;0.003;-0.003];
        
        ext_torq = ext_torq + cross(mag_mom, earth_field_b);

        [T, Y] = ode45('rigbody', tspan, att_vec, options, ext_torq, iner, invin, ...
            mag_mom, earth_field);
    end
    
    t

    att_vec = Y(3, :)';         % propagated attitude vector
    quat = att_vec(1:4);        % propagated quaternion
	w_ang = att_vec(5:7);       % propagated angular velocity

    %eulzxz = quatezxz(quat);    % euler angles

    % attitude control torque logic (if any)
 	cont_torq = [0; 0; 0];
 
    % Store data to be plotted
    time = cat(2, time, t);
    %euler = cat(2, euler, eulzxz*180/pi);
    omeg = cat(2, omeg, w_ang);
    orbit = cat(2, orbit, stat');
    keorb = cat(2, keorb, kep2');
    vdq = [vdq dq];
    vdqs = [vdqs dqs];
    vdw = [vdw dw]; 
    vmag_mom = [vmag_mom mag_mom];
    vext_torq = [vext_torq ext_torq];

end

close all

% Output visualization
figure(1)
plot(time/orb_period, vdq(1,:),'r');hold on;
plot(time/orb_period, vdq(2,:),'g');hold on;
plot(time/orb_period, vdq(3,:),'b');hold on;
xlabel('Time (orbits)')
ylabel('Quaternion 1-3')
title('Attitude in quaternion vector');
grid on;

figure(10)
plot(time/orb_period, vdqs(1,:),'r');hold on;
plot(time/orb_period, vdqs(2,:),'g');hold on;
plot(time/orb_period, vdqs(3,:),'b');hold on;
xlabel('Time (orbits)')
ylabel('Quaternion 1-3')
title('Attitude in partial (Sun) quaternion vector');
grid on;

figure(2)
plot(time/orb_period, omeg(1,:),'r');hold on;
plot(time/orb_period, omeg(2,:),'g');hold on;
plot(time/orb_period, omeg(3,:),'b');hold on;
plot(time/orb_period,sqrt(omeg(1,:).^2+omeg(2,:).^2+omeg(3,:).^2),'k');hold on;
hold on;
xlabel('Time (orbits)')
ylabel('Angular velocity (rad/s)')
title('Attitude angular velocity')
grid on;

figure(22);clf;
plot(time/orb_period,vmag_mom(1,:),'r');hold on;
plot(time/orb_period,vmag_mom(2,:),'g');hold on;
plot(time/orb_period,vmag_mom(3,:),'b');hold on;
plot(time/orb_period,maxmagmom*ones(size(vext_torq(3,:))),'k--');hold on;
plot(time/orb_period,-maxmagmom*ones(size(vext_torq(3,:))),'k--');hold on;
title('Magnetorquers Control Magnetic Moment [Nm]')
grid on;