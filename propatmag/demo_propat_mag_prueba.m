% demo_propat
%
%       Program to show how to use PROPAT
%       Type demo_propat from the Matlab prompt

% Inputs:
% Orbit keplerian elements:
%kepel = [6786000, 0.0005, 1.7, 0, 0, 0];   % see function delkep

clear all
clc

timestamp = datestr(now,'yyyymmdd_HHMMSS');

%---------------------------CONFIGURACIÓN----------------------------------

ECLIPSE = 0;
RMM = 0;
GRAV_GRAD = 0;
SOLAR_TORQ = 0;
DRAG = 0;

SUN_POINTING = 0;

RMM_ESTIMATE = 0;
RMM_COMPENSATE = 0;
Q_ESTIMATE = 0;

SAVE_FIG = 0;


%-------------------------------ORBITA-------------------------------------

kepel = [7000000, 0.01, 95*pi/180, 0, 0, 0];

% Orbit state vector:
stat = kepel_statvec(kepel);

% Compute the variations in keplerian elements due to the Earth oblateness
delk = delkep(kepel);

% Velocidad Orbital
GM = 3.986004418e14;
orb_period = sqrt(4*pi^2*kepel(1)^3/GM);
omeg_0 = 2*pi/orb_period;

%--------------------------VALORES INICIALES-------------------------------

% Attitude elements in Euler angles of a 3-1-3 (z-x-z) rotation
eulzxz = [30, 50, 20]'*pi/180;   % converted from degrees to radians

% Attitude in quaternions
quat = ezxzquat(eulzxz);        % converted from Euler angles
quat = quat/norm(quat);

% Angular velocity vector in body frame:
w_ang = [10, 10, 10]'*pi/180;           % in radians/sec

% Initial control torque:
contq = [0 0 0]';
    

%-----------------------------EFEMERIDES-----------------------------------

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
tstep = 2;               % step time (sec)
tend = 6*orb_period;    % end time (10 minutes)


%-----------------------------DINAMICA-------------------------------------

% Inertia matrix of axis-symetric rigid body:
iner = [0.0022 0 0; 0 0.0022 0; 0 0 0.0022];         % in kg*m*m
%iner = [27 0 0; 0 17 0; 0 0 25];       % in kg*m*m

% Inverse inertia matrix:
invin = inv(iner);

% Magnetic moment torque flag and moment:
flag_mag = 1;   % 1=compute magnetic moment / 0=discard magnetic moment
mag_mom = [0; 0; 0];      % in A.m

% ODE solver precision:
options = odeset('abstol', 1e-4, 'reltol', 1e-4);

earthradius = 6371000;

%------------------------------CONTROLADOR---------------------------------

k_p = 0.000025;
k_v = 0.3; 
eps = 0.01;            

controller = controller(k_p, k_v, eps, iner);

%--------------------------------KALMAN------------------------------------

H = [eye(3) zeros(3)];
Q = [0.0001^2*eye(3) zeros(3)
    zeros(3) 0.001^2*eye(3)];
Q_min = Q;

alpha = 0;
if Q_ESTIMATE
    alpha = 0.001;
end

R = 0.01*eye(3);

x_pred = [0; 0; 0; 0; 0; 0];
sigma = [[eye(3)*0.00001^2 zeros(3)];[zeros(3) eye(3)*0.001^2]]*2000;

ekf_rmm = rmm_estimator(x_pred, sigma, Q, R, H, Q_min, alpha, iner, controller, tstep);

%------------------------------OBSERVABILIDAD------------------------------

obs_1 = observability_calculator(H, 2);
obs_2 = observability_calculator(H, 3);
obs_5 = observability_calculator(H, 5);
obs_10 = observability_calculator(H, 10);
obs_50 = observability_calculator(H, 50);

%------------------------------PRESION SOLAR------------------------------

T_sun = 5.772e3; %Kelvin
S_sat = 0.01; %m^2 superficie la cual le da el sol
C_r = 1; %Coeficiente de relexion del satelite
r_cp = [0.001 0 0];%Distancia del centro de masas al centro de presion

%-----------------------VECTORES PARA GRAFICAR-----------------------------

% Initial vectors
time = tstart;          % to store time
euler = eulzxz*180/pi;  % Euler angles
omeg = w_ang;           % Angular velocity
orbit = stat';          % Orbit elements (state vector)
keorb = kepel';         % Orbit elements (keplerian)

gamavg = 0*eye(3);      % Gamma avg inicial
maxmagmom = 0.02;       % Tope momento magnético en Am2
vdq = [0;0;0];
vdqs = [0;0;0];
vdw = [0;0;0];
u = [0;0;0];
vu = u;
signq4 = 0;
vext_torq = u;
vmag_mom = u;
vrmm_hat = [0;0;0];
vdw_hat = [0;0;0];
vrmm_diag_cov = [sigma(4,4); sigma(5,5); sigma(6,6)];
vQ = [Q(4,4); Q(5,5); Q(6,6)];
vsingularM_1 = 0;
vsingularM_2 = 0;
vsingularM_5 = 0;
vsingularM_10 = 0;
vsingularM_50 = 0;

%------------------------------SIMULACION----------------------------------

for t = tstart:tstep:tend
%for t = 1:tstep

    % Orbit propagation
    kep2 = kepel + delk*t;
    % To convert from keplerian elements to state vector (if needed)
    stat = kepel_statvec(kep2);

    %-------------------------PERTURBACIONES-------------------------------
    grad_grav = 0;
    if GRAV_GRAD == 1
        A    = quatrmx(quat);
        exb  = A(:,1);
        grad_grav = 3*omeg_0^2*cross(exb,iner*exb);
    end
    
    mom_res = 0;
    if RMM == 1
        mom_res = maxmagmom*[0.005; 0.005; -0.005];
    end
    
    % Peor caso
    %ambt = 2e-10*[1 1 1]';
    ambt = 0;
    
    % External torques (perturbation)
    ext_torq = ambt + grad_grav;
    
    % Initial attitude vector:
    att_vec = [quat; w_ang]';         % Transposed

    % ODE Solver parameters
    tspan = [t, t+tstep/2, t+tstep];
    
    % Numeric integration (ODE45)
    if flag_mag == 0
        [T, Y] = ode45('rigbody', tspan, att_vec, options, ext_torq, iner, invin);
    else
        
        %-----------------------CAMPO MAGNETICO----------------------------
        
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
        
        %--------------------------ECLIPSE---------------------------------

        
        eclipse = 0;
        if (orbit(1:3,end)'*sun_dir(mjd, dfra+t)<0 && ECLIPSE == 1)
            perpsun = (eye(3) - sun_dir(mjd, dfra+t)*sun_dir(mjd, dfra+t)')*orbit(1:3,end);
            if (norm(perpsun)<earthradius)
                eclipse = 1;
            end
        end

        %---------------------CONTROL MAGNETICO----------------------------
        
        
        k_p = 0.00025;% ganancia proporcional
        k_v = 0.3; % ganancia derivativa
        eps = 0.01;% epsilon

        dq = quat(1:3);
        dw = w_ang;
        
        dqs = cross(quatrmx(quat)*sun_dir(mjd, dfra+t),[0;0;1]);
        %dqs = cross(quatrmx(quat)*stat(1:3)'/norm(stat(1:3)),[0;0;1])
        angles = asin(norm(dqs));
        versors = dqs/norm(dqs);
        dqs = sin(angles/2)*versors;
        dqs4 = cos(angles/2); 
        
        if (norm(dw)<0.001 && dqs4*dqs4>0.1 && signq4==0) 
           signq4=sign(dqs4);
        end
        
        % Control Derivativo
        %u = - eps*k_v*iner*dw;
        pointing = 0;
        % Sun pointing: término proporcional
        if (eclipse == 0 && signq4*signq4>0 && SUN_POINTING)
            pointing = 1;
            %u = u - eps*eps*k_p*inv(iner)*dqs*signq4;
            %u = u - eps*eps*k_p*dqs*signq4;
        end
        
        u = get_control_action(controller, dqs, signq4, dw, pointing);
        
        
        %----------------ESTIMACION DE MOMENTO RESIDUAL--------------------
        phik = 0;
        rmm_avas = [0;0;0];
        if RMM_ESTIMATE == 1
            %[x_pred, sigma, phikm1, Q] = ekf_rmm(x_pred, (mag_mom - mom_res), iner, earth_field_b, sigma, dw, Q, Q_min, alpha, R, tstep);
            [ekf_rmm, phik] = update(ekf_rmm, mag_mom - mom_res, earth_field_b, dw);
            [x_pred, sigma, Q] = get_estimates(ekf_rmm);
            avas_sigma = eig(sigma);
            rmm_avas = avas_sigma(4:6);
        end
        %---------------------MOMENTO MAGNÉTICO----------------------------
        
        normb2 = norm(earth_field)^2;

        mag_mom = 1/normb2 * cross(earth_field_b, u) + mom_res - RMM_COMPENSATE*x_pred(4:6);
        
        
        if max(abs(mag_mom)) > maxmagmom             % torqrod saturation with bisection
          mag_mom = mag_mom*maxmagmom/max(abs(mag_mom));
        end
        
        ext_torq = ext_torq + cross(mag_mom, earth_field_b);
        
        %---------------------MOMENTO SOLAR----------------------------
        
        sun_torq = 0;
        if (SOLAR_TORQ == 1 && eclipse == 0)
            s_eci = sun_dir(mjd, dfra + t); % sun vector en ECI
            v_sun_body = quatrmx(quat) * s_eci(:); %sun vector en body
            norm_sun = norm(v_sun_body);
            r_presion = (v_sun_body/norm_sun)*0.001;
            [F_sun, sun_torq] = solar_pressure(v_sun_body, T_sun, S_sat, r_presion', C_r);
        end
        ext_torq = ext_torq + sun_torq';
        
        %-------------------------OBSERVABILIDAD---------------------------

        obs_1 = update(obs_1, phik);
        obs_2 = update(obs_2, phik);
        obs_5 = update(obs_5, phik);
        obs_10 = update(obs_10, phik);
        obs_50 = update(obs_50, phik);
        
        O_1 = get_M(obs_1);
        O_2 = get_M(obs_2);
        O_5 = get_M(obs_5);
        O_10 = get_M(obs_10);
        O_50 = get_M(obs_50);
        sing_O_1 = log10(svds(O_1, 1, 'smallest'));
        sing_O_2 = log10(svds(O_2, 1, 'smallest'));
        sing_O_5 = log10(svds(O_5, 1, 'smallest'));
        sing_O_10 = log10(svds(O_10, 1, 'smallest'));
        sing_O_50 = log10(svds(O_50, 1, 'smallest'));
        
        %------------------------------------------------------------------
        
        [T, Y] = ode45('rigbody', tspan, att_vec, options, ext_torq, iner, invin, ...
            mag_mom, earth_field);
    end
    
    t

    att_vec = Y(3, :)';         % propagated attitude vector
    quat = att_vec(1:4);        % propagated quaternion
	w_ang = att_vec(5:7);       % propagated angular velocity

    %eulzxz = quatezxz(quat);    % euler angles
 
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
    vrmm_hat = [vrmm_hat x_pred(4:6)];
    vdw_hat = [vdw_hat x_pred(1:3)];
    vrmm_diag_cov = [vrmm_diag_cov rmm_avas];
    vQ = [vQ [Q(4,4); Q(5,5); Q(6,6)]];
    vsingularM_1 = [vsingularM_1; sing_O_1];
    vsingularM_2 = [vsingularM_2; sing_O_2];
    vsingularM_5 = [vsingularM_5; sing_O_5];
    vsingularM_10 = [vsingularM_10; sing_O_10];
    vsingularM_50 = [vsingularM_50; sing_O_50];

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
if SAVE_FIG == 1
    print(gcf, ['Quat' timestamp '.png'], '-dpng')
end

figure(10)
plot(time/orb_period, vdqs(1,:),'r');hold on;
plot(time/orb_period, vdqs(2,:),'g');hold on;
plot(time/orb_period, vdqs(3,:),'b');hold on;
xlabel('Time (orbits)')
ylabel('Quaternion 1-3')
title('Attitude in partial (Sun) quaternion vector');
grid on;
if SAVE_FIG == 1
    print(gcf, ['Quat_sun' timestamp '.png'], '-dpng')
end

figure(2)
plot(time/orb_period, omeg(1,:),'r');hold on;
plot(time/orb_period, omeg(2,:),'g');hold on;
plot(time/orb_period, omeg(3,:),'b');hold on;
plot(time/orb_period,sqrt(omeg(1,:).^2+omeg(2,:).^2+omeg(3,:).^2),'k');hold on;

plot(time/orb_period,vdw_hat(1,:),'r--');hold on;
plot(time/orb_period,vdw_hat(2,:),'g--');hold on;
plot(time/orb_period,vdw_hat(3,:),'b--');hold on;

xlabel('Time (orbits)')
ylabel('Angular velocity (rad/s)')
title('Attitude angular velocity')
grid on;
if SAVE_FIG == 1
    print(gcf, ['Velocidad_angular' timestamp '.png'], '-dpng')
end
hold on;


figure(22);clf;
plot(time/orb_period,vmag_mom(1,:),'r');hold on;
plot(time/orb_period,vmag_mom(2,:),'g');hold on;
plot(time/orb_period,vmag_mom(3,:),'b');hold on;
plot(time/orb_period,maxmagmom*ones(size(vext_torq(3,:))),'k--');hold on;
plot(time/orb_period,-maxmagmom*ones(size(vext_torq(3,:))),'k--');hold on;
title('Magnetorquers Control Magnetic Moment [Am2]')
grid on;
if SAVE_FIG == 1
    print(gcf, ['Momento_mag' timestamp '.png'], '-dpng')
end

figure(23);clf;
plot(time/orb_period,vrmm_hat(1,:),'r');hold on;
plot(time/orb_period,vrmm_hat(2,:),'g');hold on;
plot(time/orb_period,vrmm_hat(3,:),'b');hold on;
%plot(time/orb_period,maxmagmom*ones(size(vext_torq(3,:))),'k--');hold on;
%plot(time/orb_period,-maxmagmom*ones(size(vext_torq(3,:))),'k--');hold on;
title('Estimated RMM [Am2]')
grid on;
if SAVE_FIG == 1
    print(gcf, ['RMM_hat' timestamp '.png'], '-dpng')
end


figure(24);clf;
plot(time/orb_period,vrmm_diag_cov(1,:),'r');hold on;
plot(time/orb_period,vrmm_diag_cov(2,:),'g');hold on;
plot(time/orb_period,vrmm_diag_cov(3,:),'b');hold on;
%plot(time/orb_period,maxmagmom*ones(size(vext_torq(3,:))),'k--');hold on;
%plot(time/orb_period,-maxmagmom*ones(size(vext_torq(3,:))),'k--');hold on;
title('Covariance od estimated RMM [Am2]')
grid on;
if SAVE_FIG == 1
    print(gcf, ['Cov_RMM_hat' timestamp '.png'], '-dpng')
end


figure(25);clf;
plot(time/orb_period,vQ(1,:),'r');hold on;
plot(time/orb_period,vQ(2,:),'g');hold on;
plot(time/orb_period,vQ(3,:),'b');hold on;
%plot(time/orb_period,maxmagmom*ones(size(vext_torq(3,:))),'k--');hold on;
%plot(time/orb_period,-maxmagmom*ones(size(vext_torq(3,:))),'k--');hold on;
title('Estimated Q [Am2]')
grid on;
if SAVE_FIG == 1
    print(gcf, ['Cov_RMM_hat' timestamp '.png'], '-dpng')
end

figure(26);clf;
plot(time/orb_period,vsingularM_1);hold on;
plot(time/orb_period,vsingularM_2);hold on;
plot(time/orb_period,vsingularM_5);hold on;
plot(time/orb_period,vsingularM_10);hold on;
plot(time/orb_period,vsingularM_50);hold on;
legend('$O_1$', '$O_2$', '$O_5$', '$O_{10}$', '$O_{50}$' )
title('Minimum singular values of $O_{\tau}$', 'Interpreter', 'latex')
grid on;
if SAVE_FIG == 1
    print(gcf, ['Cov_RMM_hat' timestamp '.png'], '-dpng')
end