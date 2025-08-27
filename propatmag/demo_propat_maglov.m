% demo_propat
%
%       Program to show how to use PROPAT
%       Type demo_propat from the Matlab prompt

% Inputs:
% Orbit keplerian elements:
kepel = [7000000, 0.01, 98, 0, 0, 0];   % see function delkep
kepel = [6.378135e06 + 5.26043425e+05, 0.0011, 1.70171324, 0, pi/2, 0];
kepel = [6.378135e06 + 450000, 0.0, 87*pi/180, 0, 0, 0];   % see function delkep

%H1=200000;
%kepel = [(70550000+6371000*2+H1)/2, (70550000-H1)/(70550000+6371000*2+H1), 39*pi/180, 0, pi/4, 1.6*pi];   % see function delkep


% SSO frozen repeat every 2 days
D       = 1;
R       = 15*D-1;
mu      = 59800000000.0*6669.0;
J2      = 1.0826268362e-3;
J3      = -2.538810e-6;
ae      = 6378139.;
To      = 86400*D/R;
ao      = (mu*(To/(2*pi))^2)^(1/3);
isso    = acos(-(sqrt(ao/12352000))^7);
s       = 365.242199;
we      = 360/86400*(1+1/s)*pi/180;
dL      = 360/86400*To*pi/180;
hoe     = ao-ae;
issodeg = isso*180/pi;
dlm     = 2*pi/R*ae;      % separación de dos cruces sobre el ecuador
efr     = -J3/J2*sin(isso)*ae/(2*ao);

kepel = [ao, efr, isso, 0, pi/2, 0];

% Orbit state vector:
stat = kepel_statvec(kepel);

% Velocidad Orbital
GM = 3.986004418e14;
orb_period = sqrt(4*pi^2*kepel(1)^3/GM);

omeg_0 = 2*pi/orb_period;

dlyap = 0.0001*0;

% Attitude elements in Euler angles of a 3-1-3 (z-x-z) rotation
eulzxz = [10, 0, 0]'*pi/180;   % converted from degrees to radians

% Attitude in quaternions
quat = ezxzquat(eulzxz)  + sqrt(2)*dlyap*[1;-1;1;0]/sqrt(3);        % converted from Euler angles
quat = quat/norm(quat);
qv = [-0.4;0.55;-0.2];
quat = [qv;sqrt(1-qv'*qv)];
quat = quat/norm(quat);

% Angular velocity vector in body frame:
w_ang = [0, 0, 0]' + sqrt(2)*dlyap*[1;-1;1]/sqrt(3);           % in radians/sec

w_ang = [0.000000001;0;0]+sqrt(2)*dlyap*[1;-1;1]/sqrt(3);
quat = [0;0;0;1];

floquet = 1:6; % 0=NO armar matriz fundamental / 1:6=calc. col de Matriz Fundamental.

%for floquet=1:6

floquet = 0;

if 1

if floquet
    w_ang = [0 0 0]';
    quat  = [0 0 0 1]';
    if floquet >= 1 && floquet < 4
        quat(floquet)  = 0.75;
        quat(4)         = sqrt(1-quat(floquet)^2);
    elseif floquet < 7
        w_ang(floquet - 3) = 0.01;
    end
    quat   = quat/norm(quat);
    x_floq = [quat(1:3); w_ang];
    tend   = (86400*D);
end


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
tend = 14*orb_period;    % end time (10 minutes)

% Inertia matrix of axis-symetric rigid body:
iner = [9 0.2 0.2; 0.2 10 0;0.2 0 11];      % in kg*m*m
%iner = [27 0 0; 0 17 0; 0 0 25];           % in kg*m*m ROSSA

% Inverse inertia matrix:
invin = inv(iner);

% Initial control torque:
contq = [0 0 0]';

% Magnetic moment torque flag and moment:
flag_mag = 1;               % 1=compute magnetic moment / 0=discard magnetic moment
mag_mom = [0; 0; 0.1];      % in A.m

% ODE solver precision:
options = odeset('abstol', 1e-5, 'reltol', 1e-5);

% Initial vectors
time = tstart;          % to store time
euler = eulzxz*180/pi;  % Euler angles
omeg = w_ang;           % Angular velocity
orbit = stat';          % Orbit elements (state vector)
keorb = kepel';         % Orbit elements (keplerian)

gamavg = 0*eye(3);  % Gamma avg inicial
maxmagmom = 2*10000000;      % Tope momento magnético en Am2
vdq = [0;0;0];
vdw = [0;0;0];
u = [0;0;0];
vu = u;
signq4 = 0;
vext_torq = u;
vmag_mom = u;

delta=[0;0;0;0];
alfa = 0.1;
lambda = 0.1;
vdelta=delta;
veKA = [vdq];
eKA = [vdq];
vambt = [0;0;0];
vu = [0;0;0];

ggenable = 0;

mnb2 = 0;

% Attitude and orbit propagation
for t = tstart:tstep:tend

    % Orbit propagation
    kep2 = kepel + delk*t;
    
    % To convert from keplerian elements to state vector (if needed)
    stat = kepel_statvec(kep2);

    % Perturbation torques:
    ambt = [0 0 0]';
    A    = quatrmx(quat);
    ambt = ggenable * 3*omeg_0^2*cross(A*stat(1:3)',iner*A*stat(1:3)')/norm(stat(1:3))^2;
    
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
        earth_field_ned = 1.e-9*igrf_field(year_frac, alt, colat, elong)';
        earth_field_ecef = ned2ecef(earth_field_ned',[sphe(2);sphe(1)]);
        earth_field_eci = terrestrial_to_inertial(gst(mjd, dfra+t), [earth_field_ecef 0 0 0])';
        earth_field_eci = earth_field_eci(1:3);
        earth_field_b = quatrmx(quat)*earth_field_eci;
        
        
        Aq = (quat(4)^2-quat(1:3)'*quat(1:3))*eye(3)+2*quat(1:3)*quat(1:3)'+2*Skew(quat(1:3));
        om = 0.001017846057324;
        b0 = [[-1.8806 -6.2168 -36.0046];[-1.0787 -0.4385 -1.0345];[-11.7084 35.6763 -6.4530]]*[1;cos(om*t);sin(om*t)];
        b0 = b0*4.315553469313078e-05/48.183430390760840;
        b  = Aq*b0;
        earth_field_b = b;
        earth_field_eci=b0;

        dq = quat(1:3);
        dw = w_ang;

        if (norm(dw)<0.005 && signq4==0)
          signq4=sign(quat(4));
        end

        % Magnetic Control
        k_p = 0.02;
        k_p = 0.002*max(max(iner));            % ganancia proporcional
        k_v = 0.05;            % ganancia derivativa
        eps = 0.1*10;            % epsilon

        u = -eps*eps*k_p*inv(iner)*dq*signq4 - eps*k_v*iner*dw - ambt + 0*cross(dw,iner*dw);
        
        %u = -eps*eps*k_p*dq*signq4 - eps*k_v*dw;
        %u = inv(iner)*u;
        
        normb2 = norm(earth_field_b)^2;

        mag_mom = 1/normb2 * cross(earth_field_b, u);

        if max(abs(mag_mom)) > maxmagmom             % torqrod saturation with bisection
          mag_mom = mag_mom*maxmagmom/max(abs(mag_mom));
        end

        ext_torq = ext_torq + cross(mag_mom, earth_field_b);

        [T, Y] = ode45('rigbody', tspan, att_vec, options, ext_torq, iner, invin, ...
            mag_mom, earth_field_eci);
    
        if normb2>mnb2
            mnb2=normb2;
        end;
        
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
    vdelta=[vdelta delta];
    veKA = [veKA eKA];
    vambt = [vambt ambt];
    vu = [vu u];

end

if floquet
    x_floq = [x_floq att_vec([1:3 5:7])];
    save(['GGfloquet_' num2str(floquet) '.mat'], 'x_floq')
    %return
end

end % floquet

if floquet

[X0, XT] = floquet_fundamental;
B = inv(X0)*XT;
mul_floquet = eig(B);
exp_floquet = log(mul_floquet)/(14*orb_period)

save(['GGepsilon_' num2str(eps) '.mat'], 'exp_floquet')

end;

%return;

if floquet

vexpo = [];
load(['epsilon_0_025.mat']);    vexpo=[vexpo exp_floquet];
load(['epsilon_0_05.mat']);     vexpo=[vexpo exp_floquet];
load(['epsilon_0_075.mat']);    vexpo=[vexpo exp_floquet];
load(['epsilon_0_1.mat']);      vexpo=[vexpo exp_floquet];
load(['epsilon_0_125.mat']);    vexpo=[vexpo exp_floquet];
load(['epsilon_0_185.mat']);    vexpo=[vexpo exp_floquet];
load(['epsilon_0_225.mat']);    vexpo=[vexpo exp_floquet];
load(['epsilon_0_235.mat']);    vexpo=[vexpo exp_floquet];
load(['epsilon_0_25.mat']);     vexpo=[vexpo exp_floquet];

veps = [0.025 0.05 0.075 0.1 0.125 0.185 0.225 0.235 0.25];

figure(99);clf;
plot(veps,real(vexpo));hold on;
title('Parte real de los exponentes característicos');hold on;
xlabel('\epsilon');hold on;
ylabel('Real(\mu_i)');hold on;
grid on;
axis([veps(1) veps(9) -6e-4 1e-4]);


vexpo = [];
load(['GGepsilon_0_025.mat']);    vexpo=[vexpo exp_floquet];
load(['GGepsilon_0_0125.mat']);     vexpo=[vexpo exp_floquet];
load(['GGepsilon_0_005.mat']);    vexpo=[vexpo exp_floquet];

veps = [0.005 0.0125 0.025];

figure(999);clf;
plot(veps,real(vexpo));hold on;
title('Parte real de los exponentes característicos');hold on;
xlabel('\epsilon');hold on;
ylabel('Real(\mu_i)');hold on;
grid on;
axis([veps(1) veps(3) -6e-4 1e-4]);

end;

close all


figure(222222);
plot(time(1,tini:end)/orb_period, sqrt(vdq(1,tini:end).^2+vdq(2,tini:end).^2+vdq(3,tini:end).^2+vdw(1,tini:end).^2+vdw(2,tini:end).^2+vdw(3,tini:end).^2), 'm');hold on;
title('|x|');hold on;xlabel('Órbitas');grid on;hold on;


% Output visualization
figure(1)
plot(time/orb_period, vdq(1,:),'r');hold on;
plot(time/orb_period, vdq(2,:),'g');hold on;
plot(time/orb_period, vdq(3,:),'b');hold on;
%plot(time/orb_period, sqrt(vdq(1,:).^2+vdq(2,:).^2+vdq(3,:).^2),'k');hold on;
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

%figure(4)
%subplot(2, 1, 1);
%plot(time/orb_period, orbit(1:3,:)/1000);
%xlabel('Time (orbits)')
%ylabel('Position (km)')
%title('Satellite inertial position ')

%subplot(2, 1, 2);
%plot(time/orb_period, orbit(4:6,:)/1000);
%xlabel('Time (orbits)')
%ylabel('Velocity (km/s)')
%title('Satellite velocity')

tini=ceil(length(time)*.5/4);
tini=1;
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

figure(222);clf;
subplot(2,1,1);
plot(time(1,tini:end)/orb_period, 1000*vext_torq(1,tini:end), 'r');hold on;
plot(time(1,tini:end)/orb_period, 1000*vext_torq(2,tini:end), 'g');hold on;
plot(time(1,tini:end)/orb_period, 1000*vext_torq(3,tini:end), 'b');hold on;
title('Magnetorquers Control Magnetic Torque [mNm]')
subplot(2,1,2);
plot(time(1,tini:end)/orb_period, vmag_mom(1,tini:end), 'r');hold on;
plot(time(1,tini:end)/orb_period, vmag_mom(2,tini:end), 'g');hold on;
plot(time(1,tini:end)/orb_period, vmag_mom(3,tini:end), 'b');hold on;
title('Magnetorquers Control Magnetic Moment [Am2]')
grid on


figure(2222);
subplot(2,1,1);
plot(time(1,tini:end)/orb_period, 1000*sqrt(vext_torq(1,tini:end).^2+vext_torq(2,tini:end).^2+vext_torq(3,tini:end).^2), 'm');hold on;
title('Torque de Control [mNm]');hold on;xlabel('Órbitas');grid on;hold on;
axis([0 2 0 0.027]);
subplot(2,1,2);
plot(time(1,tini:end)/orb_period, asin(sqrt(vdq(1,tini:end).^2+vdq(2,tini:end).^2+vdq(3,tini:end).^2))*2*180/pi, 'm');hold on;
title('Ángulo de error [grados]');hold on;xlabel('Órbitas');grid on;
axis([0 2 0 90]);


figure(22222);
subplot(2,1,1);
plot(time(1,tini:end)/orb_period, log(sqrt(vdw(1,tini:end).^2+vdw(2,tini:end).^2+vdw(3,tini:end).^2)), 'm');hold on;
title('Logaritmo de la norma de la velocidad angular');hold on;xlabel('Órbitas');grid on;hold on;
%axis([0 2 0 0.027]);
subplot(2,1,2);
plot(time(1,tini:end)/orb_period, log(asin(sqrt(vdq(1,tini:end).^2+vdq(2,tini:end).^2+vdq(3,tini:end).^2))*2*180/pi), 'm');hold on;
title('Logaritmo del ángulo de error [grados]');hold on;xlabel('Órbitas');grid on;
%axis([0 2 0 90]);


figure(222223);clf;
subplot(2,1,1);
plot(time(1,tini:end)/orb_period, (vdw(1,tini:end).^2+vdw(2,tini:end).^2+vdw(3,tini:end).^2).^0.5, 'm');hold on;
title('Norma de la velocidad angular');hold on;xlabel('Órbitas');grid on;hold on;
subplot(2,1,2);
plot(time(1,tini:end)/orb_period, asin(sqrt(vdq(1,tini:end).^2+vdq(2,tini:end).^2+vdq(3,tini:end).^2))*2*180/pi, 'm');hold on;
title('Ángulo de error [grados]');hold on;xlabel('Órbitas');grid on;



figure(222224);clf;
plot(time(1,tini:end)/orb_period, log((vdq(1,tini:end).^2+vdq(2,tini:end).^2+vdq(3,tini:end).^2+vdw(1,tini:end).^2+vdw(2,tini:end).^2+vdw(3,tini:end).^2).^0.5), 'm');hold on;
title('Norma del estado');hold on;xlabel('Órbitas');grid on;hold on;


%figure(25);clf;
%plot(time,vdelta(1,:),'r');hold on;
%plot(time,vdelta(2,:),'g');hold on;
%plot(time,vdelta(3,:),'b');hold on;
%plot(time,vdelta(4,:),'k');hold on;
%title('delta')
%grid on

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

tini=1;
logdif = log( (1/dlyap)*sqrt( (newomeg(1,tini:end)-vdw(1,tini:end)).^2 + (newomeg(2,tini:end)-vdw(2,tini:end)).^2 + (newomeg(3,tini:end)-vdw(3,tini:end)).^2 + (newvdq(1,tini:end)-vdq(1,tini:end)).^2 + (newvdq(2,tini:end)-vdq(2,tini:end)).^2 + (newvdq(3,tini:end)-vdq(3,tini:end)).^2 ));

figure(23);clf;
%title('Logarithm of attitude state differences norm');hold on;
title('Logaritmo de la norma de las diferencias entre soluciones inicialmente próximas');hold on;
xlabel('Órbitas');
plot(time/orb_period,logdif,'k');hold on;grid on;

end