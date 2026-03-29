% demo_propat
%
%       Program to show how to use PROPAT
%       Type demo_propat from the Matlab prompt

% Inputs:
% Orbit keplerian elements:
RE = 6378000;
kepel = [RE+450000*0+522000+621000*0+721000*0, 0.0013, 97.5*pi/180, 1*pi/180, 0, 0];     % see function delkep
% Orbit state vector:
stat = kepel_statvec(kepel);

% Velocidad Orbital
GM = 3.986004418e14;
orb_period = sqrt(4*pi^2*kepel(1)^3/GM);
omeg_0 = 2*pi/orb_period;

% Attitude elements in Euler angles of a 3-1-3 (z-x-z) rotation
eulzxz = [30, 50, 20]'*pi/180;              % converted from degrees to radians

% Attitude in quaternions
quat = ezxzquat(eulzxz);                    % converted from Euler angles
quat = quat/norm(quat);

% Floquet
floquet = 0;

% Angular velocity vector in body frame:    
w_ang = (1-floquet)*[10*0, 10*0, 10*0]'*pi/180 + floquet * [0;0;1]*pi/180;           % in radians/sec

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
tstep = 1;               % step time (sec)
tend = 12*orb_period;    % end time (10 minutes)

% Inertia matrix of axis-symetric rigid body:
iner = [0.005 0 0; 0 0.005 0; 0 0 0.005];           % in kg*m*m
iner = [0.0021 0.000001 -0.000001; 0.000001 0.0018 -0.0000001; -0.000001 -0.0000001 0.002];         % in kg*m*m
%iner = [27 0 0; 0 17 0; 0 0 25];       % in kg*m*m
factorI = 5/2.2;

% Inverse inertia matrix:
invin = inv(iner);

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
maxmagmom = 0.02*5;     % Tope momento magnético en Am2
vdq = [0;0;0];
vdqs = [0;0;0];
vdw = [0;0;0];
u = [0;0;0];
vu = u;
signq4 = 0;
vext_torq = u;
vmag_mom = u;
exu = u;

earthradius = 6371000;

xhatkmi     = zeros(6,1);
xhatk       = zeros(6,1);
vxhat       = [xhatk];
Pk          = [[eye(3)*0.00001^2 zeros(3)];[zeros(3) eye(3)*0.001^2]];
Pk          = Pk*2000;
Akm         = zeros(6);
Gamk        = tstep*eye(6);
Qk          = [[eye(3)*0.0001^2 zeros(3)];[zeros(3) eye(3)*0.001^2]];
Qmin        = Qk;
alphaQ      = 0.001;

vQk         = [];

wref        = zeros(3,1);

simularmm   = 1;
compensarmm = 1;
ekfullmodel = 1;
partialomeg = 1;
ekfrungekutta4 = 0; % No se justifica para el paso de un segundo

exmag_mom   = zeros(3,1);
mombias     = maxmagmom*0.003*[-1;0;1];
%mombias     = maxmagmom*0.003*[-0.5945;0.5229;-0.6108]*sqrt(2);
vdP         = zeros(6,1);

tnadirpointing      = 0;
nadirpointing       = floquet;
factorcompensado    = 1;

Phik=zeros(6);
exPhik=zeros(6);
exexPhik=zeros(6);
exexexPhik=zeros(6);
exexexexPhik=zeros(6);
exexexexexPhik=zeros(6);
vsigmaM     = [0;0;0;0;0;0];
vcondM      = [0;0;0;0;0;0];
H           = [eye(3) zeros(3)];

R           = eye(3)*0.1^2;
Rmin        = R;
alphaR      = 0.;
vRk         = [];
vunovM      = [];
vanglexv    = [];
vangles_sun = [];
vangles_nad = [];
angles_sun  = 0;
angles_nad  = 0;

V           = eye(6);

eclipsereal = 0;

vmagnetic   = [];

vsgam       = [];
vsrptorque  = [];
vrmmtorque  = [];
vcrmmtorque = [];
vggtorque   = [];
vaerotorque = [];
vmagtorque  = [];

vlatlon     = [];
vqos        = [1;0;0;0];
vqss        = [1;0;0;0];
vsunb       = [0;0;1];
vnadirb     = [0;0;-1];

% Attitude and orbit propagation
for t = tstart:tstep:tend

    % Orbit propagation
    kep2 = kepel + delk*t;

    % To convert from keplerian elements to state vector (if needed)
    stat = kepel_statvec(kep2);

    % Perturbation torques:
    A    = quatrmx(quat);
    exb  = A(:,1);
    ggtorque = 3*omeg_0^2*cross(exb,iner*exb);

    % RTN Frame
    Radial = stat(1:3)'/norm(stat(1:3));
    Normal = cross(stat(1:3),stat(4:6));
    Normal = Normal'/norm(Normal);
    Transverse = cross(Normal,Radial);
    Transverse = Transverse/norm(Transverse);
    NADIRP = [Normal Transverse -Radial];
    Sun = sun_dir(mjd, dfra+t);
    Xs = [0;0;1];
    Ys = cross(Xs,Sun);
    Ys = Ys/norm(Ys);
    Xs = cross(Ys,Sun);
    SUNP = [Xs Ys Sun];
    qm0 = rmxquat(SUNP);
    qm1 = rmxquat(NADIRP);
    qm2 = rmxquat(A);
    qos = quat_prod(qm1,qm2);
    qss = quat_prod(qm0,qm2);
    
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
        vlatlon = [vlatlon [colat; elong]];

        earth_field_ned = 1.e-9*igrf_field(year_frac, alt, colat, elong)';
        earth_field_ecef = ned2ecef(earth_field_ned',[sphe(2);sphe(1)]);
        earth_field_eci = terrestrial_to_inertial(gst(mjd, dfra+t), [earth_field_ecef 0 0 0])';
        earth_field_eci = earth_field_eci(1:3);
        earth_field_b = quatrmx(quat)*earth_field_eci;
        
        % Magnetic Control
        k_p = 0.00001*factorI;              % ganancia proporcional
        k_v = 0.75/factorI;                 % ganancia derivativa
        eps = 0.01;                         % epsilon

        if nadirpointing == 1
            tnadirpointing = tnadirpointing+1;
            if tnadirpointing>orb_period
                eps=2*eps; % 2*
            end;
        end;

        dq      = quat(1:3);
        dw      = w_ang + rand(3,1)*0.0000001+0*1e-06*[1;1;1];

        % Máquina de estados para apuntamiento:
        % Apuntamiento al Sol
        dqs     = cross(quatrmx(quat)*sun_dir(mjd, dfra+t),[0;0;1]);
        suni    = sun_dir(mjd, dfra+t);
        angles_sun  = asin(norm(dqs));
        if (quatrmx(quat)*suni)'*[0;0;1]>0.9 && t>1.25*orb_period && abs(-colat+pi/2-asin(suni(3)))*180/pi < 1. % && abs(stat(3))<100000
            nadirpointing=1;
        end;
        % Apuntamiento a Nadir
        if nadirpointing==1 
            dqs    = cross(quatrmx(quat)*stat(1:3)'/norm(stat(1:3)),[0;0;-1]);
            angles_nad  = asin(norm(dqs));
            %if angles_nad>45*pi/180;
            %    dqs    = 0*dqs;
            %end;
        end;
        angles      = asin(norm(dqs));
        versors     = dqs/norm(dqs);
        dqs         = sin(angles/2)*versors;
        dqs4        = cos(angles/2);
        
        nadirb      = quatrmx(quat)*stat(1:3)'/norm(stat(1:3));
        vnadirb     = [vnadirb nadirb];
        
        % SRP Torque
        sunb        = quatrmx(quat)*sun_dir(mjd, dfra+t);       % Sun vector in body frame
        Asun        = 0.1*0.1;                                  % SRP Cubesat Area
        C_refsun    = 1.3;                                      % Coefficient of Reflectivity
        Isun        = 1358/3e8*Asun*C_refsun;                   % Solar Flux intensity         1358 W/m2 (power density on Earth)  /  3e8 m/s (speed of light)
        cg_cp       = 0.01 * [1;1;0]/sqrt(2);                   % Simple model of deviation
        srptorque   = cross(sunb*Isun, cg_cp);                  % SRP torque model
        vsunb       = [vsunb sunb];

        % Aero Torque
        % https://www.spaceacademy.net.au/watch/debris/atmosmod.htm#:~:text=The%20simple%20isothermal%20model%20is,in%20the%20table%20at%20right.
        F10         = 150;
        Ap          = 0;
        T           = 900 + 2.5*(F10 - 70) + 1.5*Ap;
        hkm         = (kepel(1)-RE)/1000;
        muar        = 27 - 0.012*(hkm - 200);
        Har         = T / muar;
        rho         = (6e-10)*exp( - (hkm - 175) / Har );
        %rho        = 1.13e-13;                                 % Vallado Tabla 8-4, rango 600km - 700km
        velb        = quatrmx(quat)*stat(4:6)'/norm(stat(4:6)); % Velocity vector in body frame
        Aaero       = 0.1*0.1;                                  % Aero Cubesat Area
        C_drag      = 2.0;                                      % Coefficient of drag
        Pdyn        = 0.5*rho*norm(stat(4:6))^2;                % Dynamic Pressure
        cg_cp       = 0.01 * [1;-1;1]/sqrt(3);                  % Simple model of deviation
        aerotorque  = cross(velb*Pdyn*Aaero*C_drag, cg_cp);     % Aerodynamic torque model

        % Ángulo de yaw
        dqv     = cross(quatrmx(quat)*stat(4:6)'/norm(stat(4:6)),[0;1;0]);
        anglexv = asin(norm(dqv));
        %if tnadirpointing>3*orb_period
        %  dqs(3) = 0.001*sin(anglexv/2);
        %end;

        if (norm(dw)<0.002/4 && dqs4*dqs4>0.1 && signq4==0)
          signq4 = sign(dqs4);
        end;

        % Control Derivativo
        u       = - eps*k_v*iner*dw;

        % Eclipse
        eclipse = 0;
        if (orbit(1:3,end)'*sun_dir(mjd, dfra+t)<0)
            perpsun = (eye(3) - sun_dir(mjd, dfra+t)*sun_dir(mjd, dfra+t)')*orbit(1:3,end);
            if (norm(perpsun)<earthradius)
                eclipse = 1;
            end;
        end;
        eclipsereal = eclipse;
        if nadirpointing==1
            eclipse=0;
        end;

        % Sun pointing: término proporcional
        if (eclipse == 0 && signq4*signq4>0)
            u = u - eps*eps*k_p*inv(iner)*dqs*signq4;
        end;
        
        if (eclipse == 1 && signq4*signq4>0)
            u = 3*u;
        end;

        gamavg = gamavg*(t/(t+1)) + Skew(earth_field_b)*Skew(-earth_field_b)/(t+1)/norm(earth_field_b)/norm(earth_field_b);

        earth_field_bx = earth_field_b;

        normb2 = norm(earth_field_b)^2;

            % Filtro de Kalman del momento magnético residual
            if ekfullmodel == 0
                wref    = 0*dw; % 0*?
            else
                wref    = dw;
            end;
            Akm     = [[inv(iner)*( Skew(iner*wref) - Skew(wref)*iner + 1/normb2 * [cross(cross(earth_field_b, - partialomeg*eps*k_v*[iner(:,1)]),earth_field_b), cross(cross(earth_field_b, - partialomeg*eps*k_v*[iner(:,2)]),earth_field_b), cross(cross(earth_field_b, - partialomeg*eps*k_v*[iner(:,3)]),earth_field_b)] ) [-inv(iner)*Skew(earth_field_b)]];zeros(3,3) zeros(3)];
            
            Phik    = expm(Akm*tstep);
            Gamk    = tstep*Phik;
            if ekfullmodel == 0
              xhatkmi = Phik*xhatk + Gamk*[inv(iner)*1/normb2 * cross(cross(earth_field_b, exu),earth_field_b) + cross(xhatk(4:6), earth_field_b);zeros(3,1)];
            else
              f = [ inv(iner)*( cross(1/normb2 * cross(earth_field_b, exu) - factorcompensado*xhatk(4:6)*compensarmm, earth_field_b) + cross(xhatk(4:6), earth_field_b)                                   - cross(dw, iner*dw))                       ; zeros(3,1) ];
              k1 = tstep*f;
              f = [ inv(iner)*( cross(1/normb2 * cross(earth_field_b, 0.5*(exu+u)) - factorcompensado*(xhatk(4:6)+k1(4:6)/2)*compensarmm, earth_field_b) + cross((xhatk(4:6)+k1(4:6)/2), earth_field_b) - cross(dw+k1(1:3)/2, iner*(dw+k1(1:3)/2))) ; zeros(3,1) ];
              k2 = tstep*f;
              f = [ inv(iner)*( cross(1/normb2 * cross(earth_field_b, 0.5*(exu+u)) - factorcompensado*(xhatk(4:6)+k2(4:6)/2)*compensarmm, earth_field_b) + cross((xhatk(4:6)+k2(4:6)/2), earth_field_b) - cross(dw+k2(1:3)/2, iner*(dw+k2(1:3)/2))) ; zeros(3,1) ];
              k3 = tstep*f;
              f = [ inv(iner)*( cross(1/normb2 * cross(earth_field_b, u) -           factorcompensado*(xhatk(4:6)+k3(4:6))*compensarmm,   earth_field_b) + cross((xhatk(4:6)+k3(4:6)),   earth_field_b) - cross(dw+k3(1:3)  , iner*(dw+k3(1:3))))   ; zeros(3,1) ];
              k4 = tstep*f;
              if ekfrungekutta4==1
                xhatkmi = xhatk + (1/6)*(k1+2*k2+2*k3+k4);
              else
                xhatkmi = xhatk + k1;
              end;
            end;

            Pkmi    = Phik*Pk*Phik' + Gamk*Qk*inv(Gamk);
            Hk      = [eye(3) zeros(3)];
            %K       = Pkmi*Hk'*inv(R);
            K       = Pkmi*Hk'*inv(Hk*Pkmi*Hk'+R);
            xhatk   = xhatkmi + K*(dw-Hk*xhatkmi);
            Pk      = Pkmi - Pkmi*Hk' * inv(Hk*Pkmi*Hk'+R) * Hk*Pkmi;

            Qcorr   = K*(dw-Hk*xhatkmi)*(K*(dw-Hk*xhatkmi))';
            Qk = Qk*(1-alphaQ) + alphaQ*(Qcorr + Qmin);

            Rcorr   = (dw-Hk*xhatkmi)*(dw-Hk*xhatkmi)';
            R  = R*(1-alphaR) + alphaR*(Rcorr + Rmin);

            if ekfullmodel == 0
                exu = u + eps*k_v*iner*dw;
            else
                exu = u;
            end;

        % Compensa el momento magnético
        mag_mom = 1/normb2 * cross(earth_field_b, u) - factorcompensado*xhatk(4:6)*compensarmm;

        % pseudoSMC
        %if norm(cross(mag_mom, earth_field_b))<norm(aerotorque)
        %    mag_mom=mag_mom*norm(aerotorque)/norm(mag_mom);
        %end;
        
        if max(abs(mag_mom)) > maxmagmom             % torqrod saturation with bisection
          mag_mom = mag_mom*maxmagmom/max(abs(mag_mom));
        end

        exmag_mom = mag_mom;
       
        % Simula el momento magnético
        mag_mom = mag_mom + mombias*simularmm;
        
        controltorque = cross(mag_mom, earth_field_b);

        ext_torq = ggtorque + srptorque + aerotorque + controltorque;

        vmagtorque = [vmagtorque controltorque];

        [T, Y] = ode45('rigbody', tspan, att_vec, options, ext_torq, iner, invin, mag_mom, earth_field);
        
    end;

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
    vxhat = [vxhat xhatk];
    vdP   = [vdP [sqrt(Pkmi(1,1));sqrt(Pkmi(2,2));sqrt(Pkmi(3,3));sqrt(Pkmi(4,4));sqrt(Pkmi(5,5));sqrt(Pkmi(6,6))]];
    vQk = [vQk [min([Qk(1,1) Qk(2,2) Qk(3,3)]);min([Qk(4,4) Qk(5,5) Qk(6,6)])]];
    vRk = [vRk [R(1,1);R(2,2);R(3,3)]];
    if (vqos(:,end)'*qos<0) 
        qos=-qos; 
    end;
    vqos = [vqos qos];
    if (vqss(:,end)'*qss<0) 
        qss=-qss; 
    end;
    vqss = [vqss qss];
    
    exexexexexPhik = exexexexPhik;
    exexexexPhik = exexexPhik;
    exexexPhik  = exexPhik;
    exexPhik    = exPhik;
    exPhik      = Phik;
    O1          = [H;H*Phik];
    O2          = [H;H*exPhik;        H*Phik*exPhik];
    O3          = [H;H*exexPhik;      H*exPhik*exexPhik;            H*Phik*exPhik*exexPhik];
    O4          = [H;H*exexexPhik;    H*exexPhik*exexexPhik;        H*exPhik*exexPhik*exexexPhik;            H*Phik*exPhik*exexPhik*exexexPhik];
    O5          = [H;H*exexexexPhik;  H*exexexPhik*exexexexPhik;    H*exexPhik*exexexPhik*exexexexPhik;      H*exPhik*exexPhik*exexexPhik*exexexexPhik;          H*Phik*exPhik*exexPhik*exexexPhik*exexexexPhik];
    O6          = [H;H*exexexexexPhik;H*exexexexPhik*exexexexexPhik;H*exexexPhik*exexexexPhik*exexexexexPhik;H*exexPhik*exexexPhik*exexexexPhik*exexexexexPhik;  H*exPhik*exexPhik*exexexPhik*exexexexPhik*exexexexexPhik;H*Phik*exPhik*exexPhik*exexexPhik*exexexexPhik*exexexexexPhik];
    M1          = O1'*O1;
    M2          = O2'*O2;
    M3          = O3'*O3;
    M4          = O4'*O4;
    M5          = O5'*O5;
    M6          = O6'*O6;
    vsigmaM     = [vsigmaM [min(eig(M1));min(eig(M2));min(eig(M3));min(eig(M4));min(eig(M5));min(eig(M6))]];
    vcondM      = [vcondM [max(eig(M1))/min(eig(M1));max(eig(M2))/min(eig(M2));max(eig(M3))/min(eig(M3));max(eig(M4))/min(eig(M4));max(eig(M5))/min(eig(M5));max(eig(M6))/min(eig(M6))]];
    [V D]       = eig(M6);
    vunovM      = [vunovM V(:,1)];
    vanglexv    = [vanglexv anglexv];
    vmagnetic   = [vmagnetic earth_field_b];
    vangles_sun = [vangles_sun angles_sun];
    vangles_nad = [vangles_nad angles_nad];

    vsgam       = [vsgam svd(gamavg)];
    vsrptorque  = [vsrptorque srptorque];
    vrmmtorque  = [vrmmtorque cross(earth_field_b, mombias)];
    vcrmmtorque = [vcrmmtorque cross(earth_field_b, mombias-xhatk(4:6))];
    vggtorque   = [vggtorque ggtorque];
    vaerotorque = [vaerotorque aerotorque];

end;

% 1) Comparar estado con dispersión de estimación (diagonal de P) 
% 2) Evaluar la matriz de observabilidad de tiempo variante 
% 3) Simular efectos de bias de gyro y magnetómetro
% 4) Documentar con la función nolineal
% 5) Analizar con eclipse

close all

figure(1010);clf;
title('Singular values of \int_0^T\Gamma(t)dt/T');hold on;
plot(time(1:length(vsgam))/orb_period, vsgam(1,:),'r');hold on;
plot(time(1:length(vsgam))/orb_period, vsgam(2,:),'g');hold on;
plot(time(1:length(vsgam))/orb_period, vsgam(3,:),'b');hold on;
grid on;
xlabel('Time (orbits)');

figure(101010);clf;
title('Magnitude of Disturbance and Control Torques [Nm, Logarithmic]');hold on;
plot(time(1:length(vsgam))/orb_period, log(sqrt(vrmmtorque(1,:).^2+vrmmtorque(2,:).^2+vrmmtorque(3,:).^2)),'r');hold on;
plot(time(1:length(vsgam))/orb_period, log(sqrt(vcrmmtorque(1,:).^2+vcrmmtorque(2,:).^2+vcrmmtorque(3,:).^2)),'k');hold on;
plot(time(1:length(vsgam))/orb_period, log(sqrt(vggtorque(1,:).^2+vggtorque(2,:).^2+vggtorque(3,:).^2)),'g');hold on;
plot(time(1:length(vsgam))/orb_period, log(sqrt(vsrptorque(1,:).^2+vsrptorque(2,:).^2+vsrptorque(3,:).^2)),'b');hold on;grid on;
plot(time(1:length(vsgam))/orb_period, log(sqrt(vaerotorque(1,:).^2+vaerotorque(2,:).^2+vaerotorque(3,:).^2)),'m');hold on;grid on;
%plot(time(1:length(vsgam))/orb_period, log(sqrt(vmagtorque(1,:).^2+vmagtorque(2,:).^2+vmagtorque(3,:).^2)),'c');hold on;grid on;
%legend('RMM','RMM compensation error','Gravity Gradient','Solar Radiation Pressure','Aerodynamic','Control','Location','NorthEast')
legend('RMM','RMM compensation error','Gravity Gradient','Solar Radiation Pressure','Aerodynamic','Location','NorthEast')
xlabel('Time (orbits)');


figure(101011);clf;
title('Magnitude of Disturbance Torques [Nm]');hold on;
plot(time(1:length(vsgam))/orb_period, (sqrt(vrmmtorque(1,:).^2+vrmmtorque(2,:).^2+vrmmtorque(3,:).^2)),'r');hold on;
plot(time(1:length(vsgam))/orb_period, (sqrt(vcrmmtorque(1,:).^2+vcrmmtorque(2,:).^2+vcrmmtorque(3,:).^2)),'k');hold on;
plot(time(1:length(vsgam))/orb_period, (sqrt(vggtorque(1,:).^2+vggtorque(2,:).^2+vggtorque(3,:).^2)),'g');hold on;
plot(time(1:length(vsgam))/orb_period, (sqrt(vsrptorque(1,:).^2+vsrptorque(2,:).^2+vsrptorque(3,:).^2)),'b');hold on;grid on;
plot(time(1:length(vsgam))/orb_period, (sqrt(vaerotorque(1,:).^2+vaerotorque(2,:).^2+vaerotorque(3,:).^2)),'m');hold on;grid on;
plot(time(1:length(vsgam))/orb_period, (sqrt(vmagtorque(1,:).^2+vmagtorque(2,:).^2+vmagtorque(3,:).^2)),'c');hold on;grid on;
legend('RMM','RMM compensation error','Gravity Gradient','Solar Radiation Pressure','Aerodynamic','Control','Location','NorthEast')
xlabel('Time (orbits)');

% Output visualization
figure(1)
plot(time/orb_period, vdq(1,:),'r');hold on;
plot(time/orb_period, vdq(2,:),'g');hold on;
plot(time/orb_period, vdq(3,:),'b');hold on;
xlabel('Time (orbits)')
ylabel('Quaternion 1-3')
title('Attitude in quaternion vector');
grid on;

inii=10000;
figure(11);clf;
plot3(vdq(1,inii:end),vdq(2,inii:end),vdq(3,inii:end),'k');hold on;
xlabel('q_1');
ylabel('q_2');
zlabel('q_3');
title('Full quaternion vector on Sun Pointing (two days)');
grid on;


figure(10)
plot(time/orb_period, vdqs(1,:),'r');hold on;
plot(time/orb_period, vdqs(2,:),'g');hold on;
plot(time/orb_period, vdqs(3,:),'b');hold on;
xlabel('Time (orbits)')
ylabel('Quaternion 1-3')
title('Attitude in partial quaternion vector');
grid on;


in0=24800*0+1;
figure(111);clf;
plot3(vqos(1,in0:end),vqos(2,in0:end),vqos(3,in0:end),'k');hold on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Spacecraft attitude quaternion relative to orbital frame');
grid on;


figure(1110);clf;
subplot(2,1,1);
plot(time/orb_period, vqos(1,:),'r');hold on;
plot(time/orb_period, vqos(2,:),'g');hold on;
plot(time/orb_period, vqos(3,:),'b');hold on;
%plot(time/orb_period, vqos(4,:),'k');hold on;
xlabel('Time (orbits)')
ylabel('Quaternion qos')
title('Attitude relative to an orbital frame');
grid on;
subplot(2,1,2);
plot(time/orb_period, vqss(1,:),'r');hold on;
plot(time/orb_period, vqss(2,:),'g');hold on;
plot(time/orb_period, vqss(3,:),'b');hold on;
%plot(time/orb_period, vqss(4,:),'k');hold on;
xlabel('Time (orbits)')
ylabel('Quaternion qss')
title('Attitude relative to a Sun pointing frame');
grid on;



figure(11101);clf;
subplot(2,1,1);
plot(time/orb_period, vqos522n(1,:),'r');hold on;
plot(time/orb_period, vqos522n(2,:),'g');hold on;
plot(time/orb_period, vqos522n(3,:),'b');hold on;
plot(time(1:length(vqos521n))/orb_period, vqos521n(1,:),'r--');hold on;
plot(time(1:length(vqos521n))/orb_period, vqos521n(2,:),'b--');hold on;
plot(time(1:length(vqos521n))/orb_period, vqos521n(3,:),'b--');hold on;
%plot(time/orb_period, vqos(4,:),'k');hold on;
xlabel('Time (orbits)')
ylabel('Quaternion qos')
title('Sensitivity to initial condition of the attitude in nadir pointing');
grid on;
subplot(2,1,2);
plot(time/orb_period, vqss522(1,:),'r');hold on;
plot(time/orb_period, vqss522(2,:),'g');hold on;
plot(time/orb_period, vqss522(3,:),'b');hold on;
plot(time(1:length(omeg521))/orb_period, vqss521(1,:),'r--');hold on;
plot(time(1:length(omeg521))/orb_period, vqss521(2,:),'g--');hold on;
plot(time(1:length(omeg521))/orb_period, vqss521(3,:),'b--');hold on;
%plot(time/orb_period, vqss(4,:),'k');hold on;
xlabel('Time (orbits)')
ylabel('Quaternion qss')
title('Sensitivity to initial condition of the attitude in Sun pointing');
grid on;

%omeg521n=omeg;
%vqss521n=vqss;

%omeg522n=omeg;
%vqos522n=vqos;


figure(1111);clf;
nadir=30000;
title('Partial Error vs. Latitude');hold on;
plot(vlatlon(1,nadir:end)*180/pi, sqrt(vdqs(1,(1+nadir):end).^2+vdqs(2,(1+nadir):end).^2)*180/pi*2.,'k');hold on;
grid on;
ylabel('Nadir Pointing Error [deg]');
xlabel('Latitude [deg]');

figure(10101010);
plot(vlatlon(1,nadir:end)*180/pi-90, vlatlon(2,nadir:end)*180/pi,'k');hold on;

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


figure(211)
subplot(2,1,1);
plot(time/orb_period, vsunb(1,:),'r');hold on;
plot(time/orb_period, vsunb(2,:),'g');hold on;
plot(time/orb_period, vsunb(3,:),'b');hold on;
hold on;
xlabel('Time (orbits)')
ylabel('Sun vector components')
title('Sun vector in body frame')
grid on;
subplot(2,1,2);
plot(time/orb_period, vnadirb(1,:),'r');hold on;
plot(time/orb_period, vnadirb(2,:),'g');hold on;
plot(time/orb_period, vnadirb(3,:),'b');hold on;
hold on;
xlabel('Time (orbits)')
ylabel('Nadir vector components')
title('Nadir vector in body frame')
grid on;



figure(22);clf;
plot(time/orb_period,vmag_mom(1,:),'r');hold on;
plot(time/orb_period,vmag_mom(2,:),'g');hold on;
plot(time/orb_period,vmag_mom(3,:),'b');hold on;
plot(time/orb_period,maxmagmom*ones(size(vext_torq(3,:))),'k--');hold on;
plot(time/orb_period,-maxmagmom*ones(size(vext_torq(3,:))),'k--');hold on;
title('Magnetorquers Control Magnetic Moment [Nm]')
grid on;

figure(999);clf;
subplot(3,1,1);
plot(time/orb_period,factorcompensado*vxhat(4,:)-mombias(1)*ones(1,length(time)),'r');hold on;
plot(time/orb_period,factorcompensado*vxhat(4,:)+vdP(4,:),'k');hold on;
plot(time/orb_period,factorcompensado*vxhat(4,:)-vdP(4,:),'k');hold on;
title('X Magnetic Moment Residual Estimation Error');
subplot(3,1,2);
plot(time/orb_period,factorcompensado*vxhat(5,:)-mombias(2)*ones(1,length(time)),'g');hold on;
plot(time/orb_period,factorcompensado*vxhat(5,:)+vdP(5,:),'k');hold on;
plot(time/orb_period,factorcompensado*vxhat(5,:)-vdP(5,:),'k');hold on;
title('Y Magnetic Moment Residual Estimation Error');
subplot(3,1,3);
plot(time/orb_period,factorcompensado*vxhat(6,:)-mombias(3)*ones(1,length(time)),'b');hold on;
plot(time/orb_period,factorcompensado*vxhat(6,:)+vdP(6,:),'k');hold on;
plot(time/orb_period,factorcompensado*vxhat(6,:)-vdP(6,:),'k');hold on;
title('Z Magnetic Moment Residual Estimation Error');
grid on;


figure(9991);clf;
subplot(2,1,1);
plot(time/orb_period,factorcompensado*vxhat(4,:),'r');hold on;
plot(time/orb_period,factorcompensado*vxhat(5,:),'g');hold on;
plot(time/orb_period,factorcompensado*vxhat(6,:),'b');hold on;
axis([0 10 -0.0006 0.0006]);
title('Magnetic Moment Residual Estimation');hold on;
xlabel('Time [orbits]');
grid on;
subplot(2,1,2);
plot(time/orb_period,vdP(4,:),'r');hold on;
plot(time/orb_period,vdP(5,:),'g');hold on;
plot(time/orb_period,vdP(6,:),'b');hold on;
title('Magnetic Moment Residual Estimation Covariance');hold on;
xlabel('Time [orbits]');
grid on;


figure(9999);clf;
subplot(3,1,1);
plot(time/orb_period,vxhat(1,:)-omeg(1,:),'r');hold on;
plot(time/orb_period,vxhat(1,:)-omeg(1,:)+vdP(1,:),'k');hold on;
plot(time/orb_period,vxhat(1,:)-omeg(1,:)-vdP(1,:),'K');hold on;
title('X Angular Velocity Estimation Error');
subplot(3,1,2);
plot(time/orb_period,vxhat(2,:)-omeg(2,:),'g');hold on;
plot(time/orb_period,vxhat(2,:)-omeg(2,:)+vdP(2,:),'k');hold on;
plot(time/orb_period,vxhat(2,:)-omeg(2,:)-vdP(2,:),'k');hold on;
title('Y Angular Velocity Estimation Error');
subplot(3,1,3);
plot(time/orb_period,vxhat(3,:)-omeg(3,:),'b');hold on;
plot(time/orb_period,vxhat(3,:)-omeg(3,:)+vdP(3,:),'k');hold on;
plot(time/orb_period,vxhat(3,:)-omeg(3,:)-vdP(3,:),'k');hold on;
title('Z Angular Velocity Estimation Error');
grid on;


figure(99991);clf;
subplot(2,1,1);
plot(time(1:(end-1))/orb_period,vQk(1,:),'k');hold on;
title('Adaptive process noise covariance for \omega');hold on;
xlabel('Time {orbits]');
grid on;
subplot(2,1,2);
plot(time(1:(end-1))/orb_period,vQk(2,:),'k');hold on;
title('Adaptive process noise covariance for m_o');hold on;
xlabel('Time [orbits]');
grid on;
xlabel('Time {orbits]');

figure(99992);clf;
title('Log of min singular value of O_{k,1} (red), O_{k,2} (green), O_{k,3} (blue), O_{k,4} (mag), O_{k,5} (cyan) and O_{k,6} (black)');hold on;
plot(time/orb_period,log(sqrt(vsigmaM(1,:))),'r');hold on;
plot(time/orb_period,log(sqrt(vsigmaM(2,:))),'g');hold on;
plot(time/orb_period,log(sqrt(vsigmaM(3,:))),'b');hold on;
plot(time/orb_period,log(sqrt(vsigmaM(4,:))),'m');hold on;
plot(time/orb_period,log(sqrt(vsigmaM(5,:))),'c');hold on;
plot(time/orb_period,log(sqrt(vsigmaM(6,:))),'k');hold on;
axis([0 10 -26 -4]);
grid on;
xlabel('Time [orbits]');

figure(99993);clf;
title('Log of condition number of M_1 (r), M_2 (g), M_3 (b), M_4 (m), M_5 (c) and M_6 (k)');hold on;
plot(time/orb_period,log(vcondM(1,(1:length(time)))),'r');hold on;
plot(time/orb_period,log(vcondM(2,(1:length(time)))),'g');hold on;
plot(time/orb_period,log(vcondM(3,(1:length(time)))),'b');hold on;
plot(time/orb_period,log(vcondM(4,(1:length(time)))),'m');hold on;
plot(time/orb_period,log(vcondM(5,(1:length(time)))),'c');hold on;
plot(time/orb_period,log(vcondM(6,(1:length(time)))),'k');hold on;
axis([0 10 0 50]);
grid on;
xlabel('Time [orbits]');

figure(99994);clf;
title('Adapted measurement covariance');hold on;
plot(time(1:length(vRk))/orb_period,vRk(1,:),'r');hold on;
plot(time(1:length(vRk))/orb_period,vRk(2,:),'g');hold on;
plot(time(1:length(vRk))/orb_period,vRk(3,:),'b');hold on;
grid on;
xlabel('Time [orbits]');

figure(99995);clf;
title('Least observable eigenvector components {\hat m}_x, {\hat m}_y, {\hat m}_z');hold on;
plot(time(1:length(vRk))/orb_period,vunovM(4,:),'r');hold on;
plot(time(1:length(vRk))/orb_period,vunovM(5,:),'g');hold on;
plot(time(1:length(vRk))/orb_period,vunovM(6,:),'b');hold on;
grid on;
xlabel('Time [orbits]');

figure(99996);clf;
hold on;
title('\beta_M: angle between the RMM and the least observable eigenvector of M_{k,6}');hold on;
plot(time(1:length(vRk))/orb_period,(180/pi)*acos((vunovM(4,:)*mombias(1)+vunovM(5,:)*mombias(2)+vunovM(6,:)*mombias(3))/norm(mombias)),'k');hold on;
ylabel('[deg]');
grid on;
xlabel('Time [orbits]');


minle=min([length(exvanglexv621) length(exvanglexv620)]);

figure(99997);clf;
title('Angle between +Y_b and Transverse axis');hold on;
plot(time(1:minle)/orb_period,2*(180/pi)*exvanglexv620(1:minle),'k');hold on;
plot(time(1:minle)/orb_period,2*(180/pi)*exvanglexv621(1:minle),'r');hold on;
ylabel('[deg]');
grid on;


minle=min([length(exvanglexv621) length(exvanglexv620) length(time)]);
Lexp = log(180/pi*abs(exvanglexv621(1:minle)-exvanglexv620(1:minle)))./((1:minle)/orb_period);
figure(999971);clf;
subplot(3,1,3);
plot(time(1:minle)/orb_period,Lexp(1:minle),'r');grid on;
axis([1 10 -0.1 3]);
minle=min([length(exvanglexv521) length(exvanglexv522) length(time)]);
Lexp = log(180/pi*abs(exvanglexv521(1:minle)-exvanglexv522(1:minle)))./((1:minle)/orb_period);
subplot(3,1,2);hold on;
plot(time(1:minle)/orb_period,Lexp(1:minle),'r');grid on;
axis([1 10 -0.1 3]);
minle=min([length(exvanglexv721) length(exvanglexv722) length(time)]);
Lexp = log(180/pi*abs(exvanglexv721(1:minle)-exvanglexv722(1:minle)))./((1:minle)/orb_period);
subplot(3,1,1);hold on;
plot(time(1:minle)/orb_period,Lexp(1:minle),'r');grid on;
axis([1 10 -0.1 3]);


figure(99997);clf;
hold on;
title('Angles +Z/Sun (r) and -Z/Nadir (g) +Y_b/Transverse (b)');hold on;
plot(time(1:length(vanglexv))/orb_period,2*(180/pi)*vanglexv,'b');hold on;
plot(time(1:length(vanglexv))/orb_period,2*(180/pi)*vangles_sun,'r');hold on;
plot(time(1:length(vanglexv))/orb_period,2*(180/pi)*vangles_nad,'g');hold on;
ylabel('[deg]');
grid on;
xlabel('Time [orbits]');


figure(99998);clf;
hold on;
subplot(4,1,1);
title('Angular velocity norm |\omega|');hold on;
plot(time(1:length(omeg))/orb_period,(180/pi)*sqrt(omeg(1,:).^2+omeg(2,:).^2+omeg(3,:).^2),'k');hold on;
ylabel('[deg/sec]');
grid on;
xlabel('Time [orbits]');
subplot(4,1,2);
title('Angle beween +Z_b and Sun vector in body frame');hold on;
plot(time(1:length(vanglexv))/orb_period,2*(180/pi)*vangles_sun,'r');hold on;
ylabel('[deg]');
grid on;
xlabel('Time [orbits]');
subplot(4,1,3);
title('Angle between -Z_b vs. Orbit Nadir in body frame');hold on;
plot(time(1:length(vanglexv))/orb_period,2*(180/pi)*vangles_nad,'g');hold on;
ylabel('[deg]');
grid on;
xlabel('Time [orbits]');
subplot(4,1,4);
title('Angle between +Y_b and Orbit Transverse in body frame');hold on;
plot(time(1:length(vanglexv))/orb_period,2*(180/pi)*vanglexv,'b');hold on;
ylabel('[deg]');
grid on;
xlabel('Time [orbits]');


%X6 = [vdqs;omeg;vxhat];
%Xf6 = X6(:,1:ceil(orb_period):end);
%save('flo6.mat','X6','Xf6');

%[X1,Xf1]=load('flo1.mat');

if 1==0

corb_period = ceil(orb_period);
XX00 = [X1(1:6,2) X2(1:6,2) X3(1:6,2) X4(1:6,2) X5(1:6,2) X6(1:6,2)];
XX0  = [Xf1(1:6,1) Xf2(1:6,1) Xf3(1:6,1) Xf4(1:6,1) Xf5(1:6,1) Xf6(1:6,1)];
XXT  = [Xf1(1:6,2) Xf2(1:6,2) Xf3(1:6,2) Xf4(1:6,2) Xf5(1:6,2) Xf6(1:6,2)];
XX2T = [Xf1(1:6,3) Xf2(1:6,3) Xf3(1:6,3) Xf4(1:6,3) Xf5(1:6,3) Xf6(1:6,3)];
XX3T = [Xf1(1:6,4) Xf2(1:6,4) Xf3(1:6,4) Xf4(1:6,4) Xf5(1:6,4) Xf6(1:6,4)];
XX4T = [Xf1(1:6,5) Xf2(1:6,5) Xf3(1:6,5) Xf4(1:6,5) Xf5(1:6,5) Xf6(1:6,5)];
XX5T = [Xf1(1:6,6) Xf2(1:6,6) Xf3(1:6,6) Xf4(1:6,6) Xf5(1:6,6) Xf6(1:6,6)];
XX6T = [Xf1(1:6,7) Xf2(1:6,7) Xf3(1:6,7) Xf4(1:6,7) Xf5(1:6,7) Xf6(1:6,7)];
XX15T = [Xf1(1:6,14) Xf2(1:6,14) Xf3(1:6,14) Xf4(1:6,14) Xf5(1:6,14) Xf6(1:6,14)];

B = inv(XX00([1 2 4 5 6],1:5))*XXT([1 2 4 5 6],1:5);
mul_floquet = eig(B);
exp_floquetoneday = log(mul_floquet)/orb_period

B = inv(XX0([1 2 4 5 6],1:5))*XXT([1 2 4 5 6],1:5);
mul_floquet = eig(B);
exp_floquetdetumbling = log(mul_floquet)/orb_period

B = inv(XXT([1 2 4 5 6],1:5))*XX2T([1 2 4 5 6],1:5);
mul_floquet = eig(B);
exp_floquetsunpointing = log(mul_floquet)/orb_period

B = inv(XX5T([1 2 4 5],[1 2 4 5]))*XX6T([1 2 4 5],[1 2 4 5]);
mul_floquet = eig(B);
exp_floquetnadirpointing = log(mul_floquet)/orb_period


figure(1010);clf;
xlabel('Time [Orbits]');hold on;
ylabel('Am^2');hold on;
title('RMM estimation error on x (r), y(g) and z(b) axes');
plot((1:length(X1))/orb_period,X1(10,:)-mombias(1)*ones(1,length(X1)),'r');hold on;
plot((1:length(X2))/orb_period,X2(10,:)-mombias(1)*ones(1,length(X1)),'r');hold on;
plot((1:length(X3))/orb_period,X3(10,:)-mombias(1)*ones(1,length(X1)),'r');hold on;
plot((1:length(X4))/orb_period,X4(10,:)-mombias(1)*ones(1,length(X1)),'r');hold on;
plot((1:length(X5))/orb_period,X5(10,:)-mombias(1)*ones(1,length(X1)),'r');hold on;
plot((1:length(X6))/orb_period,X6(10,:)-mombias(1)*ones(1,length(X1)),'r');hold on;
plot((1:length(X1))/orb_period,X1(11,:)-mombias(2)*ones(1,length(X1)),'g');hold on;
plot((1:length(X2))/orb_period,X2(11,:)-mombias(2)*ones(1,length(X1)),'g');hold on;
plot((1:length(X3))/orb_period,X3(11,:)-mombias(2)*ones(1,length(X1)),'g');hold on;
plot((1:length(X4))/orb_period,X4(11,:)-mombias(2)*ones(1,length(X1)),'g');hold on;
plot((1:length(X5))/orb_period,X5(11,:)-mombias(2)*ones(1,length(X1)),'g');hold on;
plot((1:length(X6))/orb_period,X6(11,:)-mombias(2)*ones(1,length(X1)),'g');hold on;
plot((1:length(X1))/orb_period,X1(12,:)-mombias(3)*ones(1,length(X1)),'b');hold on;
plot((1:length(X2))/orb_period,X2(12,:)-mombias(3)*ones(1,length(X1)),'b');hold on;
plot((1:length(X3))/orb_period,X3(12,:)-mombias(3)*ones(1,length(X1)),'b');hold on;
plot((1:length(X4))/orb_period,X4(12,:)-mombias(3)*ones(1,length(X1)),'b');hold on;
plot((1:length(X5))/orb_period,X5(12,:)-mombias(3)*ones(1,length(X1)),'b');hold on;
plot((1:length(X6))/orb_period,X6(12,:)-mombias(3)*ones(1,length(X1)),'b');hold on;
axis([0 15 -0.0002 0.0002]);
grid on;

end;

exvanglexv722=vanglexv;