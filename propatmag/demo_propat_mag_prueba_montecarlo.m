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

ECLIPSE = 1;
RMM = 1;
GRAV_GRAD = 1;
SOLAR_TORQ = 1;
DRAG = 1;
GRAPH_PERTURBATIONS = 1;

J_DIAGONAL = 0;
SUN_POINTING = 1;
NADIR_POINTING = 0;

RMM_ESTIMATE = 1;
RMM_COMPENSATE = 1;
Q_ESTIMATE = 1;
GRAPH_ESTIMATES = 1;

MONTECARLO_ITERATIONS = 100;

RE = 6378000;
kepel = [RE + 400000, 0.01, 98*pi/180, 0, 0, 0];
num_orbits = 10;

info = [kepel, num_orbits];

mkdir("./mc", timestamp)

writematrix(info, "mc/" + timestamp + "/dq.csv")
writematrix(info, "mc/" + timestamp + "/dqs.csv")
writematrix(info, "mc/" + timestamp + "/dw.csv")
writematrix(info, "mc/" + timestamp + "/mag_mom.csv")
writematrix(info, "mc/" + timestamp + "/ext_torq.csv")
writematrix(info, "mc/" + timestamp + "/rmm_hat.csv")

for i=0:MONTECARLO_ITERATIONS

    %-------------------------------ORBITA-------------------------------------
    
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
    eulzxz = [0, 0, 0]'*pi/180;   % converted from degrees to radians
    
    % Attitude in quaternions
    quat = ezxzquat(eulzxz);        % converted from Euler angles
    quat = quat/norm(quat);
    
    % Angular velocity vector in body frame:
    w_ang = [normrnd(0, 3.3), normrnd(0, 3.3), normrnd(0, 3.3)]'*pi/180;           % in radians/sec
    
    % Initial control torque:
    contq = [0 0 0]';
    angles = 0;
    angles_ = 0;
    
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
    tend = num_orbits*orb_period;    % end time (10 minutes)
    
    %-----------------------------DINAMICA-------------------------------------
    
    % Inertia matrix of axis-symetric rigid body:
    iner = [0.0021 1e-6 -1e-6; 1e-6 0.0018 -1e-7; -1e-6 -1e-7 0.002];
    if J_DIAGONAL
        iner = [0.0021 0 0; 0 0.0018 0; 0 0 0.002];
    end
             % in kg*m*m
    
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
    
    %---------------------------DRAG ATMOSFERICO------------------------------
    
    C_d = 2.2; %Coeficiente de rozamiento
    p_atmos = 3.8e-12; %Densidad atmosferica
    
    %-----------------------VECTORES PARA GRAFICAR-----------------------------
    
    % Initial vectors
    time = tstart;          % to store time
    euler = eulzxz*180/pi;  % Euler angles
    omeg = w_ang;           % Angular velocity
    orbit = stat';          % Orbit elements (state vector)
    keorb = kepel';         % Orbit elements (keplerian)
    
    gamma = 0*eye(3);
    gamma_acum = 0*eye(3);
    gamavg = 0*eye(3);      % Gamma avg inicial
    maxmagmom = 0.1;       % Tope momento magnético en Am2
    vdq = [0;0;0];
    vdqs = [0;0;0];
    vdw = w_ang;
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
    vsun_torq = [0; 0; 0];      
    vdrag_torq = [0; 0; 0];
    vgrav_grad = [0; 0; 0];
    vrmm_torq = [0; 0; 0];
    vdqs_true = [0; 0; 0];
    vgamma_avas = [0; 0; 0];
    %------------------------------SIMULACION----------------------------------
    
    for t = tstart:tstep:tend
    %for t = 1:tstep
    
        % Orbit propagation
        kep2 = kepel + delk*t;
        % To convert from keplerian elements to state vector (if needed)
        stat = kepel_statvec(kep2);
    
        %-------------------------PERTURBACIONES-------------------------------
        grad_grav = [0; 0; 0];
        if GRAV_GRAD == 1
            A    = quatrmx(quat);
            exb  = A(:,1);
            grad_grav = 3*omeg_0^2*cross(exb,iner*exb);
        end
        
        mom_res = [0; 0; 0];
        if RMM == 1
            mom_res = maxmagmom*[0.005; 0.005; -0.005];
        end
        
        %Peor caso
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
    
            %--------------------------GAMMA-----------------------------------
    
            gamavg = gamavg*(t/(t+tstep)) + Skew(earth_field_b)*Skew(-earth_field_b)/(t+tstep)/norm(earth_field_b)/norm(earth_field_b);
    
            gamma_avas = eig(gamavg);
            
            %--------------------------ECLIPSE---------------------------------
    
            
            eclipse = 0;
            if (orbit(1:3,end)'*sun_dir(mjd, dfra+t)<0 && ECLIPSE == 1)
                perpsun = (eye(3) - sun_dir(mjd, dfra+t)*sun_dir(mjd, dfra+t)')*orbit(1:3,end);
                if (norm(perpsun)<earthradius)
                    eclipse = 1;
                end
            end
    
            %-------------------------HORIZONTE--------------------------------
    
            no_horizon = 0;
            if angles_ > pi/4 || angles_ < -pi/4
                no_horizon = 1;
            end
    
            %---------------------CONTROL MAGNETICO----------------------------
            
    
            dq = quat(1:3);
            dw = w_ang;
            
            % Error en base al modo
            if SUN_POINTING  %Se configura que apunte solo al sol, dqs se actualiza siempre
                dqs = cross(quatrmx(quat)*sun_dir(mjd, dfra+t),[0;0;1]);
            end
            if NADIR_POINTING && no_horizon==0 %Se configura apuntamiento solo a nadir, se actualiza solo si se ve el horizonte
                dqs = cross(quatrmx(quat)*stat(1:3)'/norm(stat(1:3)),[0;0;-1]);
            end
            if SUN_POINTING==0 && NADIR_POINTING==0 %No hay apuntamiento, solo detumbling
                dqs = [0; 0; 0];
            end
            angles = asin(norm(dqs));
            versors = dqs/norm(dqs);
            dqs = sin(angles/2)*versors;
            dqs4 = cos(angles/2); 
    
    
            if SUN_POINTING % Se calculan los errores verdaderos
                dqs_true = cross(quatrmx(quat)*sun_dir(mjd, dfra+t),[0;0;1]);
            end
            if NADIR_POINTING
                dqs_true = cross(quatrmx(quat)*stat(1:3)'/norm(stat(1:3)),[0;0;-1]);
            end
            if SUN_POINTING==0 && NADIR_POINTING==0
                dqs_true = [0; 0; 0];
            end
            angles_ = asin(norm(dqs_true));
            versors_ = dqs_true/norm(dqs_true);
            dqs_true = sin(angles_/2)*versors_; 
           
    
            % Marcar fin del detumbling
            if (norm(dw)<0.0005 && dqs4*dqs4>0.1 && signq4==0) 
               signq4=sign(dqs4);
            end
    
            % Acualizar ganancias según el modo
            if SUN_POINTING && signq4*signq4>0
                k_p = 0.000025;% ganancia proporcional
                k_v = 0.3; % ganancia derivativa
                if eclipse 
                    eps = 0.03;% epsilon
                else
                    eps = 0.01;% epsilon
                end
            end 
            if NADIR_POINTING && signq4*signq4>0
                if no_horizon
                    k_p = 0.0000025;% ganancia proporcional
                    k_v = 0.4; % ganancia derivativa
                    eps = 0.01;%epsilon
                else
                    k_p = 0.000025;% ganancia proporcional
                    k_v = 0.3; % ganancia derivativa
                    eps = 0.02;%epsilon
                end
            end
            controller = update(controller, k_v, k_p, eps);
    
            % Determinar apuntamiento
            pointing = 0;
            if (((eclipse == 0 && SUN_POINTING) || (NADIR_POINTING)) && signq4*signq4>0)
                pointing = 1;
            end
            
            % Generar acción de control
            u = get_control_action(controller, dqs, signq4, dw, pointing);
            
            %----------------ESTIMACION DE MOMENTO RESIDUAL--------------------
            phik = 0;
            rmm_avas = [0;0;0];
            if RMM_ESTIMATE == 1
                %[x_pred, sigma, phikm1, Q] = ekf_rmm(x_pred, (mag_mom - mom_res), iner, earth_field_b, sigma, dw, Q, Q_min, alpha, R, tstep);
                [ekf_rmm, phik] = update(ekf_rmm, mag_mom - mom_res, earth_field_b, dw);
                [x_pred, sigma, Q] = get_estimates(ekf_rmm);
                rmm_diag_cov = diag(sigma);
            end
    
            %---------------------MOMENTO MAGNÉTICO----------------------------
            
            normb2 = norm(earth_field)^2;
    
            mag_mom = 1/normb2 * cross(earth_field_b, u) + mom_res - RMM_COMPENSATE*x_pred(4:6);
            
            
            if max(abs(mag_mom)) > maxmagmom             % torqrod saturation with bisection
              mag_mom = mag_mom*maxmagmom/max(abs(mag_mom));
            end
            
            ext_torq = ext_torq + cross(mag_mom, earth_field_b);
            
            %---------------------MOMENTO SOLAR----------------------------
            
            sun_torq = [0 0 0];
            if (SOLAR_TORQ == 1 && eclipse == 0)
                s_eci = sun_dir(mjd, dfra + t); % sun vector en ECI
                v_sun_body = quatrmx(quat) * s_eci(:); %sun vector en body
                norm_sun = norm(v_sun_body);
                r_presion = (v_sun_body/norm_sun)*0.001;
                [F_sun, sun_torq] = solar_pressure(v_sun_body, T_sun, S_sat, r_presion', C_r);
            end
            ext_torq = ext_torq + sun_torq';
            
            %-------------------------DRAG ATMOSFERICO-------------------------
            drag_torq = [0 0 0];
            if (DRAG == 1)
                v_inercial = stat(4:end);
                v_rela = quatrmx(quat)*v_inercial';
                [F_drag, drag_torq] = atmosferic_drag(r_cp, v_rela, S_sat, C_d, p_atmos);
            end
            ext_torq = ext_torq + drag_torq';
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
        vdqs_true = [vdqs_true dqs_true];
        vdw = [vdw dw]; 
        vmag_mom = [vmag_mom mag_mom];
        vext_torq = [vext_torq ext_torq];
        if RMM_ESTIMATE
            vrmm_hat = [vrmm_hat x_pred(4:6)];
            vdw_hat = [vdw_hat x_pred(1:3)];
            vrmm_diag_cov = [vrmm_diag_cov rmm_diag_cov(4:6)];
        end
        vrmm_torq = [vrmm_torq cross(mom_res, earth_field_b)];
        vQ = [vQ [Q(4,4); Q(5,5); Q(6,6)]];
        vsingularM_1 = [vsingularM_1; sing_O_1];
        vsingularM_2 = [vsingularM_2; sing_O_2];
        vsingularM_5 = [vsingularM_5; sing_O_5];
        vsingularM_10 = [vsingularM_10; sing_O_10];
        vsingularM_50 = [vsingularM_50; sing_O_50];
        vsun_torq = [vsun_torq sun_torq'];
        vdrag_torq = [vdrag_torq drag_torq'];
        vgrav_grad = [vgrav_grad grad_grav];
        vgamma_avas = [vgamma_avas gamma_avas];
        vgamma_avas(:,1) = vgamma_avas(:,2);
    
    end

  

    writematrix(vdq, "mc/" + timestamp + "/dq.csv", 'WriteMode', 'append')
    writematrix(vdqs, "mc/" + timestamp + "/dqs.csv", 'WriteMode', 'append')
    writematrix(vdw, "mc/" + timestamp + "/dw.csv", 'WriteMode', 'append')
    writematrix(vmag_mom, "mc/" + timestamp + "/mag_mom.csv", 'WriteMode', 'append')
    writematrix(vext_torq, "mc/" + timestamp + "/ext_torq.csv", 'WriteMode', 'append')
    writematrix(vrmm_hat, "mc/" + timestamp + "/drmm_hat.csv", 'WriteMode', 'append')
    clear controller
    clear vdq
    clear vdqs
    clear vdw
    clear vmag_mom
    clear vext_torq
    clear vrmm_hat
    


end
    
    

