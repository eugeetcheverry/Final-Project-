clear all; %close all;

Parameters;                         % all constant parameters and initial conditions

% orbital elements
[~, mean_elem_deputy, ~] = testfastmean(ECI_deputy(1:3),ECI_deputy(4:6));
[~, mean_elem_chief , ~] = testfastmean(ECI_chief(1:3),ECI_chief(4:6));

% 1 smi-continuous
% 2 continuous
% 3 no contorl
% 4 within longitude
% 5 over poles

% Initial vectors and matrices
time = tstart;
sun_disturbance = zeros(3,1);
LLH = ecef2llh(inertial_to_terrestrial(gst(mjd, dfra), ECI_deputy'));
ROE_delta_maneuver = zeros(6,1);
Delta_V = zeros(2,1);
Control_RTN = zeros(3,1); control_RTN = zeros(3,1);
relative_position_RTN = zeros(3,1);
delta_alpha = zeros(6,1);
delta_alphaf = zeros(6,1);
delta_alpha_orig = zeros(6,1);
delta_alpha_dot = zeros(6,1);
delta_xi = zeros(6,1);
thrust_on_time = 0; thrust_seconds = 0;
condition_final = 0;
time_since_thrust = orbit_period*2;
thrust_wait = 1;
Thrust_angle = 0; thrust_angle = 0; thrust_angle_ = 0;
Thrust_angle_dot = 0; thrust_angle_dot = 0;
ROE2 = 0;
Mean_elem_chief = zeros(6,1); mean_elem_chief = zeros(6,1);
Mean_elem_deputy = zeros(6,1); mean_elem_deputy = zeros(6,1);
ECI_chief_sun = ECI_chief;
N_captures = 0; total_captures = 0;
In_capture = 0; in_capture = 0;

ROE_delta_alpha_filtered = zeros(6,1);

if control_strategy == 1
    control_state_machine = 3;
    thrust_seconds_last = minimum_thrust_time*1.5;
elseif control_strategy == 2
    control_state_machine = 2;
    thrust_seconds_last = 1;
    time_since_thrust_last = thrust_seconds_last;
    thruster_status = 1;
end


for t = tstart:orbit_control_tau:tend_orbit
    % ODE Solver parameters
    tspan = [t, t+orbit_control_tau/2, t+orbit_control_tau];
    
    % External forces
    F_drag_ECI = 0.5*drag_cd*atm_rho*surface_drag*norm(ECI_deputy(4:6))*ECI_deputy(4:6);
    F_sun_ECI = surface_sun*(1+reflectivity)*solar_flux*-sun_dir(mjd, dfra+t);
    
    % Orbit propagation
    [T,propagation] = ode45(@orbit_prop, tspan, ECI_chief,  options, mass, F_control_chief, [[0;0;0];1]);
    ECI_chief = propagation(3,:)';
    % [T,propagation] = ode45(@orbit_prop, tspan, ECI_chief_sun,  options, mass, F_control_chief-F_sun_ECI, [[0;0;0];1]);
    % ECI_chief_sun = propagation(3,:)';
    ECI_chief_sun = ECI_chief;
    % ECI_chief = ECI_chief_sun;
    [T,propagation] = ode45(@orbit_prop, tspan, ECI_deputy, options, mass, F_control_deputy-F_drag_ECI-F_sun_ECI, [[0;0;0];1]);
    ECI_deputy = propagation(3,:)';
    
    % ECI --> RTN
    orbit_radial = ECI_deputy(1:3)/norm(ECI_deputy(1:3));
    orbit_normal = cross(ECI_deputy(1:3),ECI_deputy(4:6))/norm(cross(ECI_deputy(1:3),ECI_deputy(4:6)));
    orbit_trans  = cross(orbit_normal,orbit_radial)/norm(cross(orbit_normal,orbit_radial));
    ECI2RTN      = [orbit_radial, orbit_trans, orbit_normal]';
    
    % Position and velocity --> relative orbital elements
    [~, mean_elem_deputy, mean_Ustinov_deputy] = testfastmean(ECI_deputy(1:3) + GNSS_pos_error*randn(size(ECI_deputy(1:3))), ECI_deputy(4:6) + GNSS_vel_error*randn(size(ECI_deputy(1:3))));
    
    % [~, mean_elem_chief,  mean_Ustinov_chief ] = testfastmean(ECI_chief(1:3), ECI_chief(4:6));
    % ROE_Tc = fTc(mean_elem_chief);
    % U_mean_elem_deputy = [mean_elem_deputy(1), mean_Ustinov_deputy(3), mean_Ustinov_deputy(2), mean_Ustinov_deputy(1), mean_elem_deputy(3), mean_elem_deputy(4)]';
    % U_mean_elem_chief  = [mean_elem_chief(1),  mean_Ustinov_chief(3),  mean_Ustinov_chief(2),  mean_Ustinov_chief(1),  mean_elem_chief(3),  mean_elem_chief(4) ]';
    % U_mean_elem_maneuver = U_mean_elem_chief + ROE_delta_maneuver;
    % ROE_delta_xi = U_mean_elem_deputy - U_mean_elem_maneuver;
    % ROE_delta_alpha_orig  = ROE_Tc * ROE_delta_xi;
    
    [~, mean_elem_chief,  mean_Ustinov_chief ] = testfastmean(ECI_chief_sun(1:3), ECI_chief_sun(4:6));
    ROE_Tc = fTc(mean_elem_chief);
    U_mean_elem_deputy = [mean_elem_deputy(1), mean_Ustinov_deputy(3), mean_Ustinov_deputy(2), mean_Ustinov_deputy(1), mean_elem_deputy(3), mean_elem_deputy(4)]';
    U_mean_elem_chief  = [mean_elem_chief(1),  mean_Ustinov_chief(3),  mean_Ustinov_chief(2),  mean_Ustinov_chief(1),  mean_elem_chief(3),  mean_elem_chief(4) ]';
    U_mean_elem_maneuver = U_mean_elem_chief + ROE_delta_maneuver;
    ROE_delta_xi = U_mean_elem_deputy - U_mean_elem_maneuver;
    ROE_delta_alpha  = ROE_Tc * ROE_delta_xi;
    lambda_deputy = mean_Ustinov_deputy(3);
    ROE_Bd = 1/(mean_elem_chief(1)*sqrt(MU_EARTH/mean_elem_chief(1)^3)) *[[                  0,                    2,                  0];
                                                                          [                 -2,                    0,                  0];
                                                                          [ sin(lambda_deputy), 2*cos(lambda_deputy),                  0];
                                                                          [-cos(lambda_deputy), 2*sin(lambda_deputy),                  0];
                                                                          [                  0,                    0, cos(lambda_deputy)];
                                                                          [                  0,                    0, sin(lambda_deputy)]];
    [ROE_Bd_horizon] = fMw(mean_Ustinov_deputy, mean_elem_chief, horizon_weight, 1-controlR);
    ROE_A_J2 = fAJ2(mean_Ustinov_chief, mean_elem_chief);
    
    % ROE averaging filter
    ROE_delta_alpha_filtered = forg*ROE_delta_alpha + (1-forg)*ROE_delta_alpha_filtered;
    
    %if ROE_filtering & (t/orbit_control_tau > ROE_filter_N)
    %    ROE_delta_alpha_filtered = sum(delta_alpha(:,end-ROE_filter_N+1:end),2)/ROE_filter_N;
    %    ROE_delta_alpha_filtered = ROE_delta_alpha_filtered + sum(delta_alpha_dot(:,end-ROE_filter_N+1:end).*ROE_filter_scale,2)/(ROE_filter_N^2/2);
    %end;
    
    if ROE_filtering & (t/orbit_control_tau > ROE_filter_N)
    
        if first == 1
            
            suma_alfa = sum(delta_alpha(:,end-ROE_filter_N+1:end),2)/ROE_filter_N;
            suma_dalfa = sum(ddal(:,(round(t/tstep)-Lav):round(t/tstep)).*wL,2)/(Lav*(Lav-1)/2);

            ROE_delta_alpha_filtered = sum(delta_alpha(:,end-ROE_filter_N+1:end),2)/ROE_filter_N;
            ROE_delta_alpha_filtered = ROE_delta_alpha_filtered + sum(delta_alpha_dot(:,end-ROE_filter_N+1:end).*ROE_filter_scale,2)/(ROE_filter_N^2/2);

            first = 0;
          
        else
            
            suma_alfa = suma_alfa + (dal(:,round(t/tstep)) - dal(:,round(t/tstep)-Lav-1))/Lav;
            suma_dalfa = suma_dalfa + (Lav * ddal(:,end) - sum(ddal(:,(round(t/tstep)-Lav-1):(round(t/tstep)-1))))/(Lav*(Lav-1)/2);
                        
            ROE_delta_alpha_filtered = suma_alfa + suma_dalfa; %sum(ddal(:,(round(t/tstep)-Lav):round(t/tstep)).*wL,2)/(Lav^2/2);
            
        end;
        

    end;
    
    if t/orbit_control_tau < ROE_filter_N
        control_state_machine == 3;
    end
    
    % Allowable control zones - initial position
    ECEF_deputy = inertial_to_terrestrial(gst(mjd, dfra+t), ECI_deputy');
    LLH_deputy = ecef2llh(ECEF_deputy);
    longitude_condition = LLH_deputy(2) > east_limit_prop | LLH_deputy(2) < west_limit_prop;
    
    if control_state_machine == 3
        thruster_status = 0;
        
        % Allowable control zones - end position
        [T,propagation] = ode45(@orbit_prop, [0,minimum_thrust_time], ECI_chief, options, mass, F_control_chief, [[0;0;0];1]);
        ECI_chief_future = propagation(end,:)';
        ECEF_chief_future = inertial_to_terrestrial(gst(mjd, dfra+t), ECI_chief_future');
        LLH_chief_future = ecef2llh(ECEF_chief_future);
        condition_final = LLH_chief_future(2) > east_limit_prop | LLH_chief_future(2) < west_limit_prop;
        
        % Attitude requirements
        thrust_angle = 0;
        thrust_angle_dot = 0;
        
        % Control timers
        thrust_seconds = 0;
        time_since_thrust = time_since_thrust+orbit_control_tau;
        thrust_wait = time_since_thrust > minimum_wait_time;
        
        % State switch
        if longitude_condition & condition_final & thrust_wait
            control_state_machine = 4;
            time_since_thrust_last = time_since_thrust;
            % ECI_chief_sun = ECI_chief;
        end
    elseif control_state_machine == 4
        thruster_status = 1;
        
        % Control timers
        time_since_thrust = 0;
        
        % State switch
        if ~longitude_condition
            control_state_machine = 5;
            latitude_last = LLH_deputy(1);
        end
    elseif control_state_machine == 5
        thruster_status = 1;
        
        % Allowable control zones - poles
        from_pole = (abs(LLH_deputy(1))-abs(latitude_last)) < 0;
        if LLH_deputy(1) > 0
            pole_end = LLH_deputy(1) < latitude_limit_north;
        else
            pole_end = LLH_deputy(1) > latitude_limit_south;
        end
        latitude_last = LLH_deputy(1);
        
        % Control timers
        time_since_thrust = 0;
        
        % State switch
        if from_pole & pole_end
            control_state_machine = 3;
            thrust_seconds_last = thrust_seconds+orbit_control_tau;
        end
    elseif control_state_machine == 2
        if (~longitude_condition) & (LLH_deputy(1) < latitude_limit_north) & (LLH_deputy(1) > latitude_limit_south)
            if N_captures < 3
                in_capture = rand < capture_probability;
                if in_capture
                    N_captures = N_captures + 1;
                    total_captures = total_captures + 1;
                end
            else
                in_capture = 0;
            end
        else
            N_captures = 0;
        end
    else
        disp('State machine input error!')
        break
    end
    
    if thruster_status
        % Control timers
        thrust_seconds = thrust_seconds+orbit_control_tau;
        
        % ROE control
        a_maneuver = ROE_alpha_gain*ROE_delta_alpha_filtered(2);
        if (t > tstart) & (control_strategy == 2)
            a_maneuver = a_maneuver + ROE_alphadot_gain*(ROE_delta_alpha_filtered(2)-ROE_delta_alpha2_old)/orbit_control_tau;
        end
        ROE_delta_alpha2_old = ROE_delta_alpha_filtered(2);
        ROE_delta_maneuver = [a_maneuver -ROE2*0 0.0 0 0.0 0.0]';
        control_RTN_horizon = - ROE_alpha_gain_horizon * pinv(ROE_Bd_horizon) * ROE_delta_alpha_filtered;
               
        if controlR == 0
            control_RTN_noperturbations = [0;control_RTN_horizon(1:2,1)];
        else
            control_RTN_noperturbations = control_RTN_horizon(1:3,1);
        end;
        
        % J2 compensation
        J2_effect_RTN = pinv(ROE_Bd)*ROE_A_J2*(ROE_delta_alpha_filtered+ROE_Tc*ROE_delta_maneuver);
        control_RTN = control_RTN_noperturbations - J2_effect_RTN;
        
        % Drag compensation
        control_RTN(2) = control_RTN(2) + norm(F_drag_ECI)*time_since_thrust_last/thrust_seconds_last;
        
        % Control saturation
        if norm(control_RTN)>max_thrust && bounded_control==1
            control_RTN = control_RTN/norm(control_RTN)*max_thrust;
        elseif norm(control_RTN) < min_thrust
            control_RTN = control_RTN/norm(control_RTN) * min_thrust*(norm(control_RTN) < min_thrust/2);
        end

        % In capture transverse-only control
        if in_capture
            control_RTN(1) = 0;
            control_RTN(3) = 0;
        end
        
        % Saturación cono normal
        if abs(control_RTN(3))>abs(control_RTN(2))*sin(angleNT)
            control_RTN(3) = control_RTN(3)/abs(control_RTN(3))*sin(angleNT)*abs(control_RTN(2));
        end;

        % Saturación cono radial
        if controlR==1 && abs(control_RTN(1))>abs(control_RTN(2))*sin(angleNT)
            control_RTN(1) = control_RTN(1)/abs(control_RTN(1))*sin(angleNT)*abs(control_RTN(2));
        end;

        % Cancelación inciial
        if t/orbit_control_tau < ROE_filter_N
            control_RTN = 0*control_RTN;
        end;
        
        if zerocontrol==1
            control_RTN = 0*control_RTN;
        end;
        
        % Attitude requirements
        thrust_angle = atan(control_RTN(3)./control_RTN(2));
        thrust_angle_dot = (thrust_angle-thrust_angle_)/orbit_control_tau * (thrust_seconds>0);
        thrust_angle_ = thrust_angle;
    else
        % No control
        control_RTN = 0*control_RTN;
        ROE2 = ROE_delta_alpha_filtered(2);
        ROE_delta_maneuver = zeros(6,1);
    end

    % Fcontrol en ECI
    F_control_deputy = ECI2RTN'*control_RTN;
    delta_v = abs(control_RTN(2:3))/mass*tstep_orbit;
    
    % Store data to be plotted
    if ~rem(t, tstep_orbit)
        time = [time, t];
        % keorb = [keorb, kep2]);
        sun_disturbance = [sun_disturbance, ECI2RTN*F_sun_ECI];
        LLH = [LLH, LLH_deputy];
        % U_orbit_mean_chief = [U_orbit_mean_chief, U_mean_elem_chief];
        Delta_V = [Delta_V, delta_v];
        Control_RTN = [Control_RTN, control_RTN];
        delta_alpha = [delta_alpha, ROE_delta_alpha];
        delta_alphaf = [delta_alphaf, ROE_delta_alpha_filtered];
        delta_alpha_dot = [delta_alpha_dot, (ROE_A_J2*(ROE_delta_alpha_filtered+ROE_Tc*ROE_delta_maneuver)-ROE_Bd(:,2)*norm(F_drag_ECI))*orbit_control_tau + ROE_Bd*(control_RTN-ECI2RTN*F_sun_ECI)*orbit_control_tau];
        % delta_alpha_orig = [delta_alpha_orig, ROE_delta_alpha_orig];
        delta_xi = [delta_xi, ROE_delta_xi];
        thrust_on_time = [thrust_on_time, thrust_seconds];
        Thrust_angle = [Thrust_angle, thrust_angle];
        Thrust_angle_dot = [Thrust_angle_dot, thrust_angle_dot];
        Mean_elem_chief = [Mean_elem_chief, mean_elem_chief];
        Mean_elem_deputy = [Mean_elem_deputy, mean_elem_deputy];
        % Relative positions
        rel_position_ECI = ECI_deputy(1:3)-ECI_chief(1:3);
        rel_position_RTN = ECI2RTN*rel_position_ECI;
        relative_position_RTN = [relative_position_RTN, rel_position_RTN];
        In_capture = [In_capture, in_capture];
    end
end


% Output visualization
orbits = time/orbit_period;

%delta_alphaf = delta_alphaf*kep_elem_chief(1);
% delta_alpha_orig = delta_alpha_orig*kep_elem_chief(1);

figure(1); hold on; grid on;
plot(orbits, delta_alphaf(1,:)*kep_elem_chief(1),'r');
plot(orbits, delta_alphaf(2,:)*kep_elem_chief(1),'g');
plot(orbits, delta_alphaf(3,:)*kep_elem_chief(1),'b');
plot(orbits, delta_alphaf(4,:)*kep_elem_chief(1),'m');
% plot(orbits, ones(size(delta_alpha(4,:)))*0.00001,'m--');
plot(orbits, delta_alphaf(5,:)*kep_elem_chief(1),'k');
% plot(orbits, -ones(size(delta_alpha(4,:)))*0.00001,'k--');
plot(orbits, delta_alphaf(6,:)*kep_elem_chief(1),'c');
xlabel('Orbits')
%ylabel('\delta \alpha_{cm}^d * a_c')
%title('\delta \alpha_{cm}^d * a_c')
ylabel('\delta \alpha_{cm}^d * a_c')
title('\delta \alpha_{cm}^d * a_c')

% figure(2); hold on; grid on;
% plot(orbits, delta_alpha_orig(1,:),'r');
% plot(orbits, delta_alpha_orig(2,:),'g');
% plot(orbits, delta_alpha_orig(3,:),'b');
% plot(orbits, delta_alpha_orig(4,:),'m');
% % plot(orbits, ones(size(delta_alpha_orig(4,:)))*0.00001,'m--');
% plot(orbits, delta_alpha_orig(5,:),'k');
% % plot(orbits, -ones(size(delta_alpha_orig(4,:)))*0.00001,'k--');
% plot(orbits, delta_alpha_orig(6,:),'c');
% xlabel('Orbits')
% ylabel('\delta \alpha_{cm}^d * a_c')
% title('\delta \alpha_{cm}^d * a_c')

% figure(2); hold on; grid on;
% plot(orbits, delta_xi(1,:)./U_orbit_mean_chief(1,2:end),'r');
% plot(orbits, delta_xi(2,:),'g');
% plot(orbits, delta_xi(3,:),'b');
% plot(orbits, delta_xi(4,:),'m');
% plot(orbits, delta_xi(5,:),'k');
% plot(orbits, delta_xi(6,:),'c');
% plot(orbits, ones(size(delta_alpha(4,:)))*0.00001,'m--');
% plot(orbits, -ones(size(delta_alpha(4,:)))*0.00001,'k--');
% xlabel('Orbits')
% title('\delta \xi_c^m  (con \xi(1) dividido por a_c)')
% axis([tstart/orbit_period tend_orbit/orbit_period -0.00015 0.00015]);

figure(3); hold on; grid on;
plot(orbits, cumsum(Delta_V(1,:)),'r');
plot(orbits, cumsum(Delta_V(2,:))),'g';
plot(orbits, cumsum(Delta_V(2,:))+cumsum(Delta_V(1,:)),'b');
legend('Transverse', 'Normal', 'Total', 'Location','northwest')
title('Fuel consumption')
ylabel('\Delta V [m/s]');
xlabel('Orbits')

figure(4); hold on; grid on;
plot(orbits, Control_RTN(1,:),'r');
plot(orbits, Control_RTN(2,:),'g');
plot(orbits, Control_RTN(3,:),'b');
legend('Radial', 'Transverse', 'Normal')
title('Control Thrust')
ylabel('Thrust [N]');
xlabel('Orbits')

figure(5)
subplot(2,1,1); hold on; grid on;
plot(relative_position_RTN(2,:), relative_position_RTN(1,:));
plot(relative_position_RTN(3,:), relative_position_RTN(1,:));
title('Relative Position')
legend('Radial vs. Transverse', 'Radial vs. Normal')
ylabel('Radial [m]');
xlabel('Transverse / Normal [m]')
subplot(2,1,2); hold on; grid on;
plot(orbits,relative_position_RTN(1,:),'r');hold on;
plot(orbits,relative_position_RTN(2,:),'g');hold on;
plot(orbits,relative_position_RTN(3,:),'b');hold on;
legend('Radial', 'Transverse', 'Normal')
title('Relative Position')

figure(55);hold on;
plot(orbits,relative_position_RTN(1,:),'r');hold on;
plot(orbits,relative_position_RTN(2,:),'g');hold on;
plot(orbits,relative_position_RTN(3,:),'b');hold on;
legend('Radial', 'Transverse', 'Normal')
title('Relative Position')
grid on;


figure(51);hold on;
plot(orbits,sqrt(relative_position_RTN(1,:).^2+relative_position_RTN(2,:).^2+relative_position_RTN(3,:).^2),'g');hold on;
plot(orbits,sqrt(relative_position_RTN(1,:).^2+relative_position_RTN(3,:).^2),'m');hold on;
title('Orbit Control Error [m]');
xlabel('Orbits');
grid on;

figure(511);hold on;
plot(orbits,sqrt( relative_position_RTN(1,:).^2 + ( relative_position_RTN(3,:)+relative_position_RTN(2,:)./7600*460.*cos(LLH(1,:))).^2 ),'k');hold on;
title('Effective Control Error [m]');
xlabel('Orbits');
grid on;


figure(512);hold on;
plot(orbits,sqrt( relative_position_RTN(1,:).^2*0 + ( relative_position_RTN(3,:)+relative_position_RTN(2,:)./7600*460.*cos(LLH(1,:))).^2 ),'k');hold on;
title('Effective Track Error [m]');
xlabel('Orbits');
grid on;

figure(6);
subplot(2,1,1); hold on; grid on;
plot(orbits, RADTODEG*Thrust_angle);
title('Thrust T->N Angle');
ylabel('\theta [\circ]')
xlabel('Orbits');
subplot(2,1,2); hold on; grid on;
plot(orbits, RADTODEG*Thrust_angle_dot);
title('Thrust T->N Angular Velocity');
ylabel('\omega [\circ/s]')
xlabel('Orbits');

% figure(7); hold on; grid on;
% plot(LLH(1,2:end)*RADTODEG,sun_disturbance(:,2:end))
% title('Solar pressure force')
% ylabel('Force [N]');
% xlabel('Latitude [deg]')
% legend('Radial', 'Transverse', 'Normal')

% figure(8); hold on; grid on;
% plot(orbits, thrust_on_time/60)
% title('Maneuver Length')
% xlabel('Orbits');
% ylabel('Time [min]')

% Mean_elem_chief(6,:) = Mean_elem_chief(6,:) - 2*pi*(Mean_elem_chief(6,:)>pi);
% figure(9)
% subplot(3,1,1); hold on; grid on;
% plot(orbits, LLH(2,:)*RADTODEG.*(thrust_on_time>0))
% title('Propulsion Longitude')
% xlabel('Orbits');
% ylabel('Longitude [^\circ]')
% ylim([-180,180])
% subplot(3,1,2); hold on; grid on;
% plot(orbits, LLH(1,:)*RADTODEG.*(thrust_on_time>0))
% title('Propulsion Latitude')
% xlabel('Orbits');
% ylabel('Latitude [^\circ]')
% ylim([-90,90])
% subplot(3,1,3); hold on; grid on;
% plot(orbits, Mean_elem_chief(6,:)*RADTODEG.*(thrust_on_time>0))
% title('Propulsion Mean True Anomaly Chief')
% xlabel('Orbits');
% ylabel('$\overline{\nu}~ [^\circ]$','Interpreter','latex')
% ylim([-180,180])


% figure(10)
% labels = ['semi-major axis a [m]';
%           'eccentricity e [-]';
%           'inclination i [rad]';
%           'RAAN \Omega [rad]';
%           'argument of perigee \omega [rad]';
%           'true anomaly \nu [rad]'];
% for i=1:6
%     subplot(3,2,i); hold on; grid on;
%     plot(orbits(2:end), Mean_elem_chief(i,2:end))
%     plot(orbits(2:end), Mean_elem_deputy(i,2:end))
%     xlabel('Orbits');
%     ylabel(labels(i))
%     legend('Chief', 'Deputy')
% end

% figure(11)
% labels = ['delta semi-major axis a [m]';
%           'delta eccentricity e [-]';
%           'delta inclination i [rad]';
%           'delta RAAN \Omega [rad]';
%           'delta argument of perigee \omega [rad]';
%           'delta true anomaly \nu [rad]'];
% for i=1:6
%     subplot(3,2,i); hold on; grid on;
%     plot(orbits(2:end), Mean_elem_chief(i,2:end)-Mean_elem_deputy(i,2:end))
%     xlabel('Orbits');
%     ylabel(labels(i))
% end

disp('Total delta v')
disp(sum(Delta_V(2,:))+sum(Delta_V(1,:)))
%disp('Root mean squared error')
%disp(rms(vecnorm(relative_position_RTN)))
disp('Total captures')
disp(total_captures)

% figure(12)
% subplot(2,1,1); hold on; grid on;
% plot(orbits, LLH(2,:)*RADTODEG.*In_capture)
% title('Capture Longitude')
% xlabel('Orbits');
% ylabel('Longitude [^\circ]')
% ylim([-180,180])
% subplot(2,1,2); hold on; grid on;
% plot(orbits, LLH(1,:)*RADTODEG.*In_capture)
% title('Capture Latitude')
% xlabel('Orbits');
% ylabel('Latitude [^\circ]')
% ylim([-90,90])