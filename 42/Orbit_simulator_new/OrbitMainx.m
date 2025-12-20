clear all; close all;

Parameters;                         % all constant parameters and initial conditions


% orbital elements
[~, mean_elem_deputy, ~] = testfastmean(ECI_deputy(1:3),ECI_deputy(4:6));
[~, mean_elem_chief , ~] = testfastmean(ECI_chief(1:3),ECI_chief(4:6));


% Initial vectors and matrices
time = tstart;
sun_disturbance = zeros(3,1);
LLH = ecef2llh(inertial_to_terrestrial(gst(mjd, dfra), ECI_deputy'));
ROE_delta_maneuver = zeros(6,1);
Delta_V = zeros(2,1);
Control_RTN = zeros(3,1); control_RTN = zeros(3,1);
relative_position_RTN = zeros(3,1);
delta_alpha = zeros(6,1);
delta_xi = zeros(6,1);
thrust_on_time = 0; thrust_seconds = 0;
condition_final = 0;
time_since_thrust = orbit_period*2;
minimum_wait_thrust = 1;
Thrust_angle = 0; thrust_angle = 0; thrust_angle_ = 0;
Thrust_angle_dot = 0; thrust_angle_dot = 0;
ROE2 = 0;

tpos = 1;
first = 1;

for t = tstart:tstep_orbit:tend_orbit
    % ODE Solver parameters
    tspan = [t, t+tstep_orbit/2, t+tstep_orbit];
    
    % External forces
    F_drag_ECI = 0.5*drag_cd*atm_rho*surface_drag*norm(ECI_deputy(4:6))*ECI_deputy(4:6);
    F_sun_ECI = surface_sun*(1+reflectivity)*solar_flux*-sun_dir(mjd, dfra+t);
    
    % Orbit propagation
    [T,propagation] = ode45(@orbit_prop, tspan, ECI_chief,  options, mass, F_control_chief,               [[0;0;0];1]);
    ECI_chief = propagation(3,:)';
    [T,propagation] = ode45(@orbit_prop, tspan, ECI_deputy, options, mass, F_control_deputy-F_drag_ECI-F_sun_ECI, [[0;0;0];1]);
    ECI_deputy = propagation(3,:)';
    
    % ECI --> RTN
    orbit_radial = ECI_deputy(1:3)/norm(ECI_deputy(1:3));
    orbit_normal = cross(ECI_deputy(1:3),ECI_deputy(4:6))/norm(cross(ECI_deputy(1:3),ECI_deputy(4:6)));
    orbit_trans  = cross(orbit_normal,orbit_radial)/norm(cross(orbit_normal,orbit_radial));
    ECI2RTN      = [orbit_radial, orbit_trans, orbit_normal]';
    
    % Relative positions
    rel_position_ECI = ECI_deputy(1:3)-ECI_chief(1:3);
    rel_position_RTN = ECI2RTN*rel_position_ECI;
    relative_position_RTN = [relative_position_RTN, rel_position_RTN];
    
    % Position and velocity --> relative orbital elements
    [~, mean_elem_deputy, mean_Ustinov_deputy] = testfastmean( ECI_deputy(1:3) + GNSS_pos_error*randn(size(ECI_deputy(1:3))), ECI_deputy(4:6) + GNSS_vel_error*randn(size(ECI_deputy(1:3))));
    [~, mean_elem_chief,  mean_Ustinov_chief ] = testfastmean( ECI_chief(1:3)  + GNSS_pos_error*randn(size(ECI_deputy(1:3))), ECI_chief(4:6)  + GNSS_vel_error*randn(size(ECI_deputy(1:3))));
    ROE_Tc = fTc(mean_elem_chief);
    U_mean_elem_deputy = [mean_elem_deputy(1), mean_Ustinov_deputy(3), mean_Ustinov_deputy(2), mean_Ustinov_deputy(1), mean_elem_deputy(3), mean_elem_deputy(4)]';
    U_mean_elem_chief  = [mean_elem_chief(1),  mean_Ustinov_chief(3),  mean_Ustinov_chief(2),  mean_Ustinov_chief(1),  mean_elem_chief(3),  mean_elem_chief(4) ]';
    U_mean_elem_maneuver = U_mean_elem_chief + ROE_delta_maneuver;
    ROE_delta_xi = U_mean_elem_deputy - U_mean_elem_maneuver;
    ROE_delta_alpha  = ROE_Tc * ROE_delta_xi;
    
    % Allowable control zones - initial position
    ECEF_deputy = inertial_to_terrestrial(gst(mjd, dfra+t), ECI_deputy');
    LLH_deputy = ecef2llh(ECEF_deputy);
    condition_initial = LLH_deputy(2) > east_limit_prop | LLH_deputy(2) < west_limit_prop;
    
    minimum_wait_thrust=1;
    
    if condition_initial & condition_final & minimum_wait_thrust
        
        
        % ROE control
        a_maneuver = ROE_alpha_gain*ROE_delta_alpha(2);
        ROE_delta_maneuver = [a_maneuver -ROE2 0.0 0 0.0 0.0]'*0;
        lambda_deputy = mean_Ustinov_deputy(3);
        ROE_Bd = 1/(mean_elem_chief(1)*sqrt(MU_EARTH/mean_elem_chief(1)^3)) *[[                  0,                    2,                  0];
                                                                              [                 -2,                    0,                  0];
                                                                              [ sin(lambda_deputy), 2*cos(lambda_deputy),                  0];
                                                                              [-cos(lambda_deputy), 2*sin(lambda_deputy),                  0];
                                                                              [                  0,                    0, cos(lambda_deputy)];
                                                                              [                  0,                    0, sin(lambda_deputy)]];
        [ROE_Bd_horizon] = fMw(mean_Ustinov_deputy, mean_elem_chief, horizon_weight, 1);
        control_RTN_horizon = - ROE_alpha_gain_horizon * pinv(ROE_Bd_horizon) * ROE_delta_alpha;
        control_RTN_noperturbations = [0;control_RTN_horizon(1:2,1)];
        
        % J2 compensation
        ROE_A_J2 = fAJ2(mean_Ustinov_chief, mean_elem_chief);
        J2_effect_RTN = pinv(ROE_Bd)*ROE_A_J2*(ROE_delta_alpha+ROE_Tc*ROE_delta_maneuver);
        control_RTN = control_RTN_noperturbations - J2_effect_RTN;
        
        % Drag compensation
        control_RTN(2) = control_RTN(2) + norm(F_drag_ECI)*7*97/40;

        
        if first==1
            tpos
            first=0;
        end;

        if tpos==1 && control_RTN(2)<0
            %control_RTN=0*control_RTN;
        end;
        if tpos==0 && control_RTN(2)>0
            %control_RTN=0*control_RTN;
        end;
        
        % Control saturation
        if norm(control_RTN)>max_thrust && bounded_control==1
            control_RTN = control_RTN/norm(control_RTN)*max_thrust;
        elseif norm(ROE_delta_alpha) < min_thrust
            control_RTN = 0*control_RTN;
        end;
        
        % Attitude requirements
        thrust_angle = atan(control_RTN(3)./control_RTN(2));
        thrust_angle_dot = (thrust_angle-thrust_angle_)/tstep_orbit * (thrust_seconds>0);
        thrust_angle_ = thrust_angle;
        
        % Control timers
        thrust_seconds = thrust_seconds+tstep_orbit;
        time_since_thrust = 0;
    else
        
        if first==0
        
          if tpos==1
              
              tpos=0;
              
          else
              
              tpos=1;
          
          end;

          first=1;
        
        end;
        
        
        % No control
        control_RTN = 0*control_RTN;
        ROE2 = ROE_delta_alpha(2);
        ROE_delta_maneuver = zeros(6,1);
        
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
        time_since_thrust = time_since_thrust+tstep_orbit;
        minimum_wait_thrust = time_since_thrust > 2*orbit_period;
    end;
    
    % Fcontrol en ECI
    F_control_deputy = ECI2RTN'*control_RTN;
    delta_v = abs(control_RTN(2:3))/mass*tstep_orbit;
    
    % Store data to be plotted
    time = [time, t];
    % keorb = [keorb, kep2]);
    sun_disturbance = [sun_disturbance, ECI2RTN*F_sun_ECI];
    LLH = [LLH, LLH_deputy];
    % U_orbit_mean_chief = [U_orbit_mean_chief, U_mean_elem_chief];
    Delta_V = [Delta_V, delta_v];
    Control_RTN = [Control_RTN, control_RTN];
    delta_alpha = [delta_alpha, ROE_delta_alpha];
    delta_xi = [delta_xi, ROE_delta_xi];
    thrust_on_time = [thrust_on_time, thrust_seconds];
    Thrust_angle = [Thrust_angle, thrust_angle];
    Thrust_angle_dot = [Thrust_angle_dot, thrust_angle_dot];
end;


% Output visualization
orbits = time/orbit_period;


orbitradius=norm(ECI_deputy(1:3));
figure(11); hold on; grid on;
plot(orbits, delta_alpha(1,:)*orbitradius,'r');
plot(orbits, delta_alpha(2,:)*orbitradius,'g');
plot(orbits, delta_alpha(3,:)*orbitradius,'b');
plot(orbits, delta_alpha(4,:)*orbitradius,'m');
% plot(orbits, ones(size(delta_alpha(4,:)))*0.00001,'m--');
plot(orbits, delta_alpha(5,:)*orbitradius,'k');
% plot(orbits, -ones(size(delta_alpha(4,:)))*0.00001,'k--');
plot(orbits, delta_alpha(6,:)*orbitradius,'c');
xlabel('Orbits')
ylabel('\delta \alpha_{cm}^d*a [meters]')
title('\delta \alpha_{cm}^d*a')


figure(1); hold on; grid on;
plot(orbits, delta_alpha(1,:),'r');
plot(orbits, delta_alpha(2,:),'g');
plot(orbits, delta_alpha(3,:),'b');
plot(orbits, delta_alpha(4,:),'m');
% plot(orbits, ones(size(delta_alpha(4,:)))*0.00001,'m--');
plot(orbits, delta_alpha(5,:),'k');
% plot(orbits, -ones(size(delta_alpha(4,:)))*0.00001,'k--');
plot(orbits, delta_alpha(6,:),'c');
xlabel('Orbits')
ylabel('\delta \alpha_{cm}^d')
title('\delta \alpha_{cm}^d')

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
plot(orbits, cumsum(Delta_V(1,:)));
plot(orbits, cumsum(Delta_V(2,:)));
plot(orbits, cumsum(Delta_V(2,:))+cumsum(Delta_V(1,:)));
legend('Transverse', 'Normal', 'Total')
title('Fuel consumption')
ylabel('\Delta V [m/s]');
xlabel('Orbits')

figure(4); hold on; grid on;
plot(orbits, Control_RTN(1,:));
plot(orbits, Control_RTN(2,:));
plot(orbits, Control_RTN(3,:));
legend('Radial', 'Transverse', 'Normal')
title('Control Thrust')
ylabel('Thrust [N]');
xlabel('Orbits')

figure(5); hold on; grid on;
plot(relative_position_RTN(2,:), relative_position_RTN(1,:));
plot(relative_position_RTN(3,:), relative_position_RTN(1,:));
title('Relative Position')
legend('Radial vs. Transverse', 'Radial vs. Normal')
ylabel('Radial [m]');
xlabel('Transverse / Normal [m]')

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

figure(7); hold on; grid on;
plot(LLH(1,2:end)*RADTODEG,sun_disturbance(1,2:end))
title('Solar pressure force on the Radial component')
ylabel('Force [N]');
xlabel('Latitude [deg]')

figure(8); hold on; grid on;
plot(orbits, thrust_on_time/60)
title('Maneuver Length')
xlabel('Orbits');
ylabel('Time [min]')

figure(9); hold on; grid on;
plot(orbits, LLH(2,:)*RADTODEG)
title('Longitude')
xlabel('Orbits');
ylabel('Longitude [^\circ]')