% Parameters contains all the constant values needed for the simulation
%     all parameters are defined as global variables
%     in order to use a constant within a function load it as
%     global par_name;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Physical Constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global DEGTORAD
global RADTODEG
global MU_EARTH
global REQ_EARTH
global REQ_EARTH2
global EARTHOMEGA
global EARTH_J2
global EARTH_J3
global EARTH_J4
global TROPICAL_YEAR
global RAAN_dot
global e_EARTH

DEGTORAD      = pi/180;             % deg --> rad
RADTODEG      = 180/pi;             % rad --> deg
MU_EARTH      = 3.986004418e14;     % Earth's gravitational parameter [m^3/s^2]
REQ_EARTH     = 6.378135e06;        % Earth's equatorial radius [m]
REQ_EARTH2    = REQ_EARTH*REQ_EARTH;    % Earth's equatorial radius squared [m^2]
EARTHOMEGA    = 7.2921151467e-5;    % Earth's angular velocity [rad/s]
EARTH_J2      = 1.082616e-3;        % Gravitational spherical harmonic J2 [-]
EARTH_J3      = -2.538810e-6;       % Gravitational spherical harmonic J3 [-]
EARTH_J4      = -1.655970e-6;       % Gravitational spherical harmonic J4 [-]
TROPICAL_YEAR = 365.2421897*24*3600;    % tropical year [s]
RAAN_dot      = 2*pi/TROPICAL_YEAR; % RAAN precession to ensure sun-sync [rad/s]
f_EARTH       = 1/298.257223563;    % Earth's flattening [-]
e_EARTH       = sqrt(1-(1-f_EARTH)^2);  % Earth's eccentricity [-]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orbit Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Orbit keplerian elements
kep_elem_chief  = [REQ_EARTH + 5.26043425e+05,  % semimajor axis [m]
                   0.01,            % eccentricity [-]
                   1.70171324,      % inclination [rad]
                   3.20039430,      % RAAN [rad]
                   pi/2,            % argument of perigee [rad]
                   0];              % mean anomaly [rad]
kep_elem_deputy = kep_elem_chief + [0,
                                    0,
                                    0,
                                    0,
                                    0,
                                    0];

orbit_period = 2*pi*sqrt(kep_elem_chief(1)^3/MU_EARTH); % Orbit period [s]
mean_motion  = 2*pi/orbit_period;   % mean motion [rad/s]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Satellite Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Platform
mass          = 140;                % [kg]
inertia       = [16 0 0; 0 35 0; 0 0 30];   % [kg.m^2]
inv_iner      = inv(inertia);       % [(kg.m^2)^1]
surface_drag  = 0.6*0.8;            % surface area normal to velocity [m^2]
drag_cd       = 2.2;                % drag constant [-]
surface_panel = 0.7*4.9;            % solar panel surface area [m^2]
surface_sun   = surface_drag+surface_panel; % surface area normal to sun [m^2]
reflectivity  = (0.94*surface_drag+0.1*surface_panel)/surface_sun;  % surface light reflectivity [-]
% mag_moment    = [0;0;0];            % satellite magnetic moment [A.m^2]

% Thruster --- MultiFEEP 6U, Morpheus Space
min_thrust  = 1E-6;                 % Minimum control output [N]
max_thrust  = 1E-3;                 % Maximum control output [N]

% Reaction Wheels --- Trillian-1, AAC Clyde Space
RW_torlim = 0.471;                  % torque limit [N.m]
RW_vellim = 6500*60*2*pi;           % velocity limit [rad/s]
RW_iner   = 1.763e-3*eye(4);        % [kg.m^2]
RW_axproj = diag(inertia);          % projection of the RW orientation in each axis
RW_axproj = RW_axproj/norm(RW_axproj);
RW_cR     = RW_axproj.*[1  1  1  1; % RW configuration matrix
                        1  1 -1 -1;
                        1 -1  1 -1];
RW_invcR  = RW_cR'*inv(RW_cR*RW_cR');   % pseudoinverse
RW_nullcR = null(RW_cR);            % null space
RW_Pnull  = RW_nullcR*RW_nullcR'/(RW_nullcR'*RW_nullcR);    % projection matrix on null space of RW_cR
RW_resistance = [450 950 4000 4600 6500 6800;   % RW drag + counter-electromotive force
                 4.3 7.0 16.4 18.5 47.1 53.0].*[60*2*pi; 1e-3]; % speed [rad/s]; torque [N.m]
RW_resistance = [flip(RW_resistance,2), RW_resistance];
RW_res_f  = @(omega) interp1(RW_resistance(1,:), RW_resistance(2,:), omega);

% Magnetorquers --- MTQ800, AAC Clyde Space
MTQ_mom = 15;                       % magnetic moment [A.m^2]

% GNSS
GNSS_pos_error = 0;                 % position estimator error [m]
GNSS_vel_error = 0;                 % velocity estimator error [m/s]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Controller Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Orbit control
ROE_alpha_gain         = 2*0.5*max_thrust/0.00005*1000000/50/5.5; % 0 cont�nuo, 1 pulsado
ROE_alpha_gain_horizon = 0.75;
horizon_weight         = 0.05/2;
east_limit_prop        = (155)*DEGTORAD;
west_limit_prop        = (-140)*DEGTORAD;
minimum_thrust_time    = orbit_period/2*0.85;

% Attitude control



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Environmental Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

atm_rho       = 1.56E-13;           % atmosphere density [kg/m^3] see http://braeunig.us/space/atmos.htm
solar_flux    = 5e-6;               % solar radiation flux [kg/m/s^2]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ephemerides date and time
year = 2021; month = 12; day = 24;
mjd  = djm(day, month, year);       % Modified Julian date
hour = 23; minutes = 59; seconds = 50;
dfra = time_to_dayf(hour, minutes, seconds);    % UTC time [s]
mjdo = djm(1, 1, year);             % modified julian date of 1/1/year
mjd1 = djm(1, 1, year+1);           % modified julian date of 1/1/(year+1)
year_frac = year + (mjd - mjdo)/(mjd1 - mjdo);  % year and fraction

% Propagation time
tstart         = 0;                 % initial time [s]
tstep_orbit    = 60*2;              % time step [s]
tend_orbit     = 24*3600*3;         % end time [s]
tstep_att.bdot = 5;                 % time step [s]
tstep_att.sun  = 5;                 % time step [s]
tstep_att.safe = 5;                 % time step [s]
tstep_att.nom  = 5;                 % time step [s]
tstep_att.sar  = 5;                 % time step [s]
tend_att       = orbit_period*5;    % end time [s]
t_safe         = orbit_period;      % minimum time in safe modes [s]

bounded_control = 1;                % Set to 1 for bounded control, 0 for unbounded control

% ODE solver precision
options = odeset('abstol', 1e-4, 'reltol', 1e-4);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Orbit state vectors
ECI_chief  = kepel_statvec(kep_elem_chief)';    % [m]
ECI_deputy = kepel_statvec(kep_elem_deputy)';   % [m/s]

% Initial attitude and rate
euler_zyx = [30; 50; 20]*pi/180;    % yaw, pitch, roll [rad]
quat      = ezyxquat(euler_zyx);    % quaternion
ang_vel   = [0.25; 0.3; 0.5]'*pi/180;   % angular velocity [rad/s]

% Initial control force
F_control_chief  = [0;0;0];         % [N]
F_control_deputy = [0;0;0];         % [N]

RW_mom = zeros(4,1);                % RW angular momentum