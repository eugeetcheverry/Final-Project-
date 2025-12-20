function sunsync_inclination = sunsync_inc (sma, exc)
%
%   sunsync_inclination = sunsync_inc (sma, exc)
%   This routine calculates the sunsynchronous inclination or a given 
%   orbit
%
%   Input:
%       sma  
%           Orbit semi-major axis (m)
%       exc
%           Orbit eccentricity
%
%   Output:
%     sunsync_inclination
%           Sun-Synchronous inclination of the orbit (rad).
%
%   Author:
%       Valder Matos de Medeiros        15/05/1987  Fortran version
%       Valdemir Carrara                May/2017    Matlab version
%

    earth_gravity   = 3.9860064e+14;    % Earth's gravity (m**3/s**2)
    tropic_year     = 365.24219879;		% Tropical year (days)
    earth_radius	= 6378139.;         % Earth's radius in meters
    arg     = 1.72;     % first approximated inclination
    j_2     = 1.0826268362e-3;		% J2 = 484.16544e-6 * SQRT(5.e0)
    el      = [sma, exc, arg, 0, 0, 0]; % keplerian elements

    omegap  = 2*pi/tropic_year/86400;  % rate of the ascending node
    amm     = sqrt (earth_gravity/(sma*sma*sma)); % mean mean motion
    con     = -1.5*j_2*amm*earth_radius*earth_radius/(sma*sma);
    del     = 1;
    ic      = 0;

    while (abs(del) > 1e-6 && ic < 20)
        delk    = delkep(el);
        chu     = cos(arg);
        del     = (omegap - delk(4))/con;
        chu     = chu + del;
        arg     = acos(chu);
        el(3)   = arg;
        ic      = ic + 1;
    end
    sunsync_inclination     = arg;

    if (ic > 20)
        disp(' Error in function sunsync_inclination:. Interaction did not converge');
    end
    return
    

