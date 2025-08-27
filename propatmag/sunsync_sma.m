function smaxis = sunsync_sma (exc, inc, q)
%
%   function smaxis = sunsync_recf (sma, exc, inc, q)
%   This function obtains the orbital semimajor axis of a low Earth orbit
%   given the orbital recovering factor q.
%  
%   Inputs:
%       exc
%           Orbital excentricity
%       inc
%           Orbital inclination (rad)
%       q 
%           Sun-synchronous orbit recovering factor. The recovering factor 
%           is calculated by q = N + M / I, where N, M and I are integers. 
%           N is the number of orbits per day, M is the exceeding orbits 
%           after I days, and I is the recovering period in days. The 
%           recovering factor normally lies between 13 to 16.
%
%   Output:
%       smaxis
%           Recurrent orbital semimajor axis, in meters.
%
%   Author:
%       Valdemir Carrara        Feb. 1992      V 1.0
%       Valder Matos de Medeiros    Mai. 1997
%       Valdemir Carrara        May/2017        Matlab version
% 

    earth_gravity   = 3.9860064e+14;    % Earth's gravity (m**3/s**2)
    earth_rate      = 7.2921158546819492e-5; % Earth's sidereal rotation speed (rad/s) = (1+1/TROPIC_YEAR)*PI_T2/86400

    el      = [6878000, exc, inc, 0, 0, 0];
    ant     = el(1);
    
    epx     = -1.5*sqrt(earth_gravity/(ant^5));
    del     = 1000000;
    ic      = 0;

    while (abs(del/ant) > 1e-9 && ic < 20)
        delk    = delkep (el);
        fact    = delk(5) + delk(6) - q*(earth_rate - delk(4));
        del     = fact/epx;
        sma     = ant - del;
        ant     = sma;
        el(1)   = sma;
        ic      = ic + 1;
    end
    smaxis  = sma;
    
    if (ic > 30)
        disp(' Error in routine sunsync_recf. Interaction did not converge');
    end
    return
    
