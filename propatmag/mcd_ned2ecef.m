function Ceci = mcd_ned2ecef(lat, lon)
% Matriz cambio de base de coordenadas NED a ECEF.
%
% IN:
%    lat: Latitud Geocentrica, [rad].
%    lon: Longitud Geocentrica, [rad].
% OUT:
%    Ceci: Matriz Coseno Directores.
%
% see: "Sistemas de Navegacion Integrada con Aplicaciones (2da Ed. 2020)"
%       Martin España, CONAE.
%

    Cned2enu = [0 1 0; 1 0 0; 0 0 -1]; % Nort-East-Down to East-North-Up.

    Cenu2ecef = [-sin(lon)            cos(lon)                   0;
                 -sin(lat)*cos(lon)  -sin(lat)*sin(lon)   cos(lat);
                  cos(lat)*cos(lon)   cos(lat)*sin(lat)   sin(lat) ]';

%    Cecef2eci = [ cos(Oe * t) -sin(Oe * t) 0;
%                  sin(Oe * t)  cos(Oe * t) 0;
%                           0            0  1];

    Cecef = Cenu2ecef * Cned2enu;

end