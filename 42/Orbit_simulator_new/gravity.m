function F_g = gravity(R,M)
% Gravity calculates gravity force on satellite
%   F_g = gravity(R,M) takes in satellite ECI position in [m] (column vectors) and
%   mass in [kg] and generates ECI force in [N]

global MU_EARTH
global REQ_EARTH
global REQ_EARTH2
global EARTH_J2
global EARTH_J3
global EARTH_J4

% Base inercial standard: X=equinoccio vernal Z=eje terrestre
Rmag = norm(R); Rmag2 = Rmag^2; Rmag3 = Rmag^3;
K3 = 0.5*MU_EARTH*REQ_EARTH2;
K0 = -MU_EARTH/Rmag3;
K1 = K3/(Rmag3*Rmag2);
zR = R(3)/Rmag; % versor de direccion norte en la terna centrada

% Armonic terms computation
zR2 = zR^2; zR3 = zR^3; zR4 = zR^4;
RE_Rmag = REQ_EARTH/Rmag; RE2_Rmag2 = RE_Rmag^2;
fxy2 = 3*EARTH_J2*(5*zR2-1);
fxy3 = 5*EARTH_J3*RE_Rmag*(7*zR3-3.*zR);
fxy4 = 3.75*EARTH_J4*RE2_Rmag2*(21*zR4-14*zR2+1);
fz2 = 3*EARTH_J2*(5*zR2-3);
fz3 = EARTH_J3*RE_Rmag*(R(3)*(35*zR3-30*zR)+3*Rmag);
fz4 = 1.25*EARTH_J4*RE2_Rmag2*(63*zR4-70*zR2+15);

% Simulador previo
F_g = [(K0*R(1) + K1*R(1)*(fxy2+fxy3+fxy4))*M;
       (K0*R(2) + K1*R(2)*(fxy2+fxy3+fxy4))*M;
       (K0*R(3) + K1*(R(3)*(fz2+fz4)+fz3) )*M];
