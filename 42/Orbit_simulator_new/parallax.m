% Parallax
MU_EARTH    = 3.986004418e14;     % Earth's gravitational parameter [m^3/s^2]
REQ_EARTH   = 6.378135e06;        % Earth's equatorial radius [m]

errorpos    = 0.1;
errorpix    = 0.1*20/1024*pi/180;
bline       = 5000;
bline1      = 2500;

% Definitions
bline2      = bline-bline1;
Hsat        = 600000;
Hrso        = 700000;
p1          = [(REQ_EARTH+Hsat); -bline1;     0];
p2          = [(REQ_EARTH+Hsat);  bline2;     0];
prso        = [(REQ_EARTH+Hrso);       0;     0];
u1          = (prso-p1)/norm(prso-p1);
u2          = (prso-p2)/norm(prso-p2);
d           = p1-p2;

% Ideal solution
M           = [[1 u1'*u2];[u1'*u2 1]];
v           = [d'*u1;d'*u2];
alfa        = inv(M)*v;
alfasimple  = 1/(1-u1'*u2*u1'*u2) * [d'*(eye(3)-u2*u2')*u1; d'*(eye(3)-u1*u1')*u2];
error1      = prso-p1+alfa(1)*u1;
error2      = prso-p2-alfa(2)*u2;
error11      = prso-p1+alfasimple(1)*u1;
error22      = prso-p2-alfasimple(2)*u2;


alfa0       = -d'*(u1-u2)/((u1-u2)'*(u1-u2));
error10     = prso-p1-alfa0*u1;
error20     = prso-p2-alfa0*u2;

[error1 error2]

% Noisy Measurements
p1m         = p1 + errorpos*randn(3,1);
p2m         = p2 + errorpos*randn(3,1);
dm          = p1m-p2m;
dp          = randn(3,1)*errorpix;
S1          = [[1 -dp(3) dp(2)];[dp(3) 1 -dp(1)];[-dp(2) dp(1) 1]];
u1m         = S1*u1;u1m=u1m/norm(u1m);
dp          = randn(3,1)*errorpix;
S2          = [[1 -dp(3) dp(2)];[dp(3) 1 -dp(1)];[-dp(2) dp(1) 1]];
u2m         = S2*u2;u2m=u2m/norm(u2m);
Mm          = [[1 u1m'*u2m];[u1m'*u2m 1]];
vm          = [dm'*u1m;dm'*u2m];
alfam       = inv(Mm)*vm;
alfasimplem = 1/(1-u1m'*u2m*u1m'*u2m) * [dm'*(eye(3)+u2m*u2m')*u1m; dm'*(eye(3)+u1m*u1m')*u2m];
error1m     = prso-p1m+alfam(1)*u1m;
error2m     = prso-p2m-alfam(2)*u2m;

[error1m error2m]
