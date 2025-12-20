function [oscul_elem, mean_elem, ustin_elem] = testfastmean(r,v)

%%%%%%%%%%%%
% Test of Fast Mean Elements from Cartesian Osculating (r,v):

%%%%%%%%%%%
% Inputs: Position and Velocity in Inertial Frame from Envisat Simulation

% % Sample 1
% r 	= [-6470510.19080675; -818494.420897379; -2974221.1731183];
% v 	= [-3175.7522986603;      826.27582022761; 6694.51634240686];
% % Sample 2
% r 	= [-7153725.3192728; -405436.720846004; 28922.6714521394];
% v 	= [-24.4162398636635;  1109.34746630818; 7377.27705885103];
% % Sample 3
% r 	= [2487103.44168157;-840935.608159341;-6672479.85354085];
% v 	= [  -6966.68295884504; -839.147779774002; -2492.35383938486];


%%%%%%%%%%%%
% Constants

mu 	    = 3.986004418e14;
RE      = 6.378135e06;

%%%%%%%%%%%%
% Compute Osculating Elements

h       = cross(r,v);
a       = norm(r)/(2-norm(r)*norm(v)^2/mu);
e       = sqrt(1-norm(h)^2/(a*mu));
i       = acos(h(3)/norm(h));
n       = cross([0;0;1],h);
evec    = cross(v,h)/mu - r/norm(r);
Om      = acos(n(1)/norm(n));
th      = acos(evec'*r/(norm(evec)*norm(r)));
om      = acos(n'*evec/(norm(evec)*norm(n)));

oscul_elem = [a; e; i*180/pi; Om*180/pi; om*180/pi; th*180/pi];

if (h(1)<0)         Om = 2*pi-Om;    end;
if (evec(3)<0)      om = 2*pi-om;    end;
if (dot(r,evec)<0)  th = 2*pi-th;    end;

a_0     = a;    e_0     = e;    i_0     = i;    w_0     = om;   W_0     = Om;   th_0    = th;
a_00    = a;    e_00    = e;   i_00     = i;   w_00     = om;  W_00     = Om;  th_00    = th;

% Canuto
th_0      = atan2(norm(h)*dot(r,v),norm(h)*norm(h)-mu*norm(r));
u         = atan2(r(3)*sin(i_0)+cos(i_0)*(-r(1)*sin(W_0)+r(2)*cos(W_0)),r(1)*cos(W_0)+r(2)*sin(W_0));
w_0       = u-th_0;
E_0       = atan2(real(dot(r,v)/sqrt(mu*a_0)),real(1-norm(r)/a)); 

M_0       = E_0 - e_0 * sin(E_0);

h_0       = e_0*sin(w_0);
l_0       = e_0*cos(w_0);
lambda_0  = (w_0 + M_0);

h_00      = h_0;
l_00      = l_0;
lambda_00 = lambda_0;

%%%%%%%%%%%%
% Deltas Computation

G_2  = - 1081.874E-6  * ( ( RE * RE ) / ( a_0 * a_0 ) );

beta_0  = sin(i_0);
lambda_prime = 1.0 - 1.5 * G_2 * ( 3.0 - 4.0 * beta_0 );

Delta_a = ( ( -3.0 * a_0 * G_2 ) / ( 2.0 * lambda_prime ) ) * ( ...
                ( 2.0 - 3.5 * beta_0 * beta_0 ) * l_0 * cos( lambda_0 ) + ...
                ( 2.0 - 2.5 * beta_0 * beta_0 ) * h_0 * sin( lambda_0 ) + ...
                beta_0 * beta_0 * cos( 2.0 * lambda_0 ) + ...
                3.5 * beta_0 * beta_0 * ( l_0 * cos( 3.0 * lambda_0 ) + ...
                                            h_0 * sin( 3.0 * lambda_0 ) ) ) + ...
                0.75 * a_0 * G_2 * G_2 * beta_0 * beta_0 * ( ...
                ( 14.0 - 21.0 * beta_0 * beta_0 ) * cos( 2.0 * lambda_0 ) + ...
                                beta_0 * beta_0  *  cos( 4.0 * lambda_0 ) );

Delta_h = ( ( -3.0 * G_2 ) / ( 2.0 * lambda_prime ) ) * ( ...
                ( ( 1.0 - 1.75 * beta_0 * beta_0 ) * sin( lambda_0 ) ) + ...
                ( ( 1.0 - 3.0  * beta_0 * beta_0 ) * l_0 * sin( 2.0 * lambda_0 ) ) + ...
                ( (-1.5 + 2.0  * beta_0 * beta_0 ) * h_0 * cos( 2.0 * lambda_0 ) ) + ...
                ( ( 7.0 / 12.0 ) * beta_0 * beta_0 * sin( 3.0 * lambda_0 )  ) + ...
                ( ( ( 17.0 / 8.0 ) * beta_0 * beta_0 ) * ( l_0 * sin( 4.0 * lambda_0 ) - ...
                                                            h_0 * cos( 4.0 * lambda_0 ) ) ) );

Delta_l = ( ( -3.0 * G_2 ) / ( 2.0 * lambda_prime ) ) * ( ...
                ( ( 1.0 - 1.25 * beta_0 * beta_0 ) * cos( lambda_0 ) ) + ...
                ( ( 1.5 - 2.5  * beta_0 * beta_0 ) * l_0 * cos( 2.0 * lambda_0 ) ) + ...
                ( ( 2.0 - 1.5  * beta_0 * beta_0 ) * h_0 * sin( 2.0 * lambda_0 ) ) + ...
                ( ( 7.0 / 12.0 ) * beta_0 * beta_0 * cos( 3.0 * lambda_0 ) ) + ...
                ( ( ( 17.0 / 8.0 ) * beta_0 * beta_0 ) * ( l_0 * cos( 4.0 * lambda_0 ) + ...
                                                            h_0 * sin( 4.0 * lambda_0 ) ) ) );


Delta_i = ( (-3.0 * G_2 * beta_0 * sqrt( 1.0 - beta_0 * beta_0 ) ) / ...
                      (4.0 * lambda_prime) ) * ( ...
                -l_0 * cos( lambda_0 ) + h_0 * sin( lambda_0 ) + ...
                cos( 2.0 * lambda_0 ) + ( ( 7.0 / 3.0 ) *...
                                                l_0 * cos( 3.0 * lambda_0 ) + ...
                                                h_0 * sin( 3.0 * lambda_0 ) ) ) ;

Delta_W = - ( (-3.0 * G_2 * sqrt( 1.0 - beta_0 * beta_0 ) ) / ...
                      (4.0 * lambda_prime) ) * ( ...
                7.0 * l_0 * sin( lambda_0 ) + 5.0 * h_0 * cos( lambda_0 ) - ...
                sin( 2.0 * lambda_0 ) + ( ( 7.0 / 3.0 ) * ...
                                                -l_0 * sin( 3.0 * lambda_0 ) + ...
                                                h_0 * cos( 3.0 * lambda_0 ) ) ) ;

Delta_lambda = ( ( -3.0 * G_2 ) / ( 2.0 * lambda_prime ) ) * ( ...
                ( 10.0 - ( 119.0 / 8.0 ) * beta_0 * beta_0 ) * l_0 * sin( lambda_0 ) + ...
                ( -9.0 + ( 85.0 / 8.0  ) * beta_0 * beta_0 ) * h_0 * cos( lambda_0 ) + ...
                ( -0.5 + 2.0 * beta_0 * beta_0 ) * sin( 2.0 * lambda_0 ) + ...
                ( -( 7.0 / 6.0 ) + ( 119.0 / 24.0 ) * beta_0 * beta_0 ) * ...
                        ( l_0 * sin( 3.0 * lambda_0 ) - ...
                            h_0 * cos( 3.0 * lambda_0 ) ) - ...
                ( 3.0 - ( 21.0 / 4.0 ) * beta_0 * beta_0 )   * l_0 * sin( lambda_0 ) + ...
                ( 3.0 - ( 15.0 / 4.0 ) * beta_0 * beta_0 )   * h_0 * cos( lambda_0 ) - ...
                0.75 * beta_0 * beta_0 * sin( 2.0 * lambda_0 ) - ( ( 21.0 / 12.0 ) * ...
                            beta_0 * beta_0 * ( l_0 * sin( 3.0 * lambda_0) - ...
                                                h_0 * cos( 3.0 * lambda_0) ) )       );

Delta_a = real(Delta_a);                                
Delta_h = real(Delta_h);
Delta_l = real(Delta_l);
Delta_i = real(Delta_i);
Delta_W = real(Delta_W);
Delta_lambda = real(Delta_lambda);

sig=1;
if (i_0>pi/2) sigi=-1; else sigi = 1; end;

%%%%%%%%%%%%
% Mean Elements Computation

a_0 = a_00 - Delta_a;
i_0 = i_00 - Delta_i*sigi;
W_0 = W_00 - Delta_W*sig*sigi;
h_0 = h_00 - Delta_h;
l_0 = l_00 - Delta_l*sig;

e_0 = real(sqrt( ( l_0 * l_0 ) + ( h_0 * h_0 ) ));
w_0 = atan2( h_0 , l_0 );
M_0 = lambda_00 - Delta_lambda - w_0;
lambda_0 = lambda_00 - Delta_lambda;


% Iteratively determine the eccentric anomaly through Kepler's equation
maxError = 0.0001;
error = 1000.0;
E_0 = M_0;
while( error > maxError )
    E_temp = E_0;
    E_0 = E_0 - (E_0 - e_0 * sin(E_0) - M_0 ) / ( 1.0 -  e_0 * cos(E_0) );
    error = abs( E_0 - E_temp );
end;

E_0  = real(E_0);

th_0 = 2 * atan( sqrt((1+e_0)/(1-e_0))*tan(E_0/2) );
if (th_0<0) th_0 = th_0+2*pi; end;

r_0  = a_0*(1-e_0^2)/(1+e_0*cos(th_0));

lambda_0= w_0+M_0;
a_1  = a;	e_1  = e_0;	i_1  = i;	W_1  = Om;	th_1 = th_0;	w_1  = w_0;    akepler = a;    ikepler = i;

arbar  = a_1/norm(r);

alfaon = u;

J2k = 0.001082616;	earthradius = RE;
gamma2 = - 0.5 * J2k * (earthradius*earthradius/(a_1*a_1));
akeplermeanbl = akepler + a_1*gamma2 * ( (3*cos(i_1)*cos(i_1)-1)*(arbar*arbar*arbar-1/((1-e_1*e_1)*sqrt(1-e_1*e_1))) + 3*(1-cos(i_1)*cos(i_1))*arbar*arbar*arbar*cos(2*(alfaon)) );
akeplermeanble = akeplermeanbl + 18000*(e-e_0);


%%%%%%%%%%%%
% Show Results

semimajoraxis	= akeplermeanble;
eccentricity	= e_0;
inclination		= i_0;
RAAN		    = W_0;
periapsisargument= w_0;
trueanomaly		= th_0;
Ustinov_h		= h_0;
Ustinov_l		= l_0;
Ustinov_lambda 	= lambda_0;

mean_elem = [semimajoraxis; eccentricity; inclination; RAAN; periapsisargument; trueanomaly];
ustin_elem = [Ustinov_h; Ustinov_l; Ustinov_lambda];

end