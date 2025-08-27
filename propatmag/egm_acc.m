function ac = egm_acc (x)

% ac = egm_acc (x)
% 
% 	To compute the acceleration due to Earth in geocentric terrestrial
%   coordinates. The accelaration is computed by lagrange expansion of the
%   Earth's gravity potential.
% 
% 	Inputs:
% 		x
% 			Position vector where the acceleration is calculated (m) in 
% 			geocentric terrestrial coordinates (3x1).
% 	Outputs:
% 		ac
% 			Acceleration vector in m/s2, in geocentrical coordinates (3x1).
%
%   The egm_acc function must be initialized before using it by calling
%   egm_read_dat function. The coeficient data are passed to this function
%   by global variables starting with 'egm_'. Please, avoid using variables 
%   starting with this string.
%
% 	Authors:
% 		Hélio Koiti Kuga		Fortran version
% 		Valdemir Carrara		Oct/2011	C version
%       Valdemir Carrara        July/2017    Matlab version
% */

global egm_order egm_length egm_conv_f egm_cc egm_sc
global egm_ae egm_gm egm_pn egm_qn egm_ip egm_nmax

    % auxiliary variables
	r     = norm(x);
	q     = egm_ae/r;
	t     = x(3)/r;		% sin (lat)
	u     = sqrt(1 - t*t);
	tf    = t/u;		% tan (lat)
% 	al    = atan2(x(2), x(1));
% 	sl    = sin(al);	% sin (long)
% 	cl    = cos(al);	% cos (long)
    sc    = sqrt(x(1)*x(1) + x(2)*x(2));
    if (sc == 0)
        sl  = 0;
        cl  = 1;
    else
        sl    = x(2)/sc;
        cl    = x(1)/sc;
    end
	gmr   = egm_gm/r;
      
	% summation initialization
	% omega = 0.d0
	vl    = 0.0;
	vf    = 0.0;
	vr    = 0.0;

	% store sectoral
	egm_pn(1)  = 1.0;
	egm_pn(2)  = 1.73205080756887730*u;	% sqrt(3) * cos (lat)
    egm_qn(1)  = 1.0;
	egm_qn(2)  = q;

	for m = 3: egm_nmax+1
		egm_pn(m)  = u*sqrt(1.0 + 0.50/(m-1))*egm_pn(m-1);
		egm_qn(m)  = q*egm_qn(m-1);
	end

	% initialize sin and cos recursions
	sm    = 0.0;
	cm    = 1.0;

	% outer n loop'
	for m = 1: egm_nmax+1
		% init 
		pnm    = egm_pn(m);			% m=n sectoral
		dpnm   = -(m-1)*pnm*tf;
		pnm1m  = pnm;
		pnm2m  = 0.0;

		% init  horner's scheme
		qc     = egm_qn(m)*egm_cc(egm_ip(m) + m);
		qs     = egm_qn(m)*egm_sc(egm_ip(m) + m);
		xc     = qc*pnm;
		xs     = qs*pnm;
		xcf    = qc*dpnm;
		xsf    = qs*dpnm;
		xcr    = m*qc*pnm;
		xsr    = m*qs*pnm;
        mm     = m - 1;

		% inner m loop 
        for n = m + 1: egm_nmax+1
            nn      = n - 1;
			anm     = sqrt(((nn + nn - 1.0)*(nn + nn + 1.0))/ ...
					((nn - mm)*(nn + mm)));
			bnm     = sqrt(((nn + nn + 1.0)*(nn + mm - 1.0)* ...
					(nn - mm - 1.0))/((nn - mm)*(nn + mm)*(nn + nn - 3.0)));
			fnm     = sqrt(((nn*nn - mm*mm)*(nn + nn + 1.0))/(nn + nn - 1.0));
			% recursion p and dp
			pnm     = anm*t*pnm1m - bnm*pnm2m;
			dpnm    = -nn*pnm*tf + fnm*pnm1m/u;		% signal opposite to paper
			% store
			pnm2m   = pnm1m;
			pnm1m   = pnm;

			% inner sum
            if (nn >= 2)
				qc     = egm_qn(n)*egm_cc(egm_ip(n) + m);
				qs     = egm_qn(n)*egm_sc(egm_ip(n) + m);
				xc     = (xc + qc*pnm);
				xs     = (xs + qs*pnm);
				xcf    = (xcf + qc*dpnm);
				xsf    = (xsf + qs*dpnm);
				xcr    = (xcr + (nn + 1.0)*qc*pnm);
				xsr    = (xsr + (nn + 1.0)*qs*pnm);
            end
        end

		% outer sum
		%	omega = omega + (xc*cm + xs*sm)  
		vl   = vl + mm*(xc*sm - xs*cm);
		vf   = vf + (xcf*cm + xsf*sm);
		vr   = vr + (xcr*cm + xsr*sm);

		% sin and cos recursions to next m
		cml  = cl*cm - sm*sl;
		sml  = cl*sm + cm*sl;
		cm   = cml;			% save to next m
		sm   = sml;			% save to next m
	end

	% finalization, include n=0 (p00=1), 
	% for n=1 all terms are zero: c,s(1,1), c,s(1,0) = 0

	% potential
	% omega =  gmr * (1.d0+omega)

	% gradient 
	vl    = -gmr*egm_conv_f*vl;
	vf    =  gmr*egm_conv_f*vf;
	vr    = -(gmr/r)*(1.0 + egm_conv_f*vr);
	
	% body x, y, z accelerations 
	ac      = [u*cl*vr - t*cl*vf/r - sl*vl/(u*r); ...
			u*sl*vr - t*sl*vf/r + cl*vl/(u*r); ...
			t*vr + u*vf/r];
return
        


