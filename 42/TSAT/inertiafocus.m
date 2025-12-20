% FOCUS Moment of Inertia

mb = 120;
mp = 1.157;
Lx = 0.5;
Ly = 0.75;
Lz = 0.75;
Lp = 4.9;
Ixx = 1/12*mb*(Ly^2+Lz^2) + 7*mp*(Ly/2)^2;
Iyy = 1/12*mb*(Lx^2+Lz^2) + 7/12*mp*(Lp^2+Lz^2);
Izz = 1/12*mb*(Lx^2+Ly^2) + 7/12*mp*(Lp^2);

