function [AJ2] = fAJ2(mean_Ustinov_chief, mean_elem_chief)

    l_0cc = mean_Ustinov_chief(2);
    h_0cc = mean_Ustinov_chief(1);
    i_0cc = mean_elem_chief(3);
    a_0cc = mean_elem_chief(1);

    J2k = 0.001082616;
    mu = 3.986004418e14;
    RE = 6378135;

    % Término de drift
    ex=l_0cc;
    ey=h_0cc;
    P=5*(cos(i_0cc))^2-1;
    Q=3*(cos(i_0cc))^2-1;
    R=sin(2*i_0cc);
    S=(sin(i_0cc))^2;
    eta=sqrt(1-ex^2-ey^2);
    M=1+eta;
    G=4+3*eta;
    H=1/(eta^2);
    gamma=(J2k/2)*((RE/a_0cc)^2)/eta^4;
    kappan=(gamma/1.5)/(a_0cc^(7/2)*eta^4);
    n=sqrt(mu)/(sqrt(a_0cc))^3;
    T=sin(i_0cc)^2;
    AJ2=[[0                0        0               0                0  0];
        [-3.5*M*Q          0      ex*G*H*Q       ey*G*H*Q         -G*S  0];
        [ 3.5*ey*P         0    -4*ex*ey*H*P  -(1+4*ey*ey*H)*P  5*ey*S  0];
        [-3.5*ex*P         0  (1+4*ex*ex*H)*P    4*ex*ey*H*P   -5*ex*S  0];
        [ 0                0        0               0                0  0];
        [ 3.5*S            0    -4*ex*H*S       -4*ey*H*S          2*T  0]];
    AJ2  = 1.5*kappan*AJ2;
       
    AK   = [[0 0 0 0 0 0];[-1.5*sqrt(mu/mean_elem_chief(1)^3) 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0];[0 0 0 0 0 0]];
    
    %AJ2 = AJ2 + AK;

