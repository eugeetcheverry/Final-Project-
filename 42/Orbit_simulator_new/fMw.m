function [Mw] = fMw8(mean_Ustinov, mean_elem_chief, lambdaexp, nopR)

    mu = 3.986004418e14;

    u0   = mean_Ustinov(3);
    u1   = mean_Ustinov(3) + pi/2;
    u2   = mean_Ustinov(3) + pi;
    u3   = mean_Ustinov(3) + 3*pi/2;
    u01  = mean_Ustinov(3) + pi/4;
    u11  = mean_Ustinov(3) + 3*pi/4;
    u21  = mean_Ustinov(3) + 5*pi/4;
    u31  = mean_Ustinov(3) + 7*pi/4;

    B0   = 1/(mean_elem_chief(1)*sqrt(mu/mean_elem_chief(1)^3)) *[[0 2 0];[-2 0 0];[sin(u0)  2*cos(u0)  0];[-cos(u0)  2*sin(u0)  0];[0 0 cos(u0)]; [0 0 sin(u0)]];
    B1   = 1/(mean_elem_chief(1)*sqrt(mu/mean_elem_chief(1)^3)) *[[0 2 0];[-2 0 0];[sin(u1)  2*cos(u1)  0];[-cos(u1)  2*sin(u1)  0];[0 0 cos(u1)]; [0 0 sin(u1)]];
    B2   = 1/(mean_elem_chief(1)*sqrt(mu/mean_elem_chief(1)^3)) *[[0 2 0];[-2 0 0];[sin(u2)  2*cos(u2)  0];[-cos(u2)  2*sin(u2)  0];[0 0 cos(u2)]; [0 0 sin(u2)]];
    B3   = 1/(mean_elem_chief(1)*sqrt(mu/mean_elem_chief(1)^3)) *[[0 2 0];[-2 0 0];[sin(u3)  2*cos(u3)  0];[-cos(u3)  2*sin(u3)  0];[0 0 cos(u3)]; [0 0 sin(u3)]];
    B01  = 1/(mean_elem_chief(1)*sqrt(mu/mean_elem_chief(1)^3)) *[[0 2 0];[-2 0 0];[sin(u01) 2*cos(u01) 0];[-cos(u01) 2*sin(u01) 0];[0 0 cos(u01)];[0 0 sin(u01)]];
    B11  = 1/(mean_elem_chief(1)*sqrt(mu/mean_elem_chief(1)^3)) *[[0 2 0];[-2 0 0];[sin(u11) 2*cos(u11) 0];[-cos(u11) 2*sin(u11) 0];[0 0 cos(u11)];[0 0 sin(u11)]];
    B21  = 1/(mean_elem_chief(1)*sqrt(mu/mean_elem_chief(1)^3)) *[[0 2 0];[-2 0 0];[sin(u21) 2*cos(u21) 0];[-cos(u21) 2*sin(u21) 0];[0 0 cos(u21)];[0 0 sin(u21)]];
    B31  = 1/(mean_elem_chief(1)*sqrt(mu/mean_elem_chief(1)^3)) *[[0 2 0];[-2 0 0];[sin(u31) 2*cos(u31) 0];[-cos(u31) 2*sin(u31) 0];[0 0 cos(u31)];[0 0 sin(u31)]];    

    M=[B0 B01 B1 B11 B2 B21 B3 B31];
    
    w0=exp([1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8]*lambdaexp-lambdaexp);

    Mw = M.*w0;

    %w00=exp([1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10 11 11 11 12 12 12 13 13 13 14 14 14 15 15 15 16 16 16]*lambdaexp-lambdaexp);
    %Mw = [M M].*w00;
    
    if nopR==1
      Mw = M(:,[2 3 5 6 8 9 11 12 14 15 17 18 20 21 23 24]).*w0([[2 3 5 6 8 9 11 12 14 15 17 18 20 21 23 24]]);
    end;