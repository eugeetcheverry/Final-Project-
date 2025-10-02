function [x_corr, sigma_corr, phikm1] = ekf_rmm(w, rmm, u, iner, earth_field_b, sigma, z, Q, R, dt)
    %Calculo de deriva y jacobianos
    
    %Ver que u enmascara una dependencia de w que podría verse en la matriz
    %podría reemplazarse u por la ley de control y solo tener en cuenta la
    %dependencia con x y el tiempo
    xdot_hat = [inv(iner)*(cross(u, earth_field_b) + cross(rmm, earth_field_b) - cross(w, iner*w));
                0;
                0;
                0];

    %xpunt = f(x,u) = f(x,t,n)
    %Linealizar en torno a 0
    w_skew = Skew(w);
    Jw = iner*w;
    Jw_skew = Skew(Jw);
    B_skew = Skew(earth_field_b);
    
    A = [iner\(Jw_skew - w_skew*iner) -inv(iner)*B_skew; %A
            zeros(3) zeros(3)];
        
    %phikm1 = expm(A*dt)
    phikm1 = eye(6) + A*dt;
    
    %Al eliminar la dependencia de u, B podría ser una identidad
    B = [-inv(iner)*B_skew; %B
            zeros(3)];
    gammkm1 = dt*phikm1*B;
    
    %Pasito de prediccion
    x_hat = [w; rmm] + xdot_hat*dt;
    sigma_hat = phikm1*sigma*phikm1' + gammkm1*Q*gammkm1';
    
    %Pasito de correccion
    H = [eye(3) zeros(3)];
    
    %Cambiado respecto del paper, se usa la sensibilidad (robotica) en vez
    %de solo la R
    K = sigma_hat*H'/(H*sigma_hat*H' + R);
    
    x_corr = x_hat + K*(z - H*x_hat);
    
    %Ecuacion de robotica
    sigma_corr = (eye(6) - K*H)*sigma_hat;
    
end
