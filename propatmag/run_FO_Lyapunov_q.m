function run_FO_Lyapunov_q(ne,ext_fcn,t_start,h_norm,t_end,x_start,h,q_min,q_max,n);
hold on;
q_step=(q_max-q_min)/n;
q=q_min;
while q<q_max
    [t,LE]=FO_Lyapunov_q(ne,ext_fcn,t_start,h_norm,t_end,x_start,h,q);
    q=q+q_step;
    fprintf('q=%10.4f\n %10.4f', q);
    plot(q,LE);
end