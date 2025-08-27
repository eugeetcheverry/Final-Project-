function run_Lyapunov_p(ne,ext_fcn,t_start,h_norm,t_end,x_start,h,q,p_min,p_max,n);
hold on;
p_step=(p_max-p_min)/n
p=p_min;
while p<=p_max
    LE=FO_Lyapunov_p(ne,ext_fcn,t_start,h_norm,t_end,x_start,h,q,p);
    p=p+p_step
    plot(p,LE);
end