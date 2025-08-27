function rxx = acorr(x)

n = length(x);
m = 2*n - 1;
rxx = zeros(1,m);
for i = 1:length(x)
    rxx(i) = dot(x(n-i+1:n), x(1:i).');
    rxx(m-i+1) = rxx(i)';
end

end
