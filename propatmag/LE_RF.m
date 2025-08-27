function f=LE_RF(t,x)
%Output data must be a column vector
f=zeros(size(x));
%variables allocated to the variational equations
X= [x(4), x(7), x(10);
x(5), x(8), x(11);
x(6), x(9), x(12)];
%RF equations
f(1)=x(2)*(x(3)-1+x(1)*x(1))+0.1*x(1);
f(2)=x(1)*(3*x(3)+1-x(1)*x(1))+0.1*x(2);
f(3)=-2*x(3)*(0.98+x(1)*x(2));
%Jacobian matrix
J=[2*x(1)*x(2)+0.1, x(1)*x(1)+x(3)-1, x(2);
-3*x(1)*x(1)+3*x(3)+1,0.1,3*x(1);
-2*x(2)*x(3),-2*x(1)*x(3),-2*(x(1)*x(2)+0.98)
];
%Righthand side of variational equations
f(4:12)=J*X; % To be modified if ne>3