function out = Skew(in)
%  Skew transforms between matrix and vector representations.
%      X = Skew(x) returns skew-symmetric matrix of 3D vector x.
%      
%      x = Skew(X) returns vector from skew-symmetric 3by3 matrix X.
%      
%      Q = Skew(q) returns tensor equivalent of quaternion q.
%      
%      q = Skew(Q) returns quaternion from tensor equivalent Q.
%      
%      Works with symbolic variables.

if size(in,1) == 3
    if size(in,2) == 1
        x = in;
        X = [    0 -x(3)  x(2);
              x(3)     0 -x(1);
             -x(2)  x(1)     0];
        out = X;
    else
        X = in;
        x = [X(3,2);
            X(1,3);
            X(2,1)];
        out = x;
    end
else
    if size(in,2) == 1
        [q0,q1,q2,q3] = in';
        Q = [ q0  q1  q2  q3;
             -q1  q0 -q3  q2;
             -q2  q3  q0 -q1;
             -q3 -q2  q1  q0];
        out = Q;
    else
        Q = in;
        q = [Q(1,1);
             Q(1,2);
             Q(1,3);
             Q(1,4)];
        out = q;
    end
end