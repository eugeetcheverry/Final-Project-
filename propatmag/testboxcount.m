% Test Box Count

% Plots the box-count of a vector containing randomly-distributed
% 0 and 1. This set is not fractal: one has N = R^-2 at large R,
% and N = cste at small R.
c = (rand(1,2048)<0.2);
boxcount(c);
figure, boxcount(c, 'slope');

% Plots the box-count and the fractal dimension of a 2D fractal set
% of size 512^2 (obtained by RANDCANTOR), with fractal dimension
% DF = 2 + log(P) / log(2) = 1.68 (with P=0.8).
c = randcantor(0.8, 512, 2);
boxcount(c);
figure, boxcount(c, 'slope');