function result = myCorr(X, Y)

meanX = mean(X);
meanY = mean(Y);
stdX = std(X);
stdY = std(Y);

for i = 1:1:length(X),
    X(i) = (X(i) - meanX)/stdX;
    Y(i) = (Y(i) - meanY)/stdY;
    mult = X(i) * Y(i);
end

result = sum(mult)/(length(X)-1);
end
