function y = exponentialToLinear(x,a,b,d,xTrans)

%trans should be the index of the transition point
y = zeros(size(x));

y(logical(x < xTrans)) = a * exp(b*x(x < xTrans));
yTrans = a * exp(b*xTrans);
intercept = yTrans - d*xTrans;

y(logical(x >= xTrans)) = intercept + d*x(x >= xTrans);