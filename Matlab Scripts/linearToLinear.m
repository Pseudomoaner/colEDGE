function y = linearToLinear(x,b,d,xTrans)

%trans should be the index of the transition point
y = zeros(size(x));

y(logical(x < xTrans)) = b*x(x < xTrans);
yTrans = b*xTrans;
intercept = yTrans - d*xTrans;

y(logical(x >= xTrans)) = intercept + d*x(x >= xTrans);