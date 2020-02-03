function [Magnitude, Angle, Width] = scalablehess2(im, scale_range, ridge_property)

if (nargin == 2)
    ridge_property = 'dark';
end

dim = 0;
for scale = scale_range
    dim = dim + 1;
    [Magnitude(:,:,dim), Angle(:,:,dim)] = im_hessangle2(im, scale);
end

if strcmp(ridge_property, 'bright');
    Magnitude = -Magnitude;
elseif strcmp(ridge_property, 'dark');
    ;
else
    error('ridge_property error @ scalablehess2');  
end

[M, N] = size(im);
[Magnitude, Max_Idx] = max(Magnitude, [], 3); %Maximal eigenvalue of the Hessian across spatial scales
Max_Idx2 = (Max_Idx-1)*M*N + reshape(1:M*N, M, N); %Linear index of the maximal eigenvalue
Angle = Angle(Max_Idx2);
Width = scale_range(Max_Idx);

Magnitude(Magnitude<0) = 0;
Magnitude = (Magnitude-min(Magnitude(:)))/(max(Magnitude(:))-min(Magnitude(:)) + realmin); %Set to vary between 0 and 1

