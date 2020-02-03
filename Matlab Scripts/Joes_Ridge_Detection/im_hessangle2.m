function [A, R] = im_hessangle2(im, scale)

%Run ridge detection filtering stages - g2 corresponds to the values of the Hessian matrix for each point in the greyscale intensity matrix.
g2 = im_hessstrflt2(im, scale);

tmp = sqrt((g2(:,:,1)- g2(:,:,3)).*(g2(:,:,1)- g2(:,:,3)) + 4 * g2(:,:,2).*g2(:,:,2)); 
%Components of the eigenvectors associated with the Hessian g2 (in polar coordinates)
eigvalue1 = (g2(:,:,1) + g2(:,:,3) + tmp)/2; %Formula for eigenvalues that falls out from working through the 2x2 case. Not Gamma normalized in this case (ie. is not normalized according to scale - see eqn 47 of Lindeberg 1996)
eigvalue2 = (g2(:,:,1) + g2(:,:,3) - tmp)/2;
eigangle1 = atan(g2(:,:,2)./(eigvalue1 - g2(:,:,3)+realmin));
eigangle2 = atan(g2(:,:,2)./(eigvalue2 - g2(:,:,3)+realmin));

%Calculate magnitude and angle of the maximum eigenvalues/vectors of the Hessian at each point in the image
A = eigvalue1.*(abs(eigvalue1)>=abs(eigvalue2)) + eigvalue2.*(abs(eigvalue1)<abs(eigvalue2));
R = eigangle2.*(abs(eigvalue1)>=abs(eigvalue2)) + eigangle1.*(abs(eigvalue1)<abs(eigvalue2));

%Do some angle rescaling to degrees
R(R < 0) = R(R < 0) + pi;
R = R*180/ pi;
