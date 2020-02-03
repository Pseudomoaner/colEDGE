function g2 =  hessstrflt2(init_im, scale)

%Create normalized kernals - x^2, y^2 and xy (horizontal, vertical and diagonal edge detectors)
[y,x]=meshgrid(-3:6/scale:3,-3:6/scale:3);
g2a=(2*(x.^2)-1).*exp(-(x.^2+y.^2));
% imshow(g2a)
g2a=g2a/sum(abs(g2a(:)));

g2b=2*x.*y.*exp(-(x.^2+y.^2));
% figure
% imshow(g2b)
g2b=g2b/sum(abs(g2b(:)));

g2c=(2*(y.^2)-1).*exp(-(x.^2+y.^2));
% figure
% imshow(g2c)
g2c=g2c/sum(abs(g2c(:)));

%Filter image with edge detectors
g2a_rst=imfilter(init_im, g2a, 'symmetric', 'same');
g2(:,:,1)=g2a_rst;
g2b_rst=imfilter(init_im, g2b, 'symmetric', 'same');
g2(:,:,2)=g2b_rst;
g2c_rst=imfilter(init_im, g2c, 'symmetric', 'same');
g2(:,:,3)=g2c_rst;