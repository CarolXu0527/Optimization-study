[x,y,z] = meshgrid(-4:.2:4,-4:.2:4,-4:.2:4);
v = 2*(x.^2) + 2*(y.^2) + z.^4 + x.*8 + y.*4 + z.*10 + 2*exp(-(z./4).^2);
xslice = [-2,0,2]; 
yslice = 2; 
zslice = [-2,0,2];
slice(x,y,z,v,xslice,yslice,zslice)
colormap hsv
