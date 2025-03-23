
% Vector spherical harmonic funtion modified from Hill's function
% V_lm for quantum numbers l and m. Angles theta and phi should be
% given in radians. 
%
function Wm = vsh_rscale(sphPoints,l,m)
r = sphPoints(:,1);
theta = sphPoints(:,2);
phi = sphPoints(:,3);

scale = r.^(-l-2);
scale_sph = (-1)^m*sqrt((2*l+1)*prod(1:(l-m))/(4*pi*prod(1:(l+m))));
scale_minus = 1/((-1)^(m-1)*sqrt((2*l+1)*prod(1:(l-(m-1)))/(4*pi*prod(1:(l+(m-1))))));
scale_plus = 1/((-1)^(m+1)*sqrt((2*l+1)*prod(1:(l-(m+1)))/(4*pi*prod(1:(l+(m+1))))));

 = spharm(theta,phi,l,m);
if m > -l
   Yminus = spharm(theta,phi,l,m-1);
else
   Yminus = 0;
end
if m < l
   Yplus = spharm(theta,phi,l,m+1);
else
   Yplus = 0;
end

dY = 0.5*scale_sph*((l+m)*(l-m+1)*scale_minus.*Yminus.*complex(cos(phi),sin(phi))-scale_plus.*Yplus.*complex(cos(-phi),sin(-phi)));  % dY/dtheta
if theta == 0 % To prevent division by zero
   theta = eps;
end
Wm(:,1) = scale*(-(l+1)).*Y;
Wm(:,2) = scale.*dY;
Wm(:,3) = scale.*(1i*m./sin(theta)).*Y;