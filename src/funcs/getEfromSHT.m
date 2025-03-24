function data = getEfromSHT(data)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here: 

sphpos = data.sphpos;
theta = sphpos(:,2);
phi = sphpos(:,3);
Y_N = data.basis;

[n,m] = size(sphpos);
L = data.expansionOrder;
k = (L+1)^2;
Vlm = zeros(n,m,k);
vr_N = zeros([n,k]);
vt_N = zeros([n,k]);
vp_N = zeros([n,k]);

% Loop over degrees l = 0 to L
for l = 0:L

    for m = -l:l
        Y = getSHEntry(Y_N, l,m);
        if m > -l
            Yminus = getSHEntry(Y_N, l,m-1);
        else
            Yminus = 0;
        end
        if m < l
            Yplus = getSHEntry(Y_N,l,m+1);
        else
            Yplus = 0;
        end

        scale_sph = (-1)^m*sqrt((2*l+1)*prod(1:(l-m))/(4*pi*prod(1:(l+m))));
        scale_minus = 1/((-1)^(m-1)*sqrt((2*l+1)*prod(1:(l-(m-1)))/(4*pi*prod(1:(l+(m-1))))));
        scale_plus = 1/((-1)^(m+1)*sqrt((2*l+1)*prod(1:(l-(m+1)))/(4*pi*prod(1:(l+(m+1))))));
        dY = 0.5*scale_sph*((l+m)*(l-m+1)*scale_minus.*Yminus.*complex(cos(phi),sin(phi))-scale_plus.*Yplus.*complex(cos(-phi),sin(-phi)));
       
        idx = l^2 + l + m + 1;
        vr_N(:,idx) = (-(l+1)).*Y;
        vt_N(:,idx) = dY;
        vp_N(:,idx) = (1i*m./sin(theta)).*Y;

        Vlm(:,:,idx) = [vr_N(:,idx),vt_N(:,idx),vp_N(:,idx)];
    end
end
Rscale = data.Escaling;
vr_N_scaled = vr_N ./ Rscale;
vt_N_scaled = vt_N ./ Rscale;
vp_N_scaled = vp_N ./ Rscale;
eps0 = Constants.eps0;
A_N = data.coeffs;
Er = -1/eps0 * vr_N_scaled * A_N;
Et = -1/eps0 * vt_N_scaled * A_N;
Ep = -1/eps0 * vp_N_scaled * A_N;
complexE       = struct();
complexE.r     = Er;
complexE.theta = Et;
complexE.phi   = Ep;
data.E.sphcomplex = complexE;
data.E.sph = real([Er,Et,Ep]);
data.E.cart = convertSphericalVecToCartesian(sphpos, data.E.sph);


% Vlms = getWsh(sphpos,data.expansionOrder);
% shapedScaledCoeffs = reshape(data.scaledCoeffs,1,1,[]);
% scaledCoeffs = data.scaledCoeffs;
% scaledCoeffs(abs(scaledCoeffs) < 1e-12) = 0;
% complex_vals = sum(shapedScaledCoeffs.* Vlms,3);
% E.test_check = reshape(scaledCoeffs,1,1,[]).* Vlms;
% 
% eps0 = Constants.eps0;
% E.sphComplex = -eps0* complex_vals;
% E.sphReal = -eps0 * real(complex_vals);
% E.Vlms = Vlms;
% E.scaledCoeffs = scaledCoeffs;
% E.shapedScaledCoeffs = shapedScaledCoeffs;
% E.cartReal = convertSphericalVecToCartesian(sphpos, E.sphReal);
% data.E = E;
end
