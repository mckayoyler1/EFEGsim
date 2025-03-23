function Vlm = getVlm(sphPos, L)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
[n,m] = size(sphPos);
k = (L+1)^2;
Vlm = zeros(n,m,k);
dirs = [sphPos(:,2),sphPos(:,3)];
Y_N = getSH2L(L, dirs, 'complex');
phi = dirs(:,2);
theta = dirs(:,1);

    % Loop over degrees l = 0 to L
    for l = 0:L
   
        for m = -l:l
            Wlm = zeros(size(sphPos));
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
            Wlm(:,1) = (-(l+1)).*Y;
            Wlm(:,2) = dY;
            Wlm(:,3) = (1i*m./sin(theta)).*Y;
            % index = l^2 + l + m + 1
            idx = l^2 + l + m + 1;
            Vlm(:,:,idx) = Wlm;
        end

    end
end