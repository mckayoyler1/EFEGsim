function vlm = getvlm(l,m)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

syms theta phi real

er = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)];
et = [cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)];
ep = [-sin(phi), cos(phi), 0];
Ylm = symbolicYlm(l,m);
vlm = -(l+1)*Ylm*er + diff(Ylm,theta)*et + 1i*m*Ylm*ep/sin(theta);
end