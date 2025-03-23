function V = evaluateSymbV(pos, symbV)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
syms x y z
Vfunc = matlabFunction(symbV, 'Vars', [x, y, z]);
xs = pos(:,1);
ys = pos(:,2);
zs = pos(:,3);
V = Vfunc(xs,ys,zs);
end