function p = plotSlice(pos,vals)
%UNTITLED22 Summary of this function goes here
%   Detailed explanation goes here
sphPos = toSphere(pos);
keep   = any(sphPos(:,3)==0, 2);
sphPos = sphPos(keep,:);
vals   = vals(keep,:);
p = plot(sphPos(:,2),vals);
end