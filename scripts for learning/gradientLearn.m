x = -3:0.2:3;
y = x';
f = x.^2 .* y.^3;
surf(x,y,f)
xlabel('x')
ylabel('y')
zlabel('z')
[fx,fy] = gradient(f,0.2);
%%
x0 = 1;
y0 = -2;
t = (x == x0) & (y == y0);
indt = find(t);
f_grad = [fx(indt) fy(indt)]