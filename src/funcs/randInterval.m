function nums = randInterval(N, M, range)
%randInterval generates an NxM matrix of random values on the interval (a,b)
a = range(1);
b = range(2);
nums = a + (b-a).*rand(N,M);
end