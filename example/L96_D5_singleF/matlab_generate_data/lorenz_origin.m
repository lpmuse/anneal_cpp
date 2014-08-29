function dxdt = lorenz_origin(t,x)

N = length(x)-1;

dxdt = zeros(N+1,1);
%v = 8.17;
for i = 1:N
    dxdt(i) = x(mymod(i-1,N))*(x(mymod(i+1,N))-x(mymod(i-2,N))) - x(i) + x(N+1);
end
dxdt(N+1) = 0;
end

function re=mymod(x,y)
    re = mod(x,y);
    if re==0
        re= re+y;
    end
end
        
    