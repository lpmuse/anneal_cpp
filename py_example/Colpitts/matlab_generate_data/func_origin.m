function dxdt = colpitts_origin(t,x)

dxdt = zeros(7,1);
dxdt(1) = x(4)*x(2);
dxdt(2) = -x(5)*(x(1)+x(3)) - x(6)*x(2);
dxdt(3) = x(7)*(x(2)+1-exp(-x(1)));
dxdt(4:7) = 0;
end

        
    
