
function xdot  = HH3(t, y, current,  P)
    inj= interp1(current(:,1),current(:,2),t);
    % n
    tt = tanh((y(1)-P(8))/P(9));
    ninf = 0.5 + 0.5*tt;
    tn = P(10) + P(11)*(1-tt^2);
    %m
    tt = tanh((y(1)-P(12))/P(13));
    minf = 0.5 + 0.5*tt;
    tm = P(14) + P(15)*(1-tt^2);
    %h
    tt = tanh((y(1)-P(16))/P(17));
    hinf = 0.5 + 0.5*tt;
    th = P(18) + P(19)*(1-tt^2);
    %P = [1 1.0 2 120 3 115 4 20 5 -12 6 0.3 7 10.6 8 10 9 30 10 1 11 5 12,25
    %.. 13,15 14,0.1 15,0.4 16,5 17,-15 18,1 19,7];
    xdot = [ P(1)*inj + y(3)^3 *y(4)*P(2)*(P(3)-y(1)) + P(4)*y(2)^4*(P(5)-y(1)) + P(6)*(P(7)-y(1));...
            (ninf - y(2))/tn;...
            (minf - y(3))/tm;...
            (hinf - y(4))/th];
    
end
