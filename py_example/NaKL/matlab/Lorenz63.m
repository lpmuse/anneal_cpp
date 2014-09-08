function ydot= Lorenz63(t, y)
    ydot = 0.08*[16*(y(2)-y(1))
        -y(1)*y(3) + 45.92*y(1) - y(2)
        y(1)*y(2) - 4*y(3)];

end

