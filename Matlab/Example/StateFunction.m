function dx = StateFunction(x, parameters)
%Function that describes the pendulum state
dx = zeros(2, 1);
dx(1, 1) = x(1) + x(2)*parameters(1, 1);
dx(2, 1) = x(2) - parameters(2, 1)*sin(x(1))*parameters(1, 1);

end

