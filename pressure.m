function p = pressure(U)

global gamma

p = (gamma - 1)*(U(3) - 0.5*U(2).^2./U(1));

end