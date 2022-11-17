function [a1,a2,a3] = eigenvalues(rho,rhou,E)

global gamma

u = rhou/rho;
p = (gamma - 1)*(E - 0.5*rho*u^2);
c = sqrt(abs(gamma*p/rho));

a1 = u + c; a2 = u; a3 = u - c;

end