function [G, dGdR, dGdz, dGdzeta]=Green_series_zero(R,z,zeta,h)

% nondimensional variables
A = R/h;
v1 = z/h;
v2 = zeta/h;

% initialise outputs
gamma = 0.5772156649;
G = -(2*log(A/2) + 2*gamma + 1i*pi);
dGdR = - 1./(2*A);
dGdz = zeros(size(R));
dGdzeta = zeros(size(R));

% sum series
for m = 1:6
    c1 = cos(m*pi*v1);
    c2 = cos(m*pi*v2);
    s1 = sin(m*pi*v1);
    s2 = sin(m*pi*v2);
    K0 = besselk(0,m*pi*A);
    K1 = besselk(1,m*pi*A);
    
    G = G + 4*c1.*c2.*K0;
    dGdR = dGdR - m*pi*c1.*c2.*K1;
    dGdz = dGdz - m*pi*s1.*c2.*K0;
    dGdzeta = dGdzeta - m*pi*c1.*s2.*K0;
end

% dimensionalise outputs
G = G/h;
dGdR = 4*dGdR/h^2;
dGdz = 4*dGdz/h^2;
dGdzeta = 4*dGdzeta/h^2;

end