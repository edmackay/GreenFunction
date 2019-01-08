function [G, dGdR, dGdz, dGdzeta]=Green_series_infty(R,z,zeta,h)

% nondimensional variables
A = R/h;
v1 = z/h;
v2 = zeta/h;

% initialise outputs
G = zeros(size(R));
dGdR = zeros(size(R));
dGdz = zeros(size(R));
dGdzeta = zeros(size(R));

% sum series
for m = 0.5 : 1 : 10.5
    c1 = cos(m*pi*v1);
    c2 = cos(m*pi*v2);
    s1 = sin(m*pi*v1);
    s2 = sin(m*pi*v2);
    K0 = besselk(0,m*pi*A);
    K1 = besselk(1,m*pi*A);
    
    G = G + 4*s1.*s2.*K0;
    dGdR = dGdR - 4*m*pi*s1.*s2.*K1;
    dGdz = dGdz + 4*m*pi*c1.*s2.*K0;
    dGdzeta = dGdzeta + 4*m*pi*s1.*c2.*K0;
end

% dimensionalise outputs
G = G/h;
dGdR = dGdR/h^2;
dGdz = dGdz/h^2;
dGdzeta = dGdzeta/h^2;

end