function [G, dGdR, dGdz, dGdzeta]=Green_series(K,k0,R,z,zeta,h)

% compute coefficient for terms in k0
C0=k0^2/((k0^2-K^2)*h+K);

% compute terms in k0
Y0=bessely(0,k0*R);
Y1=bessely(1,k0*R);
J0=besselj(0,k0*R);
J1=besselj(1,k0*R);
ekz=exp(k0*z);
ekzh=exp(-k0*(z+2*h));
ekze=exp(k0*zeta);
ekzeh=exp(-k0*(zeta+2*h));
a=2*pi*C0/((1+exp(-2*k0*h))^2);

G=-a*(Y0+1i*J0).*(ekz+ekzh).*(ekze+ekzeh);
dGdR=a*k0*(Y1+1i*J1).*(ekz+ekzh).*(ekze+ekzeh);
dGdz=-a*k0*(Y0+1i*J0).*(ekz-ekzh).*(ekze+ekzeh);
dGdzeta=-a*k0*(Y0+1i*J0).*(ekz+ekzh).*(ekze-ekzeh);

% compute complex wave numbers
N=20;
kn=wavenumber_imag(K,h,N);
a=kn.^2+K^2;
Cn=a./(a*h-K);

% add terms whilst difference is >1e-7
n=0;
dif=1;
while dif>1e-7 && n<N
    n=n+1;
    
    % preliminary definitions
    K0=besselk(0,kn(n)*R);
    K1=besselk(1,kn(n)*R);
    cz=cos(kn(n)*(z+h));
    sz=sin(kn(n)*(z+h));
    cze=cos(kn(n)*(zeta+h));
    sze=sin(kn(n)*(zeta+h));
    
    % terms in summation
    term_G=4*Cn(n)*cz.*cze.*K0;
    term_GR=-4*Cn(n)*kn(n)*cz.*cze.*K1;
    term_Gz=-4*Cn(n)*kn(n)*sz.*cze.*K0;
    term_Gze=-4*Cn(n)*kn(n)*cz.*sze.*K0;
    
    % Green function and derivatives
    G=G+term_G;
    dGdR=dGdR+term_GR;
    dGdz=dGdz+term_Gz;
    dGdzeta=dGdzeta+term_Gze;
    
    % maximum difference
    dif=max([abs(term_G);abs(term_GR);abs(term_Gz);abs(term_Gze)]);
end

end