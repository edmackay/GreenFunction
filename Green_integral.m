function [G,dGdR,dGdz,dGdzeta]=Green_integral(K,k,R,z,zeta,h,Xa,Y1,Y2,Flookup)

% define vertical distance terms
v1=abs(z-zeta);
v2=2*h+z+zeta;
v3=abs(z+zeta);
v4=2*h+z-zeta;
v5=2*h-z+zeta;
v6=4*h+z+zeta;

% define non-dimensional variables
V3=K*v3;
V4=K*v4;
V5=K*v5;
V6=K*v6;
X=K*R;
H=K*h;

% calculate integrals from Chebyshev approximation
[M1,MA1,MB1]=Mfun(Xa,Y1,H);
[M2,MA2,MB2]=Mfun(Xa,Y2,H);
 
% calculate infinite depth parts
[FXV3,dFdX_V3,dFdV_V3]=FXY_interp(X,V3,Flookup);
[FXV4,dFdX_V4,dFdV_V4]=FXY_interp(X,V4,Flookup);
[FXV5,dFdX_V5,dFdV_V5]=FXY_interp(X,V5,Flookup);
[FXV6,dFdX_V6,dFdV_V6]=FXY_interp(X,V6,Flookup);

% distance terms
R2=R.^2;
r1inv=(R2+v1.^2).^(-0.5);
r2inv=(R2+v2.^2).^(-0.5);
r3inv=(R2+v3.^2).^(-0.5);
r4inv=(R2+v4.^2).^(-0.5);
r5inv=(R2+v5.^2).^(-0.5);
r6inv=(R2+v6.^2).^(-0.5);
r1inv3=r1inv.^3;
r2inv3=r2inv.^3;
r3inv3=r3inv.^3;
r4inv3=r4inv.^3;
r5inv3=r5inv.^3;
r6inv3=r6inv.^3;

% Green function and derivatives
K2=K^2;
s1=sign(z-zeta);
G = (r1inv + r2inv + r3inv + r4inv + r5inv + r6inv) + K *(FXV3 + FXV4 + FXV5 + FXV6) + (M1 + M2)/h;
dGdR = -R.*(r1inv3 + r2inv3 + r3inv3 + r4inv3 + r5inv3 + r6inv3) + K2 *(dFdX_V3 + dFdX_V4 + dFdX_V5 + dFdX_V6) + (MA1 + MA2)/(h.^2);
dGdz = (-s1.*v1.*r1inv3 - v2.*r2inv3 + v3.*r3inv3 - v4.*r4inv3 + v5.*r5inv3 - v6.*r6inv3) + K2 *(- dFdV_V3 + dFdV_V4 - dFdV_V5 + dFdV_V6) + (s1.*MB1 + MB2)/(h.^2);
dGdzeta = (s1.*v1.*r1inv3 - v2.*r2inv3 + v3.*r3inv3 + v4.*r4inv3 - v5.*r5inv3 - v6.*r6inv3) + K2 *(- dFdV_V3 - dFdV_V4 + dFdV_V5 + dFdV_V6) + (-s1.*MB1 + MB2)/(h.^2);

% Imaginary part
J0=besselj(0,k*R);
J1=besselj(1,k*R);
ekz=exp(k*z);
ekzh=exp(-k*(z+2*h));
ekze=exp(k*zeta);
ekzeh=exp(-k*(zeta+2*h));
C0=k^2/((k^2-K^2)*h+K);
a=2*pi*1i*C0/((1+exp(-2*k*h))^2);

Im_G=-a*J0.*(ekz+ekzh).*(ekze+ekzeh);
Im_dGdR=a*k*J1.*(ekz+ekzh).*(ekze+ekzeh);
Im_dGdz=-a*k*J0.*(ekz-ekzh).*(ekze+ekzeh);
Im_dGdzeta=-a*k*J0.*(ekz+ekzh).*(ekze-ekzeh);

% Add in imaginary component
G=G+Im_G;
dGdR=dGdR+Im_dGdR;
dGdz=dGdz+Im_dGdz;
dGdzeta=dGdzeta+Im_dGdzeta;

end