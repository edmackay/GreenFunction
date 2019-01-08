function [G,dGdR,dGdz,dGdzeta]=Green_integral_infty(R,z,zeta,h,X,Y1,Y2)

% define vertical distance terms
v1=abs(z-zeta);
v2=2*h+z+zeta;
v3=abs(z+zeta);
v4=2*h+z-zeta;
v5=2*h-z+zeta;
v6=4*h+z+zeta;

% calculate integrals from Chebyshev approximation
[M1,MA1,MB1]=Mfun(X,Y1,Inf);
[M2,MA2,MB2]=Mfun(X,Y2,Inf);

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
s1=sign(z-zeta);
G = (r1inv + r2inv - r3inv - r4inv - r5inv - r6inv) + (M1 + M2)/h;
dGdR = R.*(-r1inv3 - r2inv3 + r3inv3 + r4inv3 + r5inv3 + r6inv3) + (MA1 + MA2)/(h.^2);
dGdz = (-s1.*v1.*r1inv3 - v2.*r2inv3 - v3.*r3inv3 + v4.*r4inv3 - v5.*r5inv3 + v6.*r6inv3) + (s1.*MB1 + MB2)/(h.^2);
dGdzeta = (s1.*v1.*r1inv3 - v2.*r2inv3 - v3.*r3inv3 - v4.*r4inv3 + v5.*r5inv3 + v6.*r6inv3) + (-s1.*MB1 + MB2)/(h.^2);

end