function [FXY,dFdX,dFdY]=FXY_interp(X,Y,lookup)

% Note that a factor of (X^2+Y^2)^(-0.5) has been removed compared to Newman's definition

% define regions
xmin=1e-6;
r00 = X==0 & Y==0;
r1 = xmin<X & X<=3 & Y<=4;
r2 = X>3 & Y<=4;
r3 = xmin<X & X<=3 & Y>4;
r4 = X>3 & Y>4;
X0_1 = X<=xmin & Y<0.25;
X0_2 = X<=xmin & Y>=0.25;

% intialise variables
FXY=0*X;
dFdX=0*X;
dFdY=0*X;

% Region 0
if sum(r00)>0
    FXY(r00)=inf;
    dFdX(r00)=-inf;
    dFdY(r00)=-inf;
end

% Region 1
if sum(r1)>0
    % pre-compute some terms
    Xr=X(r1);
    Yr=Y(r1);   
    R=sqrt(Xr.^2+Yr.^2);
    logterm=log((R+Yr)./Xr);
    J0=besselj(0,Xr);
    J1=besselj(1,Xr);
    Y0=bessely(0,Xr);
    Y1=bessely(1,Xr);
    eY=exp(-Yr);

    % interpolate look-up tables
    fxy=interp2(lookup.R1.x,lookup.R1.y,lookup.R1.fxy,Xr,Yr,'cubic');
    gxy=interp2(lookup.R1.x,lookup.R1.y,lookup.R1.gxy,Xr,Yr,'cubic');
    hxy=interp2(lookup.R1.x,lookup.R1.y,lookup.R1.hxy,Xr,Yr,'cubic');
    
    % deal with edge cases
    lowX=Xr<min(lookup.R1.x);
    Xlow=0*Xr(lowX)+min(lookup.R1.x);
    fxy(lowX)=interp2(lookup.R1.x,lookup.R1.y,lookup.R1.fxy,Xlow,Yr(lowX),'cubic');
    gxy(lowX)=interp2(lookup.R1.x,lookup.R1.y,lookup.R1.gxy,Xlow,Yr(lowX),'cubic');
    hxy(lowX)=interp2(lookup.R1.x,lookup.R1.y,lookup.R1.hxy,Xlow,Yr(lowX),'cubic');
    
    lowY=Yr<min(lookup.R1.y);
    Ylow=0*Yr(lowY)+min(lookup.R1.y);
    fxy(lowY)=interp2(lookup.R1.x,lookup.R1.y,lookup.R1.fxy,Xr(lowY),Ylow,'cubic');
    gxy(lowY)=interp2(lookup.R1.x,lookup.R1.y,lookup.R1.gxy,Xr(lowY),Ylow,'cubic');
    hxy(lowY)=interp2(lookup.R1.x,lookup.R1.y,lookup.R1.hxy,Xr(lowY),Ylow,'cubic');

    % compute FXY and partial derivatives
    FXY(r1)=-eY.*(2*J0.*logterm+pi*Y0+R.*fxy);
    dFdX(r1)=-eY.*(-2*(Yr.*J0./(Xr.*R)+J1.*logterm)-pi*Y1+gxy.*Xr./R);
    dFdY(r1)=eY.*(2*J0.*(logterm-1./R)+pi*Y0-Yr.*fxy./R+R.*hxy);
end

% Region 2
if sum(r2)>0
    % pre-compute some terms
    Xr=X(r2);
    Yr=Y(r2);
    Xrinv=1./Xr;
    Y0=bessely(0,Xr);
    Y1=bessely(1,Xr);
    eY=2*pi*exp(-Yr);
    eYY0=eY.*Y0;
    % interpolate look-up tables
    fxy=interp2(lookup.R2.xinv,lookup.R2.y,lookup.R2.fxy,Xrinv,Yr,'cubic');
    gxy=interp2(lookup.R2.xinv,lookup.R2.y,lookup.R2.gxy,Xrinv,Yr,'cubic');
    hxy=interp2(lookup.R2.xinv,lookup.R2.y,lookup.R2.hxy,Xrinv,Yr,'cubic');
    % compute FXY and partial derivatives
    FXY(r2)=fxy-eYY0;
    dFdX(r2)=gxy+eY.*Y1;
    dFdY(r2)=hxy+eYY0;
end

% Region 3
if sum(r3)>0
    % pre-compute some terms
    Xr=X(r3);
    Yr=Y(r3);
    Yrinv=1./Yr;
    FXY(r3)=interp2(lookup.R3.x,lookup.R3.yinv,lookup.R3.FXY,Xr,Yrinv,'cubic');
    dFdX(r3)=interp2(lookup.R3.x,lookup.R3.yinv,lookup.R3.dFdX,Xr,Yrinv,'cubic');
    dFdY(r3)=interp2(lookup.R3.x,lookup.R3.yinv,lookup.R3.dFdY,Xr,Yrinv,'cubic');
end

% Region 4
if sum(r4)>0
    % pre-compute some terms
    Xr=X(r4);
    Yr=Y(r4);
    Xrinv=1./Xr;
    Yrinv=1./Yr;
    Y0=bessely(0,Xr);
    Y1=bessely(1,Xr);
    eY=2*pi*exp(-Yr);
    eYY0=eY.*Y0;
    % interpolate look-up tables
    fxy=interp2(lookup.R4.xinv,lookup.R4.yinv,lookup.R4.fxy,Xrinv,Yrinv,'cubic');
    gxy=interp2(lookup.R4.xinv,lookup.R4.yinv,lookup.R4.gxy,Xrinv,Yrinv,'cubic');
    hxy=interp2(lookup.R4.xinv,lookup.R4.yinv,lookup.R4.hxy,Xrinv,Yrinv,'cubic');
    % compute FXY and partial derivatives
    FXY(r4)=fxy-eYY0;
    dFdX(r4)=gxy+eY.*Y1;
    dFdY(r4)=hxy+eYY0;
end

if sum(X0_1)>0
    Yr=Y(X0_1);
    fXY=interp1(lookup.R0_1.y,lookup.R0_1.fXY,Yr,'pchip');
   
    FXY(X0_1)=fXY-2*exp(-Yr).*log(Yr);
    dFdX(X0_1)=0;
    dFdY(X0_1)=-2./Yr-FXY(X0_1);
end

if sum(X0_2)>0
    Yr=Y(X0_2);
    Yrinv=1./Yr;
    FXYi=interp1(lookup.R0_2.yinv,lookup.R0_2.FXY,Yrinv,'pchip');

    FXY(X0_2)=FXYi;
    dFdX(X0_2)=0;
    dFdY(X0_2)=-2./Yr-FXY(X0_2);
end