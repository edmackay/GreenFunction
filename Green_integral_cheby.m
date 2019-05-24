function [G,dGdR,dGdz,dGdzeta]=Green_integral_cheby(K,k0,R,z,zeta,h,lookup)

% vertical terms
s=sign(z-zeta);
v1=abs(z-zeta);
v2=z+zeta+2*h;
v3=abs(z+zeta);
v4=z-zeta+2*h;
v5=zeta-z+2*h;
v6=z+zeta+4*h;

% distance terms
R2=R.^2;
r1inv=(R2 + v1.^2).^(-0.5);
r2inv=(R2 + v2.^2).^(-0.5);
r3inv=(R2 + v3.^2).^(-0.5);
r4inv=(R2 + v4.^2).^(-0.5);
r5inv=(R2 + v5.^2).^(-0.5);
r6inv=(R2 + v6.^2).^(-0.5);
r1inv3=r1inv.^3;
r2inv3=r2inv.^3;
r3inv3=r3inv.^3;
r4inv3=r4inv.^3;
r5inv3=r5inv.^3;
r6inv3=r6inv.^3;

% non-dimensional variables
A=R/h;
B1=v1/h;
B2=v2/h;
H=K*h;

% separate into regions depending on value of B2
ind_low=B2<1;
ind_high=~ind_low;

% map variables onto [-1,1]
x = 2*A-1;
y1 = 2*B1-1;
y2_low = 2*B2(ind_low)-1;
y2_high = 2*B2(ind_high)-3;

% Compute infinite depth part when B2>1
if K==0
    FXV=0;
    dFdX=0;
    dFdV=0;
elseif K==Inf
    FXV=-2*r3inv(ind_high);
    dFdX=2*R(ind_high).*r3inv3(ind_high);
    dFdV=2*v3(ind_high).*r3inv3(ind_high);
else
    X=K*R(ind_high);
    V3=K*v3(ind_high);
    [FXV,dFdX,dFdV]=FXY_interp(X,V3,lookup);
end

h2=h^2;
H2=H^2;

% calculate finite depth components
if H<=1
    % compute integrals L1 and L2
    [zL1,aL1]=cheby_coeff(H,'L1');
    [G1, G1dA, G1dB]=cheby_series_3D(x,y1,zL1,aL1);
    
    [zL2,aL2]=cheby_coeff(H,'L2');
    L2=cheby_series_1D(zL2,aL2);
    if H>0
        L2=L2-log(H)/2;
    end
    
    % B1 component
    G1 = G1 + L2;
    
    % account for change of variable in differentiation
    G1dA = 2*G1dA;
    G1dB = 2*G1dB;
    
    % compute components for high and low values of B2
    if sum(ind_low)>0
        % compute integrals L1 and L2
        [G2_low, G2dA_low, G2dB_low]=cheby_series_3D(x(ind_low),y2_low,zL1,aL1);
        
        % G2 low component
        G2_low = G2_low + L2;
        G2dA_low = 2*G2dA_low;
        G2dB_low = 2*G2dB_low;
    else
        G2_low=[];
        G2dA_low=[];
        G2dB_low=[];
    end
    if sum(ind_high)>0
        % compute integrals M1 and M2
        [zM1,aM1]=cheby_coeff(H,'M1');
        [G2_high, G2dA_high, G2dB_high]=cheby_series_3D(x(ind_high),y2_high,zM1,aM1);
        
        [zM2,aM2]=cheby_coeff(H,'M2');
        M2=cheby_series_1D(zM2,aM2);
        if H>0
            M2=M2-log(H)/2;
        end
        
        % account for change of variable in differentiation
        if K==Inf
            G2_high = G2_high + M2 + h*FXV;
            G2dA_high = 2*G2dA_high + h2*dFdX;
            G2dB_high = 2*G2dB_high - h2*dFdV;
        else
            G2_high = G2_high + M2 + H*FXV;
            G2dA_high = 2*G2dA_high + H2*dFdX;
            G2dB_high = 2*G2dB_high - H2*dFdV;
        end
    else
        G2_high=[];
        G2dA_high=[];
        G2dB_high=[];
    end
    
else % H>1
    % compute integral L3
    [zL3,aL3]=cheby_coeff(H,'L3');
    [G1, G1dA, G1dB]=cheby_series_3D(x,y1,zL3,aL3);
    
    % subtract Rankine terms
    G1 = G1 - 2*h*(r4inv + r5inv);
    G1dA = 2*G1dA + 2*h2*R.*(r4inv3 + r5inv3);
    G1dB = 2*G1dB + 2*h2*s.*(v4.*r4inv3 - v5.*r5inv3);
    
    % compute components for high and low values of B2
    if sum(ind_low)>0
        % compute integral L3
        [G2_low, G2dA_low, G2dB_low]=cheby_series_3D(x(ind_low),y2_low,zL3,aL3);
        
        % subtract Rankine terms
        G2_low = G2_low - 2*h*(r3inv(ind_low) + r6inv(ind_low));
        G2dA_low = 2*G2dA_low + 2*h2*R(ind_low).*(r3inv3(ind_low) + r6inv3(ind_low));
        G2dB_low = 2*G2dB_low + 2*h2*(-v3(ind_low).*r3inv3(ind_low) + v6(ind_low).*r6inv3(ind_low));
    else
        G2_low=[];
        G2dA_low=[];
        G2dB_low=[];
    end
    if sum(ind_high)>0
        % compute integrals M1 and M2
        [zM3,aM3]=cheby_coeff(H,'M3');
        [G2_high, G2dA_high, G2dB_high]=cheby_series_3D(x(ind_high),y2_high,zM3,aM3);
        
        % account for change of variable in differentiation
        if K==Inf
            G2_high = G2_high + h*FXV - 2*h*r6inv(ind_high);
            G2dA_high = 2*G2dA_high + h2*dFdX + 2*h2*R(ind_high).*r6inv3(ind_high);
            G2dB_high = 2*G2dB_high - h2*dFdV + 2*h2*v6(ind_high).*r6inv3(ind_high);
        else
            G2_high = G2_high + H*FXV - 2*h*r6inv(ind_high);
            G2dA_high = 2*G2dA_high + H2*dFdX + 2*h2*R(ind_high).*r6inv3(ind_high);
            G2dB_high = 2*G2dB_high - H2*dFdV + 2*h2*v6(ind_high).*r6inv3(ind_high);
        end
    else
        G2_high=[];
        G2dA_high=[];
        G2dB_high=[];
    end
end

% combine components for high and low values of B2
G2=zeros(size(x));
G2dA=zeros(size(x));
G2dB=zeros(size(x));
G2(ind_low)=G2_low;
G2(ind_high)=G2_high;
G2dA(ind_low)=G2dA_low;
G2dA(ind_high)=G2dA_high;
G2dB(ind_low)=G2dB_low;
G2dB(ind_high)=G2dB_high;

% Green function and derivatives
G = r1inv + r2inv + r3inv + r4inv + r5inv + r6inv + (G1 + G2)/h;
dGdR = -R.*(r1inv3 + r2inv3 + r3inv3 + r4inv3 + r5inv3 + r6inv3) + (G1dA + G2dA)/h2;
dGdz    = -( s.*v1.*r1inv3 + v2.*r2inv3 - v3.*r3inv3 + v4.*r4inv3 - v5.*r5inv3 + v6.*r6inv3) + ( s.*G1dB + G2dB)/h2;
dGdzeta = -(-s.*v1.*r1inv3 + v2.*r2inv3 - v3.*r3inv3 - v4.*r4inv3 + v5.*r5inv3 + v6.*r6inv3) + (-s.*G1dB + G2dB)/h2;

if K>0 && K<inf
    % Imaginary part
    J0=besselj(0,k0*R);
    J1=besselj(1,k0*R);
    C0=k0^2/((k0^2-K^2)*h+K);
    aM1=2*pi*1i*C0/((1+exp(-2*k0*h))^2);
    ekv3=exp(-k0*v3);
    ekv4=exp(-k0*v4);
    ekv5=exp(-k0*v5);
    ekv6=exp(-k0*v6);
    
    Im_G=-aM1*J0.*(ekv3 + ekv4 + ekv5 + ekv6);
    Im_dGdR=aM1*k0*J1.*(ekv3 + ekv4 + ekv5 + ekv6);
    Im_dGdz=-aM1*k0*J0.*(ekv3 - ekv4 + ekv5 - ekv6);
    Im_dGdzeta=-aM1*k0*J0.*(ekv3 + ekv4 - ekv5 - ekv6);
    
    % Add in imaginary component
    G=G+Im_G;
    dGdR=dGdR+Im_dGdR;
    dGdz=dGdz+Im_dGdz;
    dGdzeta=dGdzeta+Im_dGdzeta;
end