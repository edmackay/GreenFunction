function [G, dGdR, dGdz, dGdzeta]=Green_finite_depth(K,k,R,z,zeta,h,Xa,Y1,Y2,Flookup)

% computes the finite depth Green function and its partial derivatives using the method described by Mackay (2018)
% Inputs:
%   K = infinite depth wavenumber
%   k = finite depth wavenumber
%   R = horizontal separation of source and field point
%   z = vertical coordinate of source point (z=0 at free surface and is positive upwards)
%   zeta = vertical coordinate of field point
%   h = water depth
%   Xa = N x 6 array of powers of x=2*A-1, where column n is x^(n-1)
%   Y1 = N x 8 array of powers of y=B1-1, where column n is y^(n-1)
%   Y2 = N x 8 array of powers of y=B2-1, where column n is y^(n-1)
%   Flookup is the lookup table for F(X,Y)

% ensure inputs are column vectors
R=R(:);
z=z(:);
zeta=zeta(:);

% initialise output variables
G=zeros(size(R));
dGdR=zeros(size(R));
dGdz=zeros(size(R));
dGdzeta=zeros(size(R));

% divide domain into various regions
region1 = R/h<=1;
region2 = R/h>1;

% calculate Green function using series or integral method
if K==0
    if sum(region1)>0
        [G(region1),dGdR(region1),dGdz(region1),dGdzeta(region1)]=Green_integral_zero(R(region1),z(region1),zeta(region1),h,Xa(region1,:),Y1(region1,:),Y2(region1,:));
    end
    if sum(region2)>0
        [G(region2),dGdR(region2),dGdz(region2),dGdzeta(region2)]=Green_series_zero(R(region2),z(region2),zeta(region2),h);
    end
elseif K==inf
    if sum(region1)>0
        [G(region1),dGdR(region1),dGdz(region1),dGdzeta(region1)]=Green_integral_infty(R(region1),z(region1),zeta(region1),h,Xa(region1,:),Y1(region1,:),Y2(region1,:));
    end
    if sum(region2)>0
        [G(region2),dGdR(region2),dGdz(region2),dGdzeta(region2)]=Green_series_infty(R(region2),z(region2),zeta(region2),h);
    end
else
    if sum(region1)>0
        [G(region1),dGdR(region1),dGdz(region1),dGdzeta(region1)]=Green_integral(K,k,R(region1),z(region1),zeta(region1),h,Xa(region1,:),Y1(region1,:),Y2(region1,:),Flookup);
    end
    if sum(region2)>0
        [G(region2),dGdR(region2),dGdz(region2),dGdzeta(region2)]=Green_series(K,k,R(region2),z(region2),zeta(region2),h);
    end
end

end