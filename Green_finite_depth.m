function [G, dGdR, dGdz, dGdzeta]=Green_finite_depth(K,k0,R,z,zeta,h,Flookup)

% computes the finite depth Green function and its partial derivatives using the method described by Mackay (2018)
% Inputs:
%   K = infinite depth wavenumber
%   k0 = finite depth wavenumber
%   R = horizontal separation of source and field point
%   z = vertical coordinate of source point (z=0 at free surface and is positive upwards)
%   zeta = vertical coordinate of field point
%   h = water depth
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
if sum(region1)>0
    [G(region1),dGdR(region1),dGdz(region1),dGdzeta(region1)]=Green_integral_cheby(K,k0,R(region1),z(region1),zeta(region1),h,Flookup);
end
if sum(region2)>0
    if K==0
        [G(region2),dGdR(region2),dGdz(region2),dGdzeta(region2)]=Green_series_zero(R(region2),z(region2),zeta(region2),h);
    elseif K==inf
        [G(region2),dGdR(region2),dGdz(region2),dGdzeta(region2)]=Green_series_infty(R(region2),z(region2),zeta(region2),h);
    else
        [G(region2),dGdR(region2),dGdz(region2),dGdzeta(region2)]=Green_series(K,k0,R(region2),z(region2),zeta(region2),h);
    end
end