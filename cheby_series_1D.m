function F=cheby_series_1D(z,a,amin)

% cut to desired accuracy
if nargin>2
    a(abs(a)<amin)=[];
end

% compute Chebyshev series
N=length(a);
n=(0:(N-1)).';
F=sum(a.*cos(n*acos(z)));

end