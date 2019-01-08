function k=wavenumber_imag(K,h,n)

% finds the positive real roots of the K=-k*tan(k*h)
% K is the infinite-depth wavenumber
% h is water depth
% n is the number of solutions required

Kh=K*h;

kh=zeros(n,1);
for i=1:n
    % first guess 
    xlow=-atan(Kh./((i-1/2)*pi));
    xhigh=-atan(Kh./(i*pi));
    xmid=0.5*(xlow+xhigh);
    x0=i*pi+xmid;
    % iterate using newton-raphson formula
    dif=1;
    while dif>1e-7
        fx0=Kh/x0+tan(x0);
        dfx0=-Kh/x0^2+1+tan(x0)^2;
        x1=x0-fx0/dfx0;
        dif=abs(x0-x1);
        x0=x1;
    end
    kh(i)=x0;
end
k=kh/h;

