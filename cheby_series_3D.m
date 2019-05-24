function [F, Fx, Fy]=cheby_series_3D(x,y,z,a,amin)

% cut to desired accuracy
if nargin>4
    good=abs(a)>=amin;
    N1=find(sum(sum(good,2),3)>0,1,'last');
    N2=find(sum(sum(good,1),3)>0,1,'last');
    N3=find(sum(sum(good,2),1)>0,1,'last');
    a=a(1:N1,:,:);
    a=a(:,1:N2,:);
    a=a(:,:,1:N3);
    a(abs(a)<amin)=0;
end

% sum coefficients over z
n = (0:(size(a,3)-1)).';
Tz = permute(repmat(cos(n*acos(z)),1,size(a,1),size(a,2)),[2 3 1]);
b = sum(a.*Tz,3);

% convert to monomial coefficients
c = cheb2monomial(b);

% initiate outputs
F  = zeros(size(x));
Fx = zeros(size(x));
Fy = zeros(size(x));

% calculate Chebyshev series
X = 1 + zeros(size(x));
Xm = 0;
for n1 = 0:(size(c,1)-1)
    Y = 1 + zeros(size(y));
    Ym = 0;
    for n2 = 0:(size(c,2)-1)
        F  = F + c(n1+1,n2+1).*X.*Y;
        Fx = Fx + n1*c(n1+1,n2+1).*Xm.*Y;
        Fy = Fy + n2*c(n1+1,n2+1).*X.*Ym;
        Ym = Y;
        Y = Y.*y;
    end
    Xm = X;
    X = X.*x;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function c=cheb2monomial(b)
        
        % check size of input
        s=size(b);
        Nx=s(1)-1;
        Ny=s(2)-1;
        
        % define monomial coefficients
        c=zeros(Nx+1,Ny+1);
        for p=0:Nx
            Pp=chebpoly(p);
            for q=0:Ny
                Pq=chebpoly(q);
                for i=0:p
                    for j=0:q
                        c(i+1,j+1)=c(i+1,j+1)+b(p+1,q+1)*Pp(p-i+1)*Pq(q-j+1);
                    end
                end
            end
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function p=chebpoly(n)
        
        p=zeros(1,n+1);
        
        if n==0
            p=1;
        else
            ind=1;
            for k=0:floor(n/2)
                p(ind)=(n/2)*((-1)^k)*(factorial(n-k-1)/(factorial(k)*factorial(n-2*k)))*(2^(n-2*k));
                ind=ind+2;
            end
        end
    end

end
