function k=wavenumber_real(K,h)

% solves the dispersion relation K=k*tanh(k*h)
%
% K and h are vectors of wave periods and water depths
% K and h must be scalars,vectors or matricies of the same dimensions
% h=0 implies deep water

if h==0
    k=K;
else
    % initial guess
    Kh=K*h;
    kh=Kh;
    
    % iterative improvements
    error_tolerance=ones(size(h));
    while max(error_tolerance)>1e-7
        % values of function and derivative
        tanh_kh=tanh(kh);
        f1=(kh.*tanh_kh) - Kh;
        f2=tanh_kh+kh.*(sech(kh).^2);
        
        % change in estimate
        dkh= f1./f2;
        
        % next iteration
        kh=kh-dkh;
        
        % check if error tolerance is exceeded
        error_tolerance=abs(dkh./kh);
    end
    k=kh./h;
end