function zp = NP_unc(t,X)

global A B C a b c d e rho N

% ensures we have a column vector by starting with a vector of 4 rows and
% 1 column
Np = N+1;
x = X(1:Np);
y = X(Np+1:2*Np);
for lp=1:Np
    xx(lp) = a*C(lp) - x'*(b*A(:,:,lp)+rho*B(:,:,lp))*y-e*x(lp); %nutrients
    yy(lp) = x'*(c*b*A(:,:,lp)+c*rho*B(:,:,lp))*y-d*y(lp); %algae
    %xx(lp) = y'*(a*A(:,:,lp)+rho*B(:,:,lp))*x-d*x(lp); % lp=1 2.5000e-04 %old notation
    %yy(lp) = b*C(lp) - y'*(a*c*A(:,:,lp)+rho*c*B(:,:,lp))*x-e*y(lp); %0 % oldnotation
end
%
zp = [xx,yy]';
return