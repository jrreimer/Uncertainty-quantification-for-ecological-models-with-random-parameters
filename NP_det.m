function zp = NP_det(t,X)

global a b c d e

% ensures we have a column vector by starting with a vector of 2 rows and
% 1 column
zp = [0;0];
%
zp = [a-b*X(1)*X(2)-e*X(1); c*b*X(1)*X(2)-d*X(2)];
%zp = [a*X(1)*X(2)-d*X(1); b-c*a*X(1)*X(2)-e*X(2)];%old notation

%
return