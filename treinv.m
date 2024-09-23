function [Tv] = treinv(X,TRE,FLE)
%TREINV Summary of this function goes her
[K,N] = size(X);
meanX = mean(X')';
Xc = X-meanX*ones(1,N);  % demeaned X
[V, Lambda, U] = svd(Xc);

% RMS distances of markers to their three principal axes:
F1 = (Lambda(2,2)^2 + Lambda(3,3)^2)/N;
F2 = (Lambda(3,3)^2 + Lambda(1,1)^2)/N;
F3 = (Lambda(1,1)^2 + Lambda(2,2)^2)/N;

dx = sqrt((N*TRE^2/FLE^2 - 1)*3*F2*F3/(F2 + F3));
dy = sqrt((N*TRE^2/FLE^2 - 1)*3*F1*F3/(F1 + F3));
dz = sqrt((N*TRE^2/FLE^2 - 1)*3*F2*F1/(F2 + F1));

Tv = [V*[dx;0;0] ,V*[0;dy;0] , V*[0;0;dz]];
% Tv = [dx;dy;dz];

end

