% cshift2.m
%
% Circular shift of a matrix.
% Usuage : y = cshift(x, ti, tj)
% x - input vector
% ti - number of spots to shift columns
% tj - number of spots to shift rows
% Direction is determined by ti,tj being positive/negative.
%
% Written by : Justin Romberg
% Created : 3/15/99

function y = cshift(x, ti, tj)

if (ti < 0)
  diri = 'l';
  ti = abs(ti);
else
  diri = 'r';
end
if (tj < 0)
  dirj = 'l';
  tj = abs(tj);
else
  dirj = 'r';
end

N = length(x);
ti = mod(ti, N);
tj = mod(tj, N);
if (diri == 'r')
  y = [x(:,N-ti+1:N) x(:,1:N-ti)];
else
  y = [x(:,1+ti:N) x(:,1:ti)];
end
if (dirj == 'r')
  y = [y(N-tj+1:N,:) ; y(1:N-tj,:)];
else
  y = [y(1+tj:N,:) ; y(1:tj,:)];
end

