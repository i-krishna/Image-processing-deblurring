% Calculates the matrix norm as I want it i.e. norm^2 is equal to sum of
% squares of all values in the matrix.

function [norm_val] = matrix_norm(X,n);

norm_val = norm(X(:),n);
