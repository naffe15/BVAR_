function G = rand_inverse_wishart(m, v, H_inv_upper_chol)

% function G = rand_inverse_wishart(m, v, H_inv_upper_chol)
% rand_inverse_wishart  Pseudo random matrices drawn from an
% inverse Wishart distribution
% G = rand_inverse_wishart(m, v, H_inv_upper_chol)
% Returns an m-by-m matrix drawn from an inverse-Wishart distribution.
%
% INPUTS:
%     m:          dimension of G and H_inv_upper_chol.
%     v:          degrees of freedom, greater or equal than m.
%     H_inv_chol: upper cholesky decomposition of the inverse of the
%                 matrix parameter.
%                 The upper cholesky of the inverse is requested here
%                 in order to avoid to recompute it at every random draw.
%                 H_inv_upper_chol = chol(inv(H))
% OUTPUTS:
%     G:          G ~ IW(m, v, H) where H = inv(H_inv_upper_chol'*H_inv_upper_chol)
%                 or, equivalently, using the correspondence between Wishart and
%                 inverse-Wishart: inv(G) ~ W(m, v, S) where 
%                 S = H_inv_upper_chol'*H_inv_upper_chol = inv(H)
%  
% SPECIAL REQUIREMENT
%     none
%    

% Copyright (C) 2003-2009 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

X = randn(v, m) * H_inv_upper_chol; 


% At this point, X'*X is Wishart distributed
% G = inv(X'*X);

% Rather compute inv(X'*X) using the SVD
[U,S,V] = svd(X, 0);
SSi = 1 ./ (diag(S) .^ 2);
G = (V .* repmat(SSi', m, 1)) * V';