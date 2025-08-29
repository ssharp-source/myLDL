# myldl.m - LDL^T decomposition
# Copyright (C) 2025 Steven Sharp
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#
# -*- texinfo -*-
# @deftypefn {} {[@var{L}, @var{D}] =} myldl (@var{A})
#
# Compute the LDL^T decomposition of a symmetric positive definite (or
# semidefinite) matrix A.
#
# Decomposes A into L*D*L' where:
#   - L is unit lower triangular (1s on diagonal)
#   - D is diagonal (simple version)
#
# Example:
# @example
# A = [4 2 2; 2 2 0; 2 0 5];
# [L, D] = myldl(A);
# A_recon = L*D*L';
# norm(A - A_recon)   % should be close to 0
# @end example
#
# @end deftypefn

function [L, D] = myldl(A)
  if nargin != 1
    print_usage();
  endif

  [n, m] = size(A);
  if n != m
    error("Matrix must be square.");
  endif

  % Ensure symmetry
  if norm(A - A', "fro") > 1e-12
    warning("Matrix is not exactly symmetric. Using (A+A')/2.");
    A = (A + A')/2;
  endif

  L = eye(n);
  D = zeros(n);

  for j = 1:n
    % Compute D(j,j)
    sum_val = 0;
    for k = 1:j-1
      sum_val = sum_val + (L(j,k)^2) * D(k,k);
    endfor
    D(j,j) = A(j,j) - sum_val;

    % Compute L(i,j) for i>j
    for i = j+1:n
      sum_val = 0;
      for k = 1:j-1
        sum_val = sum_val + L(i,k)*D(k,k)*L(j,k);
      endfor
      L(i,j) = (A(i,j) - sum_val) / D(j,j);
    endfor
  endfor
endfunction