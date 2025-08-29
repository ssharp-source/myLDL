## Unit tests for myldl

%!test
%! A = [4 2 2; 2 2 0; 2 0 5];
%! [L, D] = myldl(A);
%! A_recon = L*D*L';
%! assert(norm(A - A_recon, "fro") < 1e-10);

%!test
%! A = [2 0.1; 0.1 2];
%! [L, D] = myldl(A);
%! A_recon = L*D*L';
%! assert(norm(A - A_recon, "fro") < 1e-10);

%!test
%! % Identity matrix should give L = I, D = I
%! A = eye(3);
%! [L, D] = myldl(A);
%! assert(norm(L - eye(3)) < 1e-12);
%! assert(norm(D - eye(3)) < 1e-12);
