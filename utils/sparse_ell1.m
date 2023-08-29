function [ubar, err] = sparse_ell1(x0, g, h, A, AT, Psi, PsiT, N, L, lambda, nmax, f, flag)

% Chambolle-Pock algorithm for solving the convex minimization problem
%       min_x 1/2*|| A PsiT (c) - g||_2^2 + alpha*||c||_1

theta = 1;
sigma = 1/L;
tau   = 1/L;

u = x0;
p = zeros(size(g));

q = zeros(size(h));


ubar = u;

n = 1;
err = zeros(1, nmax);

while n <= nmax
    p = gmultiply(gadd(p, gmultiply(sigma, gsubtract(A(PsiT(ubar)), g))),1/(1+sigma));
    if flag == 1 
        q = gmultiply(gadd(q, gmultiply(sigma, gsubtract(A(PsiT(ubar)), h))),1/(1+sigma));
        v = gsubtract(u, gmultiply(tau, Psi(AT(p + q))));
    else 
        v = gsubtract(u, gmultiply(tau, Psi(AT(p))));
    end

    uiter = soft_thresholding(v, lambda);
    ubar = gadd(uiter, gmultiply(theta, gsubtract(uiter, u)));

    u = uiter;

    err(n) = sum(abs(reshape(PsiT(ubar), (N)^2,1) - f(:)).^2)/sum(abs(f(:)).^2);
    fprintf('Sparse ell1 regularization. Iteration %d / %d, with relative l2-error %d  \n', n, nmax, err(n))
    n = n + 1;
    
end
end
