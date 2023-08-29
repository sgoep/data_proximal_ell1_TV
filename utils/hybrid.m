function [ubar, err] = hybrid(x0, g, A, AT, Psi, PsiT, N, alpha, beta, L, nmax, f)

% Chambolle-Pock algorithm for solving the convex minimization problem
%       min_x 1/2*|| Ax - g||_2^2 + alpha*||Psi(x)||_1 + beta*||x||_TV

theta = 1;
sigma = 1/L;
tau   = 1/L;

grad_scale = 10; % rescaling of gradient

u   = x0;

q_x   = zeros(N);
q_y   = zeros(N);

p = zeros(size(g));
ubar = u;

n = 1;
err = zeros(1, nmax);
while n <= nmax
    p = (p + sigma*(A(PsiT(ubar)) - g))/(1+sigma);

    [ubar_x, ubar_y] = grad(PsiT(ubar));
    if beta ~= 0
        q_x = beta*(q_x + grad_scale*sigma*ubar_x)./max(beta, abs(q_x + grad_scale*sigma*ubar_x));
        q_y = beta*(q_y + grad_scale*sigma*ubar_y)./max(beta, abs(q_y + grad_scale*sigma*ubar_y));
        uiter = soft_thresholding(gsubtract(u, gmultiply(tau, gsubtract(Psi(AT(p)), Psi(grad_scale*div(q_x, q_y))))), alpha);
    else
        uiter = soft_thresholding(gsubtract(u, gmultiply(tau, Psi(AT(p)))), alpha);
    end

    
    ubar = gadd(uiter, gmultiply(theta, gsubtract(uiter, u)));

    u = uiter;   
    tmp = PsiT(ubar);
    err(n) = sum(sum(abs(tmp - f).^2))/sum(abs(f(:)).^2);
    fprintf('Hybrid regularization. Iteration %d / %d, with relative l2-error %d  \n', n, nmax, err(n))
    n = n + 1;
    
end
end




