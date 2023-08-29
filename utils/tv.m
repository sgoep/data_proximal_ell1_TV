function [ubar, err, n] = tv(x0, g, K, KT, N, L, lambda, nmax, f, flag)


% Chambolle-Pock algorithm for solving the convex minimization problem
%       min_x 1/2*|| Ax - g||_2^2 + beta*||x||_TV

theta = 1;
sigma = 1/L;
tau   = 1/L;

grad_scale = 10;

u   = x0;
q_x   = zeros(size(x0));
q_y   = zeros(size(x0));

if flag == 1
    p = zeros(size(g));
else
    p = K(zeros(size(g)));
end

ubar = u;

n = 1;
err = zeros(1, nmax);

while n <= nmax
    p = gmultiply(gadd(p, gmultiply(sigma, gsubtract(K(ubar), g))),1/(1+sigma));
    
    if lambda ~= 0
        [ubar_x, ubar_y] = grad(ubar);
        q_x = lambda*(q_x + grad_scale*sigma*ubar_x)./max(lambda, abs(q_x + grad_scale*sigma*ubar_x));
        q_y = lambda*(q_y + grad_scale*sigma*ubar_y)./max(lambda, abs(q_y + grad_scale*sigma*ubar_y));
        v   = u - tau*(KT(p) - grad_scale*div(q_x, q_y));
    else
        v   = u - tau*(KT(p));
    end
% 
    uiter = max(0, v);
    ubar  = uiter + theta*(uiter - u);

    u = uiter;    
    err(n) = sum(abs(ubar(:) - f(:)).^2)/sum(abs(f(:)).^2);
    fprintf('TV regularization. Iteration %d / %d, with relative l2-error %d  \n', n, nmax, err(n))
    n = n + 1;
    
end

end