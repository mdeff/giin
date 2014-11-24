function [ sol ] = giin_global( G, img, gparam )
%GIIN_GLOBAL Global stage by convex optimization.
%   We inpaint again the image using the created non-local graph.

tstart = tic;

init_unlocbox();
verbose = 1;

% Observed signal (image).
M = reshape(img>=0, [], 1);
y = M .* reshape(img, [], 1);

% Data term.
% fdata.grad = @(x) 2*M.*(M.*x-y);
% fdata.eval = @(x) norm(M.*x-y)^2;
param_b2.verbose = verbose -1;
param_b2.y = y;
param_b2.A = @(x) M.*x;
param_b2.At = @(x) M.*x;
param_b2.tight = 1;
param_b2.epsilon = gparam.optim.sigma*sqrt(sum(M(:)));
fdata.prox = @(x,T) proj_b2(x,T,param_b2);
fdata.eval = @(x) eps;

% Prior.
param_prior.verbose = verbose-1;
switch(gparam.optim.prior)
    
    % Thikonov prior.
    case 'thikonov'
        fprior.prox = @(x,T) gsp_prox_tik(x,T,G,param_prior);
        fprior.eval = @(x) gsp_norm_tik(G,x);
    
    % TV prior.
    case 'tv'
        G = gsp_adj2vec(G);
        G = gsp_estimate_lmax(G);
        fprior.prox = @(x,T) gsp_prox_tv(x,T,G,param_prior);
        fprior.eval = @(x) gsp_norm_tv(G,x);
end

% Solve the convex optimization problem.
param_solver.verbose = verbose;
param_solver.tol = 1e-7;
param_solver.maxit = gparam.optim.maxit;
[sol, info] = douglas_rachford(y,fprior,fdata,param_solver);

% Execution time.
fprintf('Global optimization : %f (%d iterations)\n', toc(tstart), info.iter);

end