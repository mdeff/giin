function [ sol ] = giin_global( G, img, imgstart, gparam )
%GIIN_GLOBAL Global stage by convex optimization.
%   We inpaint again the image using the created non-local graph.

tstart = tic;

init_unlocbox();
verbose = 1;

% Observed signal (image).
M = reshape(img(:,:,1)>=0, [], 1);
M = repmat(M,1,size(img,3));
y = M .* reshape(img, [], size(img,3));
imgstart = reshape(imgstart, [], size(img,3));

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

% Optimization.
if ~gparam.optim.sigma
   fdata.prox = @(x,T) x - M.*x + y;
end

% Prior.
param_prior.verbose = verbose-1;
switch(gparam.optim.prior)
    
    % Thikonov prior (gradient is faster than proximal operator).
    case 'thikonov'
        %fprior.prox = @(x,T) gsp_prox_tik(x,T,G,param_prior);
        fprior.eval = @(x) sum(gsp_norm_tik(G,x));
        fprior.grad = @(x) 2*G.L*x;
    
    % TV prior.
    case 'tv'
        G = gsp_adj2vec(G);
        G = gsp_estimate_lmax(G);
        fprior.prox = @(x,T) gsp_prox_tv(x,T,G,param_prior);
        fprior.eval = @(x) sum(gsp_norm_tv(G,x));
        
    otherwise
        error('Unknown prior.');
end

% Solve the convex optimization problem.
param_solver.verbose = verbose;
param_solver.tol = 1e-12;
param_solver.maxit = gparam.optim.maxit;
if strcmp(gparam.optim.prior,'thikonov')
    param_solver.gamma = 0.5/G.lmax;
    [sol, info] = forward_backward(imgstart,fdata,fprior,param_solver);
else
    [sol, info] = douglas_rachford(imgstart,fdata,fprior,param_solver);
end

% Execution time.
fprintf('Global optimization : %f seconds (%d iterations)\n', toc(tstart), info.iter);

end