function [ gparam ] = giin_default_parameters( )
%GIIN_DEFAULT_PARAMETERS Default algorithm parameters
%  Suitable defaults for most use cases.

gparam.graph.psize = 7; % Patch size.
gparam.graph.knn = 10; % Patch graph minimum number of connections (KNN).
gparam.graph.sigma = 1e-1; % Variance of the distance kernel. We want the graph weights to be well spread.
gparam.graph.loc = 0.001; % Importance of local information. (default 0.001, 0.1)

gparam.connect.max_unknown_pixels = gparam.graph.psize; % Maximum number of unknown pixels to connect a patch.

gparam.priority.threshold = 1e-3; % Threshold when creating priority from diffused energy.
gparam.priority.heat_scale = 500; % Depends on sigma. 1000 for lena
gparam.priority.cheb_order = 30; % Order of the Chebyshev approximation (number of hopes).

gparam.inpainting.retrieve = 'copy'; % Average connected patches or copy the strongest.
gparam.inpainting.compose = 'overwrite'; % Keep known pixels or overwrite everything.

gparam.optim.prior = 'thikonov'; % Global optimization constraint : thikonov or tv.
gparam.optim.maxit = 100; % Maximum number of iterations.
gparam.optim.sigma = 0; % Noise level.

end