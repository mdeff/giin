function [ gparam ] = giin_default_parameters( )
%GIIN_DEFAULT_PARAMETERS Default algorithm parameters
%  Suitable defaults for most use cases.

gparam.graph.psize = 5; % Patch size.
gparam.graph.knn = 10; % Patch graph minimum number of connections (KNN).
gparam.graph.sigma = 0.05; % Variance of the distance kernel. We want the graph weights to be well spread.
gparam.graph.loc = 0.0003; % Importance of local information. (default 0.001, 0.1)
gparam.graph.symetrize_type = 'full'; % Symmetrization applied to the weight matrix (full, average or none).

gparam.connect.max_unknown_pixels = 2*gparam.graph.psize; % Maximum number of unknown pixels to connect a patch.

gparam.priority.type = 'threshold_l1'; % Priority type : 'nthreshold', 'sparsity', 'threshold_l1', 'threshold_l2'.
gparam.priority.thresholda = 0.005; % Threshold when creating priority from diffused energy (nthreshold, threshold_l1, threshold_l2)
gparam.priority.thresholdb = 0.1;
gparam.priority.heat_scale = 50; % Depends on sigma. 1000 for lena
gparam.priority.order = 30; % Order of the Chebyshev approximation (number of hopes). Linear time for gsp_cheby_op.
gparam.priority.p = 0.5; % Balance between structure and information priority. Higher the number, higher the weight of structure.

gparam.inpainting.psize = 3; % Size of the patch being inpainted. Could be smaller than the comparizon patch.
gparam.inpainting.retrieve = 'copy'; % Average connected patches (average) or copy the strongest (copy).
gparam.inpainting.compose = 'mixed'; % Keep known pixels (mixed) or overwrite everything (overwrite).

gparam.optim.prior = 'thikonov'; % Global optimization constraint : thikonov or tv.
gparam.optim.maxit = 1000; % Maximum number of iterations.
gparam.optim.sigma = 0; % Noise level.

end