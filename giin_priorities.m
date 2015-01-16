function [ Pstructure, diffused ] = giin_priorities( vertices, Pstructure, G, gparam )
%GIIN_PRIORITIES Compute the priority of a list of vertices.
%   Receive a list of vertices, update the priority signal.

% tstart = tic;

% Estimate the maximum eigenvalue of the Lapalacian (for filters).
% Done only once as it does not change much when adding new vertices.
% G = gsp_estimate_lmax(G);

% Heat kernel.
Hk = gsp_design_heat(G, gparam.priority.heat_scale);
% figure(1000)
% gsp_plot_filter(G, Hk);

% Wavelet kernel (not a filterbank). Which frequency ?
% Wk = gsp_design_mexican_hat(G, 1);
% gsp_plot_filter(G, Wk);

% Create the Kronecker deltas.
% deltas = eye(G.N);
% deltas = deltas(:,vertices);
deltas = sparse(vertices, 1:length(vertices), ones(size(vertices)), G.N, length(vertices));
% gsp_plot_signal(G, deltas(:,1));

% Graph filtering.
param.cheb_order = gparam.priority.cheb_order;
diffused = gsp_filter_analysis(G, Hk, deltas, param);
% f1_mat = gsp_vec2mat(f1_out, Nf);
% f1_img = reshape(filtered, height, width);
% imshow(f1_img);
% imagesc(f1_img);

% Update priority signal. Not normalized in [0,1].
switch(gparam.priority.type)
    
    case 'threshold'
        Pstructure(vertices) = sum(diffused > gparam.priority.threshold, 1);% / G.N^2;
        % bin = diffused > max(diffused(:)) / 10;

    % Sparsity feature.
    case 'sparsity'
        % sigma(i) = ||T_ig||_1 / || T_ig||_2 = C / || T_ig||_2
        % ||T_ig||_1 = sum(abs(diffused), 1) = 1, i.e. energy is conserved
        % ||T_ig||_2 = sum(diffused.^2, 1)
        Pstructure(vertices) = sum(diffused.^2, 1);
        
    otherwise
        error('Unknown priority type.');
end

% Execution time.
% fprintf('giin_priorities : %f seconds\n', toc(tstart));

end