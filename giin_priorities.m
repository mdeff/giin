function [ Pstructure, diffused ] = giin_priorities( vertices, Pstructure, G, gparam )
%GIIN_PRIORITIES Compute the priority of a list of vertices.
%   Receive a list of vertices, update the priority signal.

tic;

% Estimate the maximum eigenvalue of the Lapalacian (for filters).
G = gsp_estimate_lmax(G);

% Heat kernel.
Hk = gsp_design_heat(G, gparam.priority.heat_scale);
% gsp_plot_filter(G, Hk);

% Create the Kronecker deltas.
deltas = eye(G.N);
deltas = deltas(:,vertices);

% Graph filtering.
diffused = gsp_filter_analysis(G, Hk, deltas);

% Update priority signal. Normalized in [0,1].
Pstructure(vertices) = sum(diffused > gparam.priority.threshold, 1);% / G.N^2;

% Execution time.
% fprintf('giin_priorities : %f seconds\n', toc);

end