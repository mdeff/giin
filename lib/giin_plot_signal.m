function [ ] = giin_plot_signal( G, pixels, show_edges )
%GIIN_PLOT_SIGNAL Plot pixels with missing values in red.

if nargin < 3
    show_edges = false;
end

% fig = figure();

% Any negative value is a missing pixel --> red.
cmap = [1,0,0;gray];
colormap(cmap); % colormap(fig, cmap);
param.climits = [-1/(length(cmap)-1),1];

param.colorbar = 0;
% param.vertex_highlight = connected; % draw by hand in different colors instead
param.show_edges = show_edges; % very slow

gsp_plot_signal(G, pixels, param);

end