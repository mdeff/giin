function giin_plot_priorities( vertices, G, gparam, savefig )
%GIIN_PLOT_PRIORITIES Visualize how priority is constructed.
%   Given a list of vertices, show how their priority is constructed.

Nplots = sqrt(length(vertices));
width = max(G.coords(:,1));
height = max(G.coords(:,2));

Pstructure = nan(G.N, 1);
[Pstructure, diffused] = giin_priorities(vertices, Pstructure, G, gparam);

for n = 1:length(vertices)
    vertex = vertices(n);
    
    % Heat diffusion through the graph connections.
    fig1 = figure(100);
    subplot(floor(Nplots), ceil(Nplots), find(vertex==vertices));
    clear param;
    param.vertex_highlight = vertex;
    gsp_plot_signal(G, diffused(:,n), param);
    title(['Vertex ',num2str(vertex)]);

    % Binary edge image.
    fig2 = figure(101);
    subplot(floor(Nplots), ceil(Nplots), find(vertex==vertices));
    bin = diffused(:,n) > gparam.priority.threshold;
    bin = reshape(bin, height, width);
    imshow(bin);
    title(['Vertex ',num2str(vertex)]);
    xlabel(['Priority ',num2str(Pstructure(vertex))]);

    % Hough transform.
    fig3 = figure(103);
    subplot(floor(Nplots), ceil(Nplots), find(vertex==vertices));
    H = hough(bin);
    imshow(imadjust(mat2gray(H)));
    title(['Vertex ',num2str(vertex)]);
    colormap(hot);
end

if savefig
    saveas(fig1,[imtype,'_heat',num2str(heat_scale),'_diffusion.png']);
    saveas(fig2,[imtype,'_heat',num2str(heat_scale),'_bin.png']);
    saveas(fig3,[imtype,'_heat',num2str(heat_scale),'_hough.png']);
end

end