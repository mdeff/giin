function [ G, pixels, Pstructure, Pinformation ] = giin_inpaint( G, pixels, patches, gparam, liveplot )
%GIIN_INPAINT Inpaint a signal and update the corresponding graph.
%   Usage :
%       G = giin_inpaint( G, pixels, patches, gparam, liveplot )
%       [ G, pixels ] = giin_inpaint( G, pixels, patches, gparam, liveplot )
%       [ G, pixels, Pstructure, Pinformation ] = giin_inpaint( G, pixels, patches, gparam, liveplot )
%
%   Input parameters :
%       G        : GSPbox graph structure
%       pixels   : pixel signal
%       patches  : patch signal (could be created)
%       gparam   : parameter structure
%       liveplot : true to plot the process evolution
%
%   Output parameters :
%       G            : retrieved graph
%       pixels       : updated pixel signal
%       Pstructure   : priority signal according to structure
%       Pinformation : priority signal according to information / confidence
%

% Author: Michaël Defferrard
% Date: November 2014

tstart = tic;

% Each unknown pixel has a value of -1e3. A patch with 4 unknown pixels
% will end up with a value of -4. The minimum is -psize^2.
patches(isnan(patches)) = -1e3;
pixels(isnan(pixels)) = -1e3;

unknowns = (patches<0) .* patches;
unknowns = sum(unknowns,2) / 1e3;
Nc = size(pixels,2);

% List of new vertices considered for inpainting.
% news = find(unknowns<0).';
currents = [];
inpainted = [];

% Patches which contain no other information than their position cannot be
% connected in the non-local graph.

fprintf('There is %d incomplete patches :\n', sum(unknowns<0));
fprintf('  %d without any information\n', sum(unknowns==-gparam.graph.psize^2));
% fprintf('  %d considered for inpainting\n', length(news));

% List of fully known patches to which we can connect.
knowns = find(unknowns==0);
if sum(unknowns<0)+length(knowns) ~= G.N
    error('Missing vertices !');
end

% Structure priority.
Pstructure = nan(G.N, 1);

% Information priority. First column is pixels priority, second is patches.
% We don't want it to be zero to preserve sign when multiplied with Pstructure.
Pinformation = zeros(G.N,2);
Pinformation(pixels(:,1)>=0,1) = 1;
Pinformation(pixels(:,1)<0, 1) = 1e-8;
Pinformation(:,2) = nan(G.N,1);
patch_pixels = giin_patch_vertices('pixels', gparam.graph.psize, max(G.coords(:,2)));
for patch = find(unknowns<0).'
    Pinformation(patch,2) = mean(Pinformation(patch+patch_pixels,1));
end

% currents  : vertices to be inpainted
% news      : vertices to be connected, i.e. newly considered pixels
% inpainted : vertices already inpainted, to avoid an infinite loop

% Until no more vertices to inpaint.
first = true;
while ~isempty(currents) || first
    first = false;
    
    % Each unknown pixel has a value of -1e3. A patch with 4 unknown pixels
    % will end up with a value of -4. The minimum is -psize^2.
    unknowns = (patches<0) .* patches;
    unknowns = sum(unknowns,2) / 1000;
    
    news = find(unknowns<0).';

    % We only consider patches with less than some number of missing pixels.
    news = news(unknowns(news)>=-gparam.connect.max_unknown_pixels*Nc);
    % Which are not already connected.
    news = news(~ismember(news, currents));
    % Neither already visited (to prevent infinite loop and reconnections).
    news = news(~ismember(news, inpainted));
    currents = [currents, news]; %#ok<AGROW>

    if any(ismember(currents, inpainted))
        error('A vertex could be visited again !');
    end
    
    % What if some patch has no more unknown pixels ?
    % Do we always inpaint over ?

    % Connect the newly reachable vertices.
    G = giin_connect(G, news, knowns, patches, gparam);

    % Compute their priorities, i.e. update the priority signal.
    Pstructure = giin_priorities(news, Pstructure, G, gparam);

    % TODO: we also need to take into account the data priority, and normalize
    % the two to give them the same weight.

    % Highest priority patch. Negate the value so that it won't be selected
    % again while we keep the information ([0,1] --> [-1,-2]).
    [~,vertex] = max( (sign(Pstructure).*abs(Pstructure).^gparam.priority.p) .* Pinformation(:,2));
    Pstructure(vertex) = -1-Pstructure(vertex);

    if ~ismember(vertex, currents)
        error('This vertex is not in the list of vertices to be inpainted !');
    end
    
    % Update pixels and patches.
    [pixels, patches, Pinformation, ~] = giin_inpaint_patch(vertex, G, pixels, patches, Pinformation, gparam);

    % Remove the currently impainted vertex from the lists.
%     news = news(news~=vertex);
    currents = currents(currents~=vertex);
    inpainted = [inpainted, vertex]; %#ok<AGROW>
    
    % Live plot.
    if liveplot
        figure(10);
        width = max(G.coords(:,1));
        height = max(G.coords(:,2));
        imshow(reshape(pixels,height,width,Nc), 'InitialMagnification',600);
        drawnow;
    end

    fprintf('Inpainted vertices : %d (%d waiting)\n', length(inpainted), length(currents));
end

% Restore priorities.
Pstructure = -1-Pstructure;

% Execution time.
fprintf('Iterative inpainting : %f seconds\n', toc(tstart));

end