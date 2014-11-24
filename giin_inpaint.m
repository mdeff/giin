function [ pixels, patches, Pinformation, impacted_patches ] = giin_inpaint( vertex, G, pixels, patches, Pinformation, gparam )
%GIIN_INPAINT Inpaint the specified patch and returns the impacted vertices.
%   Detailed explanation goes here

tic;

%% Impacted patches, i.e. patches that should be investigated next
% List of impacted patches, i.e. patches where an unknown pixel is now
% known. If the value is solely modified, it's not impacted (otherwise we
% would never be able to stop).

height = max(G.coords(:,2));

% Pixel vertices contained in a patch.
patch_pixels = giin_patch_vertices('pixels', gparam.graph.psize, height);

impacted_pixels = patch_pixels .* (patches(vertex,1:end-2)<0);
[~,~,impacted_pixels] = find(impacted_pixels);

if isempty(impacted_pixels)
    warning('No new unknown pixel inpainted.');
    impacted_patches = [];
else
    impacted_patches = repmat(patch_pixels, length(impacted_pixels), 1);
    impacted_patches = impacted_patches + repmat(impacted_pixels.', 1, length(impacted_patches));
    impacted_patches = reshape(impacted_patches, 1, []);
    impacted_patches = unique(impacted_patches);
    impacted_patches = impacted_patches + vertex;
end

%% Inpaint the patch

% Retrieve pixel values.
switch(gparam.inpainting.retrieve)
    % Weighted average of connected patches.
    case 'average'
        new = G.W(vertex,:) * patches / norm(G.W(vertex,:),1);
    % Copy of the strongest connected patch.
    case 'copy'
        [~,strongest] = max(G.W(vertex,:));
        new = patches(strongest,:);
end

% Compose the new patch.
switch(gparam.inpainting.compose)
    % Keep valid pixels, only inpaint unknown ones.
    case 'mixed'
        M = patches(vertex,:)>=0;
    % Inpaint the entire patch, i.e. replace known pixels.
    case 'overwrite'
        M = [zeros(1,gparam.graph.psize^2), ones(1,2)];
end

% Inpaint the patch.
old = patches(vertex,:);
patches(vertex,:) = old .* M + new .* ~M;

%% Update the signals.

% Update pixel values.
pixels(vertex+patch_pixels) = patches(vertex,1:end-2).';

% Update pixels information priorities.
new = Pinformation(vertex,2);
new = repmat(new, length(patch_pixels), 1);
old = Pinformation(vertex+patch_pixels,1);
Pinformation(vertex+patch_pixels,1) = old .* M(1:end-2).' + new .* ~M(1:end-2).';

% Patch vertices affected by a patch.
patch_patches = giin_patch_vertices('patches', gparam.graph.psize, height);

% Update neighboring (not necessarily impacted) patches.
% dim = gparam.psize^2;
% margin = floor(gparam.psize / 2);
% inpainted = reshape(pixels,height,width);
for patch = vertex+patch_patches
%     h = mod(patch-1, height) + 1;
%     w = floor((patch-1) / height) + 1;
%     patches(patch, 1:dim) = reshape(inpainted(h-margin:h+margin, w-margin:w+margin), 1, dim);
    patches(patch,1:end-2) = pixels(patch+patch_pixels);
    % Update patches information priorities.
    Pinformation(patch,2) = mean(Pinformation(patch+patch_pixels,1));
end

% Execution time.
% fprintf('giin_inpaint : %f seconds\n', toc);

end