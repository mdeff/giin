function [ pixels, patches, impacted_patches ] = giin_inpaint( vertex, G, pixels, patches, gparam )
%GIIN_INPAINT Inpaint the specified patch and returns the impacted vertices.
%   Detailed explanation goes here

tic;

%% Impacted patches, i.e. patches that should be investigated next
% List of impacted patches, i.e. patches where an unknown pixel is now
% known. If the value is solely modified, it's not impacted (otherwise we
% would never be able to stop).

width = max(G.coords(:,1));
height = max(G.coords(:,2));

% Pixel vertices contained in a patch.
patch_pixels = patch_vertices('pixels', gparam.psize, height);

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
switch(gparam.inpainting_retrieve)
    % Weighted average of connected patches.
    case 'average'
        new = G.W(vertex,:) * patches / norm(G.W(vertex,:),1);
    % Copy of the strongest connected patch.
    case 'copy'
        [~,strongest] = max(G.W(vertex,:));
        new = patches(strongest,:);
end

% Compose the new patch.
switch(gparam.inpainting_compose)
    % Keep valid pixels, only inpaint unknown ones.
    case 'mixed'
        M = patches(vertex,:)>=0;
    % Inpaint the entire patch, i.e. replace known pixels.
    case 'overwrite'
        M = [zeros(1,gparam.psize^2), ones(1,2)];
end

% Inpaint the patch.
old = patches(vertex,:);
patches(vertex,:) = old .* M + new .* ~M;

%% Update the signals.

% Update pixel values.
pixels(vertex+patch_pixels) = patches(vertex,1:end-2).';

% Patch vertices affected by a patch.
patch_patches = patch_vertices('patches', gparam.psize, height);

% Update neighboring (not necessarily impacted) patches.
dim = gparam.psize^2;
margin = floor(gparam.psize / 2);
inpainted = reshape(pixels,height,width);
for patch = vertex+patch_patches
    h = mod(patch-1, height) + 1;
    w = floor((patch-1) / height) + 1;
    patches(patch, 1:dim) = reshape(inpainted(h-margin:h+margin, w-margin:w+margin), 1, dim);
end

% Execution time.
% fprintf('giin_inpaint : %f seconds\n', toc);

end

function [ vertices ] = patch_vertices( type, psize, height )
%PATCH_VERTICES Return a list of vertices impacted by a patch. Either the
%impacted pixels or the impacted patches.

switch(type)
    case 'pixels'
        size = psize;
    case 'patches'
        size = 2 * psize - 1;
end

% Pixel vertices contained in a patch.
vert = -(size-1)/2 : (size-1)/2;
horiz = vert * height;
vert = repmat(vert.', 1, size);
horiz = repmat(horiz, size, 1);
vertices = horiz + vert;
vertices = reshape(vertices, 1, size*size);

end