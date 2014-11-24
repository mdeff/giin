function [ vertices ] = giin_patch_vertices( type, psize, height )
%GIIN_PATCH_VERTICES Return a list of vertices impacted by a patch. Either the
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