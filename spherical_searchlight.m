function search_indexes = spherical_searchlight(xyzcoords,radius)
% search_indexes = volume_searchlight(xyzcoords,radius)
% 
%   returns spherical searchlights indexs within a brain volume
%
%   Input
%       xyzcoords >>> matrix of xyz coordinates in the form [n 3]: output
%                     of mat2coords
%       radius >>> the desired radius 
%
%   Output:
%       search_indexes >>> a cell indicating which voxels contributes to
%                          each of the {n} searchlights
%
%   2018 - Paolo Papale fecit

sphere_diam = radius*2+1;
mask_center = ceil(sphere_diam/2);
[X Y Z] = meshgrid(1:sphere_diam, 1:sphere_diam, 1:sphere_diam);
sphere_mask = sqrt((X-mask_center).^2+(Y-mask_center).^2+(Z-mask_center).^2)<=radius;
disp(sprintf('Computing searchlights...'))
parfor i = 1:size(xyzcoords,1)
    center = xyzcoords(i,:);
    X_temp = X+(center(1)-radius)-1;
    Y_temp = Y+(center(2)-radius)-1;
    Z_temp = Z+(center(3)-radius)-1;
    cube_voxels = [X_temp(:) Y_temp(:) Z_temp(:)];
    sphere_voxels = cube_voxels.*sphere_mask(:);
    [~,goods_temp] = intersect(xyzcoords,sphere_voxels,'rows');
    search_indexes{i} = goods_temp;
end
