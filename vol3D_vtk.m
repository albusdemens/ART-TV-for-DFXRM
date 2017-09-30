% This script compares the reconstructed 3D data with the sum of the
% intensities recorded at different projections

close all; clear;

cutout_completeness = 0;

addpath('/npy_matlab_master/');
% Read reconstructed volume. Format: X, Y, Z, param. Parameters: gamma, mu,
% completeness
V = readNPY('/u/data/alcer/DFXRM_rec/Rec_test_2/grain_ang.npy');

% Volume selected using completeness
V_th = zeros(size(V,1), size(V,2), size(V,3));
% Volume with angular values
V_th_mos = zeros(size(V,1), size(V,2), size(V,3), 3);

for ii =1:size(V,1)
    for jj = 1:size(V,2)
        for kk = 1:size(V,3)
            % The minimum completeness value for a voxel
            % to be part of the volume is 0.5
            if V(ii,jj,kk,3) > cutout_completeness
                V_th(ii,jj,kk) = V(ii,jj,kk,3);
                V_th_mos(ii,jj,kk,:) = V(ii,jj,kk,:);
            end
        end
    end
end

% Save the selected region
save('V_mos_recon3d.mat', 'V_th_mos');

% Rescale, se we can compare with the reocnstruction from ART+TV
V_resc = zeros((size(V,1) * 3) -3, (size(V,2) * 3) -3, (size(V,3) * 3) - 3);
for jj = 1:(size(V,3) - 1)
    Layer = squeeze(V_th(1:100,1:100,jj));
    Layer_resc = imresize(Layer, 3);
    V_resc(:,:,jj) = Layer_resc(:,:);
end

% Save calculated volumes
savevtk(V_th, '/u/data/alcer/DFXRM_rec/Rec_test_2/V_th.vtk');
