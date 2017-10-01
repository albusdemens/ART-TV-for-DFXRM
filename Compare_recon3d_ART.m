% Script to compare the volumes reconstructed using recon3d and ART+TV. It
% also returns the orientation distribution for the reconstrcuted volume

close all; clear;

addpath('/npy_matlab_master/'); % Command required to lead npy files. If 
% this doesn't work, add the npy folder using "Add to pach --> selected 
% folder and subfolders"

% Read reconstructed volume. Format: X, Y, Z, param. Parameters: gamma, mu,
% completeness
Summed_img = readNPY('summed_data_astra.npy');
% Load as input volume reconstructed using ART-TV in Reconstruct3D_ART.m
Vol = load('Binary_vol_ART_all_p.mat');
V = Vol.R_bin;
% Load volume reconstructed using recon3d
Vol_recon = load('V_mos_recon3d.mat');
V_recon = Vol_recon.V_th_mos;

% Resize V_recon to have the same dimensions of V
V_r3d_1 = zeros(300,300,100);
V_r3d_2 = zeros(300,300,300);
for i = 1:100
    V_r3d_1(:,:,i) = imresize(squeeze(V_recon(1:100,1:100,i,3)), [300 300]);
end
for j = 1:300
    V_r3d_2(j,:,:) = imresize(squeeze(V_r3d_1(j,:,:)), [300 300]);
end

% Rotate the ART reconstructed volume by 90 deg around X (ParaView
% geometry) to be consistent with the summed diffraction orientation
V_rot = zeros(size(V));
V_rot_1 = zeros(size(V));
for ii = 1:size(V,1)
    Slice = zeros(size(V,1), size(V,2));
    Slice = imrotate(squeeze(V(:,:,ii)), 90);
    V_rot(:,:,ii) = Slice(:,:);
end
for ii = 1:size(V,2)
    Slice = zeros(size(V,1), size(V,2));
    Slice = flipud(imrotate(squeeze(V_rot(ii,:,:)), 90));
    V_rot_1(ii,:,:) = Slice(:,:);
end

% Rotate the recon3d geometry to make it consistent with the other ones
V_recon_2 = zeros(size(V_recon));
for ii = 1:size(V_recon,1)
    Slice_3_i = zeros(size(V_recon,1), size(V_recon,2));
    Slice_3_j = zeros(size(V_recon,1), size(V_recon,2));
    Slice_3_k = zeros(size(V_recon,1), size(V_recon,2));
    Slice_3_i = imrotate(squeeze(V_recon(ii,:,:,1)), 270);
    Slice_3_j = imrotate(squeeze(V_recon(ii,:,:,2)), 270);
    Slice_3_k = imrotate(squeeze(V_recon(ii,:,:,3)), 270);
    V_recon_2(ii,:,:,1) = Slice_3_i(:,:);
    V_recon_2(ii,:,:,2) = Slice_3_j(:,:);
    V_recon_2(ii,:,:,3) = Slice_3_k(:,:);
end

V_r3d_4 = zeros(size(V));
for ii = 1:size(V,2)
    Slice = zeros(size(V,1), size(V,2));
    Slice = imrotate(squeeze(V_r3d_2(ii,:,:)), 270);
    V_r3d_4(ii,:,:) = Slice(:,:);
end

% Shift down the recon3d volume
V_r3d_3 = zeros(size(V));
for jj = 1:(size(V,2)-20)
    V_r3d_3(:,jj,:) = V_r3d_4(:,jj+20,:);
end

V_recon_3 = zeros(100,100,100,3);
for jj = 1:100-7
    V_recon_3(:,jj,:,:) = V_recon_2(1:100,jj+7,1:100,:);
end

% Save the two rotated volumes
%save('V_rot_ART.mat', 'V_rot_1');
%save('V_rot_recon3d.mat', 'V_r3d_3');

% Find the intersection between the ART reconstruction and the recon3d
% reconstruction
V_intersect = zeros(size(V));
for ii = 1:size(V,1)
    for jj = 1:size(V,2)
        for kk = 1:size(V,3)
            if V_rot_1(ii,jj,kk) > 0
                V_intersect(ii,jj,kk) = V_r3d_3(ii,jj,kk);
            end
        end
    end    
end

% Rescale V_intersect fom (300,300,300) to (100,100,100), so that it can be
% visualized in paraview
V_intersect_1 = zeros(100,100,300);
V_intersect_2 = zeros(100,100,100);
for i = 1:300
    V_intersect_1(:,:,i) = imresize(squeeze(V_intersect(1:300,1:300,i)), [100 100]);
end
for j = 1:100
    V_intersect_2(j,:,:) = imresize(squeeze(V_intersect_1(j,:,:)), [100 100]);
end

%%% First reconstrcution modality: set completeness thershold, find
%%% intersection between reconstructed volumes

% % We still require completeness to be greater than 0.5
% V_intersect_2(V_intersect_2 < 0.5) = 0;
% Using the intersetction of the volume, we find the corresponding 4-dim
% volume, with information on mosaicity too
for aa = 1:size(V_recon_3,1)
    for bb = 1:size(V_recon_3,2)
        for cc = 1:size(V_recon_3,3)
            if V_intersect_2(aa,bb,cc) == 0
                V_recon_3(aa,bb,cc,:) = 0;
            end
        end
    end
end

% Rescale the mu and gamma values to [0,1]
% Maybe also include opticolors
V_rec_3_resc = zeros(size(V_recon_3));
range_mu = max(max(max(V_recon_3(:,:,:,1)))) - min(min(min(V_recon_3(:,:,:,1))));
range_gamma = max(max(max(V_recon_3(:,:,:,2)))) - min(min(min(V_recon_3(:,:,:,2))));

list_min_mu = zeros(nnz(V_recon_3),1);
mu_el_n = 0;
for i = 1:size(V_recon_3,1)
    for j = 1:size(V_recon_3,2)
        for k = 1:size(V_recon_3,3)
            if V_recon_3(i,j,k) > 0
                mu_el_n = mu_el_n + 1;
                list_min_mu(mu_el_n) = V_recon_3(i,j,k);
            end
        end
    end
end

list_min_mu = list_min_mu(list_min_mu > 0);
min_mu = min(list_min_mu);
max_mu = max(max(max(V_recon_3(:,:,:,1))));
min_gamma = min(min(min(V_recon_3(:,:,:,2))));

%M = unique(V_recon_3(:,:,:,1));
%G = unique(V_recon_3(:,:,:,2));

%CM = opticolor(M,min_mu,max_mu);

for aa = 1:size(V_recon_3,1)
    for bb = 1:size(V_recon_3,2)
        for cc = 1:size(V_recon_3,3)
            if (V_recon_3(aa,bb,cc,3)) > 0
                V_rec_3_resc(aa,bb,cc,1) = (V_recon_3(aa,bb,cc,1) - min_mu)/range_mu;
                V_rec_3_resc(aa,bb,cc,2) = (V_recon_3(aa,bb,cc,2) - min_gamma)/range_gamma;
            end
        end
    end
end

% For each voxel, calculate the correspnding HVS color. 
V_rec_3_HVS = zeros(size(V_recon_3));
for aa = 1:size(V_recon_3,1)
    for bb = 1:size(V_recon_3,2)
        for cc = 1:size(V_recon_3,3)
            if (V_recon_3(aa,bb,cc,3)) > 0
                vr = V_rec_3_resc(aa,bb,cc,1);
                ur = V_rec_3_resc(aa,bb,cc,2);
                V_rec_3_HVS(aa,bb,cc,1) = wrapTo2Pi(atan2(vr,ur))/pi/2;
                V_rec_3_HVS(aa,bb,cc,2) = sqrt(ur^2+vr^2)/sqrt(2);
                V_rec_3_HVS(aa,bb,cc,3) = 1;
            end
        end
    end
end

for i = 1:size(V_recon_3,1)
    gg = figure;
    image(hsv2rgb(imrotate(squeeze(V_rec_3_HVS(:,i,:,:)),90)));
    saveas(gg,sprintf('Mosaicity_plot/mosaicity_layer_%03i.png', i), 'png');
    close;
end

% Flip before saving for ParaView
V_intersect_3 = zeros(size(V_intersect_2));
for i = 1:size(V_intersect_2,1)
    V_intersect_3(i,:,:) = fliplr(squeeze(V_intersect_2(i,:,:)));
end

savevtk(V_intersect_3, 'V_intersect_2.vtk');

% Test rotation of ART volume
S = zeros(size(V_rot,1), size(V_rot,2));
S_r3d = zeros(size(V_rot,1), size(V_rot,2));
S_intersect = zeros(size(V_rot,1), size(V_rot,2));
for i = 1:size(V_rot,1)
    for j = 1:size(V_rot,2)
        for k = 1:size(V_rot,3)
            S(j,k) = S(j,k) + V_rot_1(i,j,k);
            S_r3d(j,k) = S_r3d(j,k) + V_r3d_3(i,j,k);
            S_intersect(j,k) = S_intersect(j,k) + V_intersect(i,j,k);
        end
    end
end

figure;
subplot(2,2,1); h = pcolor(S); shading flat;
xlabel('Y'); ylabel('Z'); title('Projection of reconstructed volume');
hold on;
subplot(2,2,2); H = pcolor(S_r3d); shading flat;
hold on;
subplot(2,2,3); h = pcolor(squeeze(S_intersect)); shading flat;
hold on;
subplot(2,2,4); h = pcolor(squeeze(Summed_img(:,1,:))); shading flat;
xlabel('Y_d'); ylabel('Z_d'); title('Sum of diffracted signal');

% Every ten degrees, compare the projection of the reconstructed volume
% with the sum of the frames collected at the corresponding angle. When
% doing so, remember that the frames were collected with two angular steps
% (first 0.8 and then 1.6 deg), and that there is a gap in the acquisition

% Proj: list projection number and corresponding degree
Proj = zeros(21,2);
% Projection number
Proj(:,1) = [1, 6, 13, 19, 25, 31, 37, 44, 50, 56, 62, 69, 75, 124, ...
    133, 137, 142, 147, 151, 155, 160];
% Degrees
Proj(:,2) = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 140, ...
    150, 155, 160, 165, 170, 175, 180];

for i = 1:size(Proj, 1)
    Compare_shapes(-Proj(i,2), Proj(i,1), V, V_rot_1, V_r3d_3, Summed_img);
end

%%% This function compares, for a selected angle, the projection of the
% reconstructed volume with the sum of the measured diffarction signal
function Compare_shapes(angle_num, p_num, V, V_rot_1, V_r3d_3, Summed_img)

% Rotate the reconstructed volume by the selected angle
angle = angle_num;
A_1 = imrotate(squeeze(V(:,1,:)), angle);
Rot_vol = zeros(size(A_1,1), size(V_rot_1,2), size(A_1,2));
Rot_V_r3d = zeros(size(A_1,1), size(V_rot_1,2), size(A_1,2));
for ii = 1:size(V_rot_1,2)
    Rot_vol(:,ii,:) = imrotate(squeeze(V_rot_1(:,ii,:)), angle);
    Rot_V_r3d(:,ii,:) = imrotate(squeeze(V_r3d_3(:,ii,:)), angle);
end

%figure; h= pcolor(squeeze(Rot_vol(:,:,150)));

Proj = zeros(size(Rot_vol,2), size(Rot_vol,3));
Proj_r3d = zeros(size(Rot_vol,2), size(Rot_vol,3));
for aa = 1:size(Rot_vol,1)
    for bb = 1:size(Rot_vol, 2)
        for cc = 1:size(Rot_vol, 3)
            Proj(bb,cc) = Proj(bb,cc) + Rot_vol(aa,bb,cc);
            Proj_r3d(bb,cc) = Proj_r3d(bb,cc) + Rot_V_r3d(aa,bb,cc);
        end
    end
end

% Resize summed image (cut out border introduced by imrotate)
Proj_clean = zeros(size(V,1), size(V,2));
Proj_r3d_clean = zeros(size(V,1), size(V,2));
frame_sz = round((size(A_1,1) - size(V,1))/2);
for i = frame_sz + 1:size(Proj,2) - frame_sz
        Proj_clean(:,i - frame_sz) = Proj(:,i);
        Proj_r3d_clean(:,i - frame_sz) = Proj_r3d(:,i);
end

%%% Binarize volume from ART+TV and find perimeter
% Binarize projected image and find perimeter
Proj_clean_bin = zeros(size(Proj_clean));
Proj_bin_e = zeros(size(Proj_clean));
Proj_sum_p = zeros(size(Proj_clean));
Proj_clean_bin(Proj_clean > 0) = 1;
% Erode the projection, to take into account of possible XY dilations
se = strel('disk',10);
Proj_bin_e = imerode(Proj_clean_bin, se);
% Find perimeter
P_e = bwperim(Proj_bin_e);

%%% Binarize volume from recon3d and find perimeter
% Binarize projected image and find perimeter
Proj_clean_bin_r3d = zeros(size(Proj_r3d_clean));
Proj_bin_e_r3d = zeros(size(Proj_r3d_clean));
Proj_sum_p_r3d = zeros(size(Proj_r3d_clean));
Proj_clean_bin_r3d(Proj_r3d_clean > 0) = 1;
% Erode the projection, to take into account of possible XY dilations
se = strel('disk',10);
Proj_bin_e_r3d = imerode(Proj_clean_bin_r3d, se);
% Find perimeter
P_e_r3d = bwperim(Proj_bin_e_r3d);

% Draw perimeter of reconstructed projection over summed diffraction
Proj_sum = squeeze(Summed_img(:,p_num,:));
Proj_sum_p = Proj_sum;
Proj_sum_p(P_e > 0) = max(max(Proj_sum));
Proj_sum_p(P_e_r3d > 0) = max(max(Proj_sum));

% Compare the projection of the reconstructed volume with the sum of the
% diffraction signal
a = figure;
%subplot(1,,1);
%h = pcolor(Proj_sum); shading flat;
%title('Sum of diffracted signal');
%hold on;
subplot(1,3,1);
h = pcolor(Proj_clean); shading flat;
title('ART+TV projection');
hold on;
subplot(1,3,2);
h = pcolor(Proj_r3d_clean); shading flat;
title('recon3d projection');
hold on;
subplot(1,3,3);
h = pcolor(Proj_sum_p); shading flat;
title('Sum diff signal');

saveas(a, sprintf('Shape_comparison/Proj_%03i.png', angle_num),'png');
close;
end
