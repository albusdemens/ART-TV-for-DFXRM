% This scipt compares the reconstructed 3D data with the sum of the
% intensities recorded at different projections

close all; clear;

addpath('/npy_matlab_master/');
% Read reconstructed volume. Format: X, Y, Z, param. Parameters: gamma, mu,
% completeness
Summed_img = readNPY('summed_data_astra.npy');
% Load as input volume reconstructed using ART-TV in Reconstruct3D_ART.m
Vol = load('Binary_vol_ART.mat');
V = Vol.R_bin;

% Rotate the reconstructed volume by 90 deg around X (ParaView geometry) to
% be consistent with the summed diffraction orientation

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

% Test rotation

S = zeros(size(V_rot,1), size(V_rot,2));
for i = 1:size(V_rot,1)
    for j = 1:size(V_rot,2)
        for k = 1:size(V_rot,3)
            S(j,k) = S(j,k) + V_rot_1(i,j,k);
        end
    end
end

%figure;
%subplot(1,2,1);h = pcolor(S); shading flat;
%xlabel('Y'); ylabel('Z'); title('Projection of reconstructed volume');
%hold on;
%subplot(1,2,2); h = pcolor(squeeze(Summed_img(:,1,:))); shading flat;
%xlabel('Y_d'); ylabel('Z_d'); title('Sum of diffarcted signal');

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
    Compare_shapes(-Proj(i,2), Proj(i,1), V, V_rot_1, Summed_img);
end

%%% This function compares, for a selected angle, the projection of the
% reconstructed volume with the sum of the measured diffarction signal
function Compare_shapes(angle_num, p_num, V, V_rot_1, Summed_img)

% Rotate the reconstructed volume by the selected angle
angle = angle_num;
A_1 = imrotate(squeeze(V(:,1,:)), angle);
Rot_vol = zeros(size(A_1,1), size(V_rot_1,2), size(A_1,2));
for ii = 1:size(V_rot_1,2)
    Rot_vol(:,ii,:) = imrotate(squeeze(V_rot_1(:,ii,:)), angle);
end

%figure; h= pcolor(squeeze(Rot_vol(:,:,150)));

Proj = zeros(size(Rot_vol,2), size(Rot_vol,3));
for aa = 1:size(Rot_vol,1)
    for bb = 1:size(Rot_vol, 2)
        for cc = 1:size(Rot_vol, 3)
            Proj(bb,cc) = Proj(bb,cc) + Rot_vol(aa,bb,cc);
        end
    end
end

% Resize summed image (cut out border introduced by imrotate)
Proj_clean = zeros(size(V,1), size(V,2));
frame_sz = round((size(A_1,1) - size(V,1))/2);
for i = frame_sz + 1:size(Proj,2) - frame_sz
    Proj_clean(:,i - frame_sz) = Proj(:,i);
end

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

Proj_sum = squeeze(Summed_img(:,p_num,:));
Proj_sum_p = Proj_sum;
Proj_sum_p(P_e > 0) = max(max(Proj_sum));

% Compare the projection of the reconstructed volume with the sum of the
% diffraction signal
a = figure;
%subplot(1,,1);
%h = pcolor(Proj_sum); shading flat;
%title('Sum of diffracted signal');
%hold on;
subplot(1,2,1);
h = pcolor(Proj_clean); shading flat;
title('Projection of reconstructed volume');
hold on;
subplot(1,2,2);
h = pcolor(Proj_sum_p); shading flat;
title('Sum diff signal with projection perimeter');

saveas(a, sprintf('Shape_comparison/Proj_%03i.png', angle_num),'png');
close;
end
