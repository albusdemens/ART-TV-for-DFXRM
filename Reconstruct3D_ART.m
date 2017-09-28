close all; clear;

% Load as input output files from reconstr.m
F1 = load('rec_ART.mat');
R = F1.rec;

R_bin = zeros(size(R));
for i = 1:size(R,3)
    disp(i),
    layer = squeeze(R(:,:,i));
    layer1 = layer;
    % Threshold and binarize
    layer1(layer1 < 45) = 0;
    layer2 = zeros(size(layer1));
    layer2(layer1 > 0) = 1;
    % Fill holes
    layer2 = imfill(layer2);
    % To get rid of artefacts, erode and dilate
    se = strel('disk',10);
    l_ed = imdilate(imerode(layer2,se),se);

    R_bin(:,:,i) = l_ed(:,:);
    % Perimeter overlay
    % P = bwperim(l_ed);
    % figure; imagesc(imoverlay(layer, P, [.3 1 .3]));
end

X = zeros(nnz(R_bin), 3);
num_p = 0;
for aa = 1:size(R,1)
    for bb = 1:size(R,2)
        for cc = 1:size(R,3)
            if R_bin(aa,bb,cc) > 0
                num_p = num_p + 1;
                X(num_p, 1) = aa;
                X(num_p, 2) = bb;
                X(num_p, 3) = cc;
            end
        end
    end
end

%figure; scatter3(X(:,1), X(:,2), X(:,3));
%xlabel('X'); ylabel('Y'); zlabel('Z');

S = zeros(size(R_bin,1), size(R_bin,2));
for i = 1:size(R_bin,2)
    for j = 1:size(R_bin,1)
        for k = 1:size(R_bin,3)
            S(j,k) = S(j,k) + R_bin(j,i,k);
        end
    end
end

figure; h = pcolor(S); shading flat;

save('Binary_vol_ART.mat', 'R_bin');

% Save for Paraview
R_bin_vtk = zeros(size(R_bin));
% Resize volume from 300x300x300 to 100x100x100
R_bin_1 = imresize(R_bin, 1/3);
R_bin_2 = zeros(100,100,100);
for i = 1:size(R_bin_1,1)
    Lay = squeeze(R_bin_1(i,:,:));
    Lay_1 = imresize(Lay, [100 100]);
    R_bin_2(:,:,i) = Lay_1(:,:);
end
R_bin_vtk = R_bin_2;
R_bin_vtk(R_bin_2 == 0) = NaN;
savevtk(R_bin_vtk, 'Vol_Art.vtk');