% This script combines the reconstruction from the first 78 deg and from
% the last 60 deg

% Alberto Cereser, DTU Fysik
% September 2017
% alcer@fysik.dtu.dk

close all; clear;

% Load as input output files from reconstr.m
F1 = load('rec.mat');
F2 = load('rec_2.mat');

R_1 = F1.rec;
R_2 = F2.rec_2;

R_comb = zeros(size(R_1));
for i = 1:size(R_1,3)
    for j = 1:size(R_1,1)
        for k = 1:size(R_1,2)
            R_comb(j,k,i) = R_1(j,k,i) + R_2(j,k,i);
        end
    end
end

figure;
subplot(1,3,1); h = pcolor(squeeze(R_1(:,:,200))); shading flat;
title('Slice from first reconstruction', 'FontSize', 20);
hold on;
subplot(1,3,2); h = pcolor(squeeze(R_2(:,:,200))); shading flat;
title('Slice from second reconstruction', 'FontSize', 20);
hold on;
subplot(1,3,3); h = pcolor(squeeze(R_comb(:,:,200))); shading flat;
title('Combined slice', 'FontSize', 20);
hold on;
