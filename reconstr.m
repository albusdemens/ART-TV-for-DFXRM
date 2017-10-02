% Script to reconstruct a tomography dataset using ART+TV

% Matteo Busi, DTU Fysik
% September 2017
% mbusi@fysik.dtu.dk

clear; close all;

%% Data processing
slice = 150;
load data.mat;
figure, imagesc(squeeze(Fab_resc(:,slice,:)))
title('Starting sinogram')
theta = 0:0.8:180;

deadvalue = -2000;
new_sino(1:96,:,:) = Fab_resc(1:96,:,:); % first data interval
new_sino(97:147,:,:) = deadvalue ; %empty data
for i = 1:35 %second data interval
    new_sino(2*(i-1)+148,:,:) = Fab_resc(122+i,:,:);
    new_sino(2*(i-1)+148+1,:,:) = deadvalue ;
end
new_sino(217:220,:,:) = deadvalue; %last missing projection
for i = 1:4
    new_sino(2*(i-1)+221,:,:) = Fab_resc(157+i,:,:);
    new_sino(2*(i-1)+221+1,:,:) = deadvalue ;
end

%% Reconstruction
ss = 30;                                              % Target area (cm)
geostruct.range_angle = 180;                          % Angular span
geostruct.nproj = size(new_sino,1);                         % Number of simulated projectons
geostruct.ndet = size(new_sino,2);                          % Number of detector pixels
geostruct.nElem = 1;                                  % Number of detectors
geostruct.Sep = 0;                                    % Pixels' gap
geostruct.det_space = 0.1*geostruct.ndet;           % Size of detector in cm (pixels*pixel size)
geostruct.model = 'par';                              % Beam geometry (par-fan-lshape-lshape1)
geostruct.imagesize = [300 300];                      % Size of output image grid
geostruct.delta = ss/geostruct.imagesize(1);          % Reconstruction resolution
geostruct.SDD = 600;                               % Source-Detector distance
geostruct.sourceCentShift = 0;                        % Vertical source shift from perfect placement
geostruct.detectCentShift = 0;                        % Vertical detector shift from
geostruct.SAD = 300;                                % Source-AxisOfRotation distance

%algorithm parameters
param.lambda = 0.8;                                    % ART parameter
param.alpha = 0.2;                                     % TV parameter
param.niter = [250 20];                               % [ART TV] # of iterations
param.stochast = false;                               % Pick random iteration order if true
deadThresh = -2000;                                    % value to set det pixels
thresh = -3e-2;                                         % Any values lower than or equal to thresh are ignored during reconstruction
backVal = 0;                                          % Background mask

for i=1:300
    proj = squeeze(new_sino(:,i,:))';
    rec(:,:,i) = ART_TV_reconstruct_2d_new(proj,geostruct,param,deadThresh,thresh,backVal);
    %title('Reconstruction');
    %pause(0.01)
    progress =  100*i/300;
    disp(progress);
end
save rec.mat rec
