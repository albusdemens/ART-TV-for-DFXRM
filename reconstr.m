clear; close all;

load data.mat

ss = 30;                                              % Target area (cm)
geostruct.range_angle = 96*180/225;                          % Angular span
geostruct.nproj = 96;                         % Number of simulated projections
geostruct.ndet = size(Fab_resc,2);                          % Number of detector pixels
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
    proj = squeeze(Fab_resc(:,i,:))';
    %figure;
    %subplot(1,2,1);
    %imagesc(proj), colorbar
    proj(:,97:161)=[];
    rec(:,:,i) = ART_TV_reconstruct_2d_new(proj,geostruct,param,deadThresh,thresh,backVal);
    %subplot(1,2,2);
    %imagesc(rec(:,:,i)), colorbar, caxis([0 500])
    %pause(0.01)
    progress =  100*i/300;
    disp(progress);
end
save rec.mat rec
