function [f_out,itOut] = ART_TV_reconstruct_2d_new(proj_data,geoStruct, param,deadThresh,thresh,backVal,imStart)
%
%Input:
% proj_data: -ln(I/I0) data of size nPixels x nAngles
% geoStruct: Struct containing geometric X-ray information
% param: Struct containing optimization related parameters
% thresh: Any values lower than or equal to thresh are ignored during reconstruction
%
%Output:
%f_out: The reconstructed image slice
%
%Follow this paper:
% Accurate image reconstruction from few-views and limited-angle data in divergent-beam CT.
% J. X-ray Sci. Tech. 2006a;14:119–139.
%
if (~strcmp(geoStruct.model,'lshape') && ~strcmp(geoStruct.model,'lshape1'))
    
    range_angle=geoStruct.range_angle;
    num_ang=geoStruct.nproj;        % Number of projections
    num_det=geoStruct.ndet;         % Number of detector pixels
    det_space=geoStruct.det_space;  % The size of the detector in cm
    det_sample=det_space/num_det;   % Sampling of the detector line
    delta=geoStruct.delta;          % The pixel size of the input attenuation volume to reconsturcut
    
    imageSize=geoStruct.imagesize;%Size that is returned
    
    detect=zeros(num_det,2); %(x,y) coordinates of
    %if the detector is made of several 1d detectors with spacings between
    if geoStruct.nElem>1
        %take care of between detector plate spacing
        det_sample=(det_space-(geoStruct.nElem-1)*geoStruct.Sep)/num_det;   % Sampling of the detector line
        plateSize=num_det/geoStruct.nElem;
        jump=plateSize:plateSize:num_det;
        shiftVal=0;
        
        for k=1:num_det
            detect(k,1)=-det_space/2-det_sample/2 + k*det_sample+shiftVal*geoStruct.Sep;
            detect(k,2)=geoStruct.SDD-geoStruct.SAD; % just suffiently far from the rotation axis origin to +infty
            
            if(k==jump(1))
                jump(1)=[];%remove it
                shiftVal=shiftVal+1;
            end
        end
        
    else%or if we have one long detector of pixel elements
        for k=1:num_det
            detect(k,1)=-det_space/2-det_sample/2 + k*det_sample;
            detect(k,2)=geoStruct.SDD-geoStruct.SAD; % just suffiently far from the rotation axis origin to +infty
            
        end
        
    end
    
    %if the detector origin is shifted add the shift
    detect(:,1)=detect(:,1)+geoStruct.detectCentShift;
    
    if (strcmp(geoStruct.model,'par'))
        % Sources positions
        source=zeros(num_det,2); %(x,y) coordinates of
        for k=1:num_det
            source(k,1)=-det_space/2-det_sample/2 + k*det_sample;
            source(k,2)=-geoStruct.SAD; % just suffiently far from the origin to -infty
        end
        %shift the source up/down
        source(:,1)=source(:,1)+geoStruct.sourceCentShift;
        % Estimate paralel geometry for one column of the detector, all other
        % columns are translates of this column
        [pointLOR, num_data] = par_beam_geometry(range_angle, num_ang, num_det, detect, source);
    end
    
    if (strcmp(geoStruct.model,'fan'))
        
        
        % Sources positions
        source=zeros(num_det,2); %(x,y) coordinates of
        for k=1:num_det
            source(k,1)=0;
            source(k,2)=-geoStruct.SAD; % just suffiently far from the origin to -infty
        end
        %shift the source up/down
        source(:,1)=source(:,1)+geoStruct.sourceCentShift;
        % Estimate paralel geometry for one column of the detector, all other
        % columns are translates of this column
        [pointLOR, num_data] = par_beam_geometry(range_angle, num_ang, num_det, detect, source);
    end
    
else%else if we have special L-shape
     if(strcmp(geoStruct.model,'lshape'))
    num_ang=geoStruct.nproj;
    delta=geoStruct.delta;          % The pxel size of the input attenuation volume
    imageSize=geoStruct.imagesize;% the reconstruction size
    [pointLOR,num_det] = findPointLORfromgeometry(num_ang,geoStruct.SADmove,[]);
    
%         %mergefactor every nth pixel
%     nthPix=1;%
%     start=1;
%     p2lor=[];
%     for jj=1:numel(num_det)
%        adj=num_det(jj);
%        temp2=pointLOR(start:adj,:);
%        while (mod(size(temp2,1),nthPix)~=0)
%        adj=adj-1;
%     
%        temp2=pointLOR(start:adj,:);
%        end
%        
%        start=num_det(jj)+1;      
%        p2lor=[p2lor squeeze(mean(reshape(temp2,nthPix,4,[]),1))];
%     end
%     
%     pointLOR=p2lor';
    
    num_data=size(pointLOR,1);
     end
     
      if(strcmp(geoStruct.model,'lshape1'))
              num_ang=geoStruct.nproj;
    delta=geoStruct.delta;          % The pxel size of the input attenuation volume
    imageSize=geoStruct.imagesize;% the reconstruction size
    [pointLOR,num_det] = findPointLORfromgeometry2(num_ang,geoStruct.SADmove,[]);
    num_data=size(pointLOR,1);
      end
end


%The iteration index is either angular sequantial/random or fixed angular mixed
itIDX=1:num_data;
if param.stochast==1%should we go trogh geometry i random order
    itIDX=itIDX(randperm(numel(itIDX)));
end

if param.stochast==2%should we go through geometry i random order
    temp=itIDX(1:round(num_det/2):end);
    itIDX=[];
    for i=0:(round(num_det/2)-1)
        itIDX=rem([itIDX [i+temp]],num_data);
    end
    itIDX(itIDX==0)=num_data;
end
proj_data=proj_data(:);

%% Time to do ART+TV
tol=0.00001;
alpha=param.alpha;      %TV penalty smoother
nIter=param.niter(1);   %Number of ART
nTV=param.niter(2);     %Number of TV
epsTV=0.0000001;        %Constant value set

if nargin==7
    f_out = imStart;    %User supplied start
    f_inp=f_out;        %ART initial guess
else
    f_out = zeros(imageSize(1),imageSize(2));%Allows reconstruction with differing (x,y) dims
    f_inp=f_out; %ART initial guess
end
try
    v=param.moment;%if this does not exist-> no momentum used
    useMoment=1;
catch
    param.moment=0;%if not set by user we set it to 0 to emulate no moment usage
    useMoment=0;
end
lambda=param.lambda;    %Learning rate
%vMoment={};%zeros(imageSize(1),imageSize(2));
vMomentF=zeros(imageSize(1),imageSize(2));
vCount=vMomentF;
f_old=f_inp;
binMask=zeros(imageSize(1),imageSize(2));%mask use to fix certain voxels to background

%Prepare array for storing indexes after 1st iter
itIDX2=[];

for iter=1:nIter
    if (useMoment==1 && iter==2)
        lambda=lambda;
    end
    
    
    
    for k = itIDX
        if iter==1
            X1=pointLOR(k,1); Y1=pointLOR(k,2);
            X2=pointLOR(k,3); Y2=pointLOR(k,4);
            [weightS, indI, indJ, nV, S] =  projection_matrix_row_2D(X1, Y1, X2, Y2, imageSize, delta, tol);
            weightS=S*weightS;
            ind = sub2ind([imageSize(1) imageSize(2)],indI,indJ);% using different imSize(1) and imSize(2) Allows reconstruction with differing (x,y) dims
            
            
            if proj_data(k)==deadThresh
                nV=0;
            end
            if deadThresh < proj_data(k) &&  proj_data(k)<=thresh%if x-ray going through air but not a dead pixel
                nV=0;%make sure to skip on next iter
                binMask(ind)=1;
            end
            
        else
            weightS=geo(k).weightS;
            ind=geo(k).ind;
            nV=geo(k).nV;
            
            if iter==2%get rid of updates wrt. know air voxels(remain constant)
                weightS=weightS(binMask(ind)==0);
                ind=ind(binMask(ind)==0);
                geo(k).ind=ind;
                geo(k).weightS=weightS;
            end
        end
        
        if(nV>0)%Check that the geometry lines are inside the volume
            if(iter==1)
                %Is allways true after first iteration (with saving of geom)
                projValueIter=sum(weightS.*f_inp(ind));
                rowNorm=sum(weightS.*weightS);
                product=(proj_data(k)-projValueIter)/rowNorm; %ART
                
                vCount(ind)=vCount(ind)+1;
                %vSum(ind)=vSum(ind)+lambda*product*weightS;
                f_inp(ind) = f_inp(ind) +lambda*product*weightS; %ART
                
                itIDX2=[itIDX2; k];
                %save geometry for any iteration beyond 1
                geo(k).weightS=weightS;
                geo(k).ind=ind;
                geo(k).nV=nV;
            else
                
                %Is allways true after first iteration (with saving of geom)
                projValueIter=sum(weightS.*f_inp(ind));
                rowNorm=sum(weightS.*weightS);
                product=(proj_data(k)-projValueIter)/rowNorm; %ART
                
                %vMoment{k}=param.moment*vMoment{k}+lambda*product*weightS;
                %vCount(ind)=vCount(ind)+1;
                f_inp(ind) = f_inp(ind) +param.moment*vMomentF(ind)+lambda*product*weightS; %ART
                %vMomentF(ind)= vMomentF(ind)+(lambda*product*weightS)./vCount(ind);
            end
        end
    end %for k
    
    
    if iter==1
        
        itIDX=itIDX2';%on all next iteration dont need to wasted time on useless lines
        itIDX(itIDX>numel(geo))=[];
    end
    
    if rem(iter,10)==0
        itOut(:,:,round(iter/10))=f_inp;
    end
    
    [f_TV] = TV_2D(f_inp, imageSize, nTV, epsTV, alpha);
    f_inp(:,:)=f_TV(:,:);
    
    
    f_inp=max(0,f_inp);
    f_inp(binMask==1)=backVal;%make sure background mask is set
    
    vMomentF=(f_inp-f_old)./vCount;
    f_old=f_inp;
    
    
end %nIter

f_out=f_inp;
return

function [pointLOR, num_data] = par_beam_geometry(range_angle, num_ang, num_det, detect, source)
%
% a system Sources-Detectors rotates around the origin Oxy
% Input:
% range_angle in degrees
% detect - array (num_det x 2) of (x,y) initial coordinates of detectors
% source - array (num_det x 2) of (x,y) initial coor-s of the source
% detect and source differ in their Y-coordinates
% Output: pointLOR (num_det x 4) - (x,y) coordinates of point1 and
% point2 on the LOR, [x1 y1 y1 y2]
%
num_data=num_ang*num_det;
pointLOR=zeros(num_data,4);
ang_sample = (range_angle/180*pi)/num_ang;

for j=1:num_ang
    ang=ang_sample*(j-1);
    cosAng=cos(ang);
    sinAng=sin(ang);
    for i=1:num_det
        detX0=detect(i,1);
        detY0=detect(i,2);
        sourceX0 = source(i,1);
        sourceY0 = source(i,2);
        srcX= sourceX0*cosAng - sourceY0*sinAng;
        srcY= sourceX0*sinAng + sourceY0*cosAng;
        detX= detX0*cosAng - detY0*sinAng;
        detY= detX0*sinAng + detY0*cosAng;
        ind_data=(j-1)*num_det + i;
        pointLOR(ind_data,1)= srcX;
        pointLOR(ind_data,2)= srcY;
        pointLOR(ind_data,3)= detX;
        pointLOR(ind_data,4)= detY;
    end % for i
end %for j
%
return

function [weightS, indJ, indK, nV, S] = projection_matrix_row_2D(Y1, Z1, Y2, Z2, imageSize, delta, tol)
% Input: Points (Y1,Z1) and (Y2,Z2); 3D array imageSizeximageSizeximageSize
% Rimg - radius of sphere inscribed tightly into the 3D array
%
% calls: Points2Normal,  ray_tracing_along_Z
%
%indJ=zeros(1,2*imageSize);
%indK=zeros(1,2*imageSize); weightS=zeros(1,2*imageSize);
%
[~, Ind] = max([abs(Y2-Y1) abs(Z2-Z1)]);
% we change Y to Z, if Y is of maximal increment
% (1 2) --> (2 1)  Z -> Y, Y-> Z
if Ind == 1,
    Ynew1=Z1; Znew1=Y1;
    Ynew2=Z2; Znew2=Y2;
end
% (1 2) --> (1 2) nothing to change
if Ind == 2,
    Ynew1=Y1; Znew1=Z1;
    Ynew2=Y2; Znew2=Z2;
end
% by this construction, abs(Znew2-Znew1) is maximal
% |dZ| > |dY|
% we wish to have points with Znew2 > Znew1 otherwise exchange swap them
Ycur1=Ynew1; Zcur1=Znew1;
if Znew1 > Znew2,
    Ynew1= Ynew2; Znew1= Znew2;
    Ynew2= Ycur1; Znew2= Zcur1;
end
% now, Znew1 < Znew2 and proceed to normal representation of the lines
[pnYZ,cosYZ,sinYZ] = Points2Normal(Ynew1,Znew1, Ynew2, Znew2);
%
[weightS, indJ, indK, ~, S] = ray_tracing_along_Z_2D(pnYZ,cosYZ,sinYZ, max(imageSize), delta, tol);

% We need to reoreder the axes Y Z back into original order
% (1 2) -> (2 1 ) %Z -> Y, Y-> Z; %  K -> J, J -> K
if Ind == 1,
    temp=indJ;
    indJ=indK;indK=temp;
end
%subtract index offset to alloq for reconstructed volume centered on (imSize(1)/2, imSize(2)/2)
%the subtraction is needed due to ray_tracing assuming equal (x,y) dims
indK=indK-round((max(imageSize)-imageSize(2))/2);
indJ=indJ-round((max(imageSize)-imageSize(1))/2);

%get rid of indices outside image
outer=indK < 1 | indK > imageSize(2) | indJ<1 | indJ > imageSize(1);

indK(outer)=[];
indJ(outer)=[];
weightS(outer)=[];
nV=numel(indJ);


return

function [weightS, indJ, indK, nV, S] = ray_tracing_along_Z_2D(pnYZ,cosYZ,sinYZ, imageSize, delta, tol)
%
% pnYZ - normal distance of line
% pnYZ,cosYZ,sinYZ - normal parameters of the line's projection  on YZ plane
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
indJ=zeros(1,2*imageSize);
indK=zeros(1,2*imageSize);
weightS=zeros(1,2*imageSize);
%
R=imageSize*delta/2;
% initialization
Zplane=-R-delta;
% YZ plane
Syz=delta/abs(cosYZ);
STyz=Syz*sinYZ;
A1yz=(pnYZ-Zplane*sinYZ)/cosYZ;
%
S=Syz;
nV=0; % number of Visited voxels
for KsliceZ=1:imageSize
    Zplane=Zplane+delta;
    A1yz=A1yz-STyz;
    A2yz=A1yz-STyz;
    J1yz=floor((R+A1yz)/delta)+1; %  - along Y
    J2yz=floor((R+A2yz)/delta)+1;
    if J1yz == J2yz,
        nV=nV+1;
        indJ(nV)= J1yz; indK(nV)=KsliceZ; weightS(nV)=1;
        continue
    end
    %
    if J1yz ~= J2yz,
        Cyz=-R + min(J1yz,J2yz)*delta;
        STyzAbs=abs(STyz);
        if STyzAbs <= tol,
            nV=nV+1;
            indJ(nV)=J1yz;indK(nV)=KsliceZ;weightS(nV)=0.5;
            nV=nV+1;
            indJ(nV)=J2yz;indK(nV)=KsliceZ;weightS(nV)=0.5;
            continue
        end
        S1=abs(Cyz-A1yz)/STyzAbs;
        nV=nV+1;
        indJ(nV) = J1yz; indK(nV)=KsliceZ;  weightS(nV)=S1;
        nV=nV+1;
        indJ(nV) = J2yz; indK(nV)=KsliceZ; weightS(nV)=1-S1;
        continue
    end % if J1yz ~= J2yz
end % for KsliceZ=1:N
return

function [pnXZ,cosXZ,sinXZ] = Points2Normal(X1,Z1, X2, Z2)
% should be |Z2-Z1| > |X2-X1|, Z2 > Z1
%
alpha = atan((X1-X2)/(Z2-Z1));
% should be -45 <= alpha <= 45
cosXZ=cos(alpha);
sinXZ=sin(alpha);
pnXZ=X2*cosXZ+Z2*sinXZ;
return

function [f_inp] = TV_2D(f_inp, imageSize, nTV, epsTV, alphaTV)
%
% Total variation on 2D image
%
%f_neg_pos=zeros(imageSize,imageSize);
%deriv= zeros(imageSize,imageSize);
%f_inp=[f_inp(:,1) f_inp f_inp(:,end)];
%f_inp=[f_inp(1,:); f_inp; f_inp(end,:)];
%imageSize=imageSize+2;

f_neg_pos= max(0,f_inp);
beta=sqrt(sum((f_inp(:)-f_neg_pos(:)).^2));
%beta=norm2(f_TV-f_neg_pos,imageSize);
for iTV=1:nTV
    deriv = partDerivTVtwoD(f_inp, imageSize, epsTV);
    %deriv = partDerivTV_v2(f_inp, imageSize, epsTV);
    norm2deriv=sqrt(sum((deriv(:)).^2));
    deriv=deriv/norm2deriv;
    f_inp=f_inp - alphaTV*beta*deriv;
end %iTV
%f_inp=f_inp(2:end-1,2:end-1);
return

function [ f_out ] = partDerivTVtwoD(f,N,epsTV)
%
% partial derivative of TV of two-Dimensional array f of size NxN
% typicl eps=0.0000001; 10-6, 10-8
%
f_out=zeros(N(1),N(2));
id1=2:N(1)-1;%indexing central part of image(no edges)
id2=2:N(2)-1;%indexing central part of image(no edges)

sqrt1=realsqrt(epsTV+(f(id1,id2)-f(id1-1,id2)).^2 + (f(id1,id2)-f(id1,id2-1)).^2);
sqrt2=realsqrt(epsTV+(f(id1+1,id2)-f(id1,id2)).^2 + (f(id1+1,id2)-f(id1+1,id2-1)).^2);
sqrt3=realsqrt(epsTV+(f(id1,id2+1)-f(id1,id2)).^2 + (f(id1,id2+1)-f(id1-1,id2+1)).^2);

f_out(id1,id2) =(f(id1,id2)-f(id1-1,id2)+f(id1,id2)-f(id1,id2-1))./sqrt1 - (f(id1+1,id2)-f(id1,id2))./sqrt2 - (f(id1,id2+1)-f(id1,id2))./sqrt3;


return


