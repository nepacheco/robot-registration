%% Fiducial Registration
clc;clear; close all;

%% Generate Ground Truth Fiducial Locations and Proper Frame Transformations
% All the transforms and locations here assume we don't have any fiducial
% localization error (FLE = 0)

N = 7; % Number of fiducials 
rng(2)
T_R_W = [rotz(120)*roty(50) [13 -20 3]'; 0 0 0 1]; % This transformation shows the proper location of the world in the robot frame
Fid_W = [(rand(3,N) - 0.5)*20]; % These are the exact locations of the fiducials in the world frame;
Fid_R = T_R_W(1:3,1:3)*Fid_W + T_R_W(1:3,4); % These are the exact locations of the fiducials in the robot frame
[R,t,FRE] = point_register(Fid_W,Fid_R);
FRE2 = FRE^2;
fprintf("Assuming the FLE = 0, we see that the resulting FRE2 = %0.4f\n",FRE2);

%% Camera setup
checkRows =8;
checkCols = 11;
checkPixSize = 35; % [pixels]
checkX = 0.5:0.5:(checkCols)*0.5; % [cm]
checkY = (checkRows)*0.5:-0.5:0.5; % [cm]
[checkX,checkY] = meshgrid(checkX,checkY);
checkX = reshape(checkX,[1,checkRows*checkCols]);
checkY = reshape(checkY,[1,checkRows*checkCols]);

% Assumed checkerboard locations in the world frame
pix_W_assumed = [checkX;checkY;zeros(1,checkRows*checkCols)]; 

 % Actual location of the checkerboard in the world frame
pix_W = rotz(0)*roty(0)*rotz(0)*pix_W_assumed + [0;0;0];

figImages = figure(1);
tiledlayout(figImages,'flow')
% The perspective transform first rotates the image by 180 degrees, then
% centeres it in the x axis and finally adds a tilt into the page by 45
% degrees
theta = deg2rad(0);

% To build the homography matrix, we multiple two matrices.
%   1. First flip about the x axis, center the image on the y axis, and
%       tilt into the z axis
%   2. Shift the frame so that x and y axis start at corner of the image. 
A = [cos(theta), -sin(theta), 0.00;
     sin(theta),  -cos(theta), 1/((checkRows+1)*checkPixSize)*(2/sqrt(2)-1); % the negative on the cosine is for a verticle flip
     -(checkCols+1)*checkPixSize/2,  0, 1]*[1 0 0; 0 1 0; (checkCols+1)*checkPixSize/2, (checkRows+1)*checkPixSize*sqrt(2)/2 1];
tform = projective2d(A); % Because I am in R2021a I have to use projective2d(A)
% which uses the post multiply convection [x y z]^T = [u v 1]^T * A

% Creating the actual checkerboard image is very cumbersome with MATLAB's
% checkerboard function
checker = checkerboard(checkPixSize,round(checkRows/2)+1,round(checkCols/2)+1);
% checker = checker((checkPixSize+1):end,1:end-checkPixSize);
checker = checker(1:end-checkPixSize,1:end-checkPixSize);
if mod(checkRows,2) == 1
    checker = checker(1:end-checkPixSize,:);
end
if mod(checkCols,2) == 1
    checker = checker(:,1:end-checkPixSize);
end
checker(1:checkPixSize,1:checkPixSize) = 0.5;
nexttile()
imshow(checker);
set(gca,'YDir','normal')
title("Checkerboard in World Frame")
xlabel("X axis")
ylabel("Y Axis")
axis on
% Warp the flat checkerboard to into the tilted position. Then detect the
% corners of the checkerboard to define the "proper" location of the
% checkerboard in the camera frame
[warpChecker, RB] = imwarp(checker,tform);
nexttile()
imshow(warpChecker)
axis on
title("Perfect Checkerboard as seen in Camera")


% Apply noise and filtering to the perfect camera image and then identify
% the corners of the checkerboard
noiseWarpChecker = imnoise(warpChecker,"gaussian",0,0.1); % Add Gaussian noise
noiseWarpChecker = imgaussfilt(noiseWarpChecker,2); % Apply a blurring to the image
[imagePoints,boardSize] = detectCheckerboardPoints(noiseWarpChecker);
markedChecker = insertMarker(noiseWarpChecker, imagePoints, 'o', 'MarkerColor', 'red', 'Size', 2);
nexttile()
imshow(markedChecker)
axis on
title("Noise added to Camera Image")

% Apply a projective transform between the identified checkerboard points
% and the location of those points in the world and 'dewarp' the image
tform3 = fitgeotrans(imagePoints,checkPixSize*2*pix_W_assumed(1:2,:)','projective');
[warpChecker3,RB3] = imwarp(noiseWarpChecker,tform3);
nexttile()
imshow(warpChecker3)
axis on
axis equal
title("Warped with fitted transform from found checkerboard points")

% Calculate the pixel reprojection error based on the transform from the
% noisy image
warpedWorld = [imagePoints,ones(checkRows*checkCols,1)]*tform3.T;
warpedWorld = warpedWorld./warpedWorld(:,3);
error = warpedWorld-[checkPixSize*2*pix_W_assumed(1:2,:)',ones(checkRows*checkCols,1)];
mse = mean(sum(error.*error,2));
fprintf("\nPixel Reprojection Error: %0.3f pixels\n",sqrt(mse));

%% Introduce FLE 
% We will now introduce Fiducicial localization error in both the world
% frame and the robot frame by adding zero mean gaussian noise to the
% fiducial positions in the world and in the robot frames. 

rng(1)
sdR = 0.1;
sdW = 0.1;
Fid_W_noise = Fid_W + normrnd(0,sdW,3,N);
Fid_R_noise = Fid_R + normrnd(0,sdR,3,N);

[R,t,FRE] = point_register(Fid_W_noise,Fid_R_noise);
FRE2 = FRE^2;

fprintf("\nAssuming there is a standard deviation on the world fiducials\n" + ...
"and robot fiducials of %0.2f and %0.2f respectively, we get\n" + ...
"an FRE^2 of %0.2f in a singular instance\n",sdW,sdR,FRE2);

%% FLE population statistics
% Once we begin to use randomness, we must create an ensemble of data from
% which we can extract means and variances. If we remove the seed selection
% from the previous section, we would see we get a different FRE each time.
% Therefore, we must perform the calculation numerous times to get a mean
% FRE. 
numSamples = 100000;
msFREex = 0;
msTREex = 0;
msFLEex = 0;
mFLE = 0;
FidR_Vec = zeros(3,N*numSamples);
FidW_Vec = zeros(3,N*numSamples);
targetR_vec = zeros(3,numSamples);
sdR = [0.2;0.2;0.2]; % S.D. x y and z directions in the robot frame
sdW = [0.5;0.5;0.5]; % S.D. x y and z directions in the world frame
for i = 1:numSamples
    % Generate the fiducial in the world after FLE is applied
    Fid_W_noise = Fid_W + normrnd(0,sdW*ones(1,N));
    % Generate fiducial location in robot frame after FLE is applied
    Fid_R_noise = Fid_R + normrnd(0,sdR*ones(1,N));
    FidR_Vec(:,(i-1)*N+1:i*N) = Fid_R_noise;
    FidW_Vec(:,(i-1)*N+1:i*N) = Fid_W_noise;
    
    % Calculate an experimental FLE
    nFLE2 = trace((Fid_R_noise - Fid_R)'*(Fid_R_noise-Fid_R))/N;
    nFLE2 = nFLE2 + trace((Fid_W_noise - Fid_W)'*(Fid_W_noise-Fid_W))/N;
    msFLEex = (msFLEex*(i-1) + nFLE2)/i;
    mFLE = (mFLE*(i-1) + sqrt(nFLE2))/i;

    % Perform rigid registration between fiducial pairs
    [R,t,FRE] = point_register(Fid_W_noise,Fid_R_noise);
    FRE2 = FRE^2;
    
    msFREex = (msFREex*(i-1) + FRE2)/i;
    
    % Choose a target in the world frame
    target_W = [5,4,10]';
    % Get the exact location of target in the robot frame
    target_R = T_R_W*[target_W;1];
    % Get the location of the target in the world frame after applying the
    % transform from the point registration
    target_R_noise = R*target_W + t;
    targetR_vec(:,i) = target_R_noise;
    % Calculate the TRE^2
    TRE2 = sum((target_R(1:3) - target_R_noise).^2);
    msTREex = (msTREex*(i-1) + TRE2)/i;

    % Get the world location in pixels
    target_Cam = [checkPixSize*2*target_W(1:2)',1]*tform.T; % this is where the pixel should be
    target_Cam_noise = [checkPixSize*2*target_W(1:2)',1]*inv(tform3.T);
    
end
msFLE = N/(N-2) * msFREex;

fprintf("\nAfter %d samples, with a world SD = [%0.2f,%0.2f,%0.2f] and a robot SD = [%0.2f,%0.2f,%0.2f]\n" + ...
    "the <FRE^2> = %0.3f\n",numSamples,sdW,sdR,msFREex);

fprintf("\nUsing <FRE^2> to compute <FLE^2> gives us a value of %0.3f\n",msFLE)
fprintf("Experimentally we find <FLE^2> to be %0.3f\n",msFLEex);
fprintf("Theoretically the <FLE^2> = %0.3f based on our noise\n",sdR'*sdR + sdW'*sdW)

fprintf("\nTarget Located at (%0.2f,%0.2f,%0.2f)\n",target_W(1:3));
fprintf("Experimentally, we find <TRE^2> = %0.3f\n",msTREex)
fprintf("Theoretically, <TRE^2> = %0.3f\n",treapprox(Fid_R,target_R(1:3),sqrt(msFLEex))^2)


%% Plotting Entire Setup
numPlotSamples = min(100,numSamples);
transparency = 0.2;
figure(2)
clf;
hold on;
T_W_R = inv(T_R_W); % location of robot in world's reference frame
FidR_Vec_Trans = T_W_R(1:3,1:3)*FidR_Vec + T_W_R(1:3,4);
% Fiducials Samples in Robot Frame
s2 = scatter3(FidR_Vec_Trans(1,1:numPlotSamples*N),...
    FidR_Vec_Trans(2,1:numPlotSamples*N),FidR_Vec_Trans(3,1:numPlotSamples*N),...
    15,'ro','filled',...
    'DisplayName',"Robot Fiducials",'MarkerEdgeAlpha',transparency,...
    'MarkerFaceAlpha',transparency);
% Fiducical Samples in World Frame
s3 = scatter3(FidW_Vec(1,1:numPlotSamples*N),FidW_Vec(2,1:numPlotSamples*N),...
    FidW_Vec(3,1:numPlotSamples*N),15,'bo','filled',...
    'DisplayName',"World Fiducials",'MarkerEdgeAlpha',transparency,...
    'MarkerFaceAlpha',transparency);
% Target in world frame
s4 = scatter3(target_W(1),target_W(2),target_W(3),100,'bx','DisplayName',"Target",...
    'LineWidth',2,'MarkerEdgeAlpha',1.0);

% Target samples after registration
target_W_noise = T_W_R(1:3,1:3)*targetR_vec + T_W_R(1:3,4);
s5 = scatter3(target_W_noise(1,1:numPlotSamples),target_W_noise(2,1:numPlotSamples),...
    target_W_noise(3,1:numPlotSamples),15,'ro','filled',...
    'DisplayName',"Target Samples",'MarkerEdgeAlpha',transparency,...
    'MarkerFaceAlpha',transparency);

tf = plotTransforms(zeros(1,3),[1 zeros(1,3)],'FrameSize',10);
plotTransforms(T_W_R(1:3,4)',rotm2quat(T_W_R(1:3,1:3)),'FrameSize',10)
% Plot Checkerboard in actual world frame
scatter3(pix_W(1,:),pix_W(2,:),pix_W(3,:),'k.','DisplayName',"Checkerboard Corners")
% Plot the Location of the assumed checkerboard
scatter3(pix_W_assumed(1,:),pix_W_assumed(2,:),pix_W_assumed(3,:),'r.','DisplayName',"Assumed Checkerboard Corners")
hold off;
grid on;
xlabel("X Axis (cm)")
ylabel("Y Axis (cm)")
zlabel("Z Axis (cm)")
legend('Location','best')
title("Fiducial Registration")
axis equal
set(gca,'FontSize',18)

%% Functions

function [Cov_TRE,Cov_FRE] = TRE(X,W,Cov_FLE,r)
T = W(:,1,:);U = W(:,2,:);V = W(:,3,:);
Xr = reshape([repmat(X(1,:),3,1);repmat(X(2,:),3,1);repmat(X(3,:),3,1)],size(W));
X1 = Xr(:,1,:);X2 = Xr(:,2,:);X3 = Xr(:,3,:);
C = reshape(permute([-U.*X3+V.*X2,T.*X3-V.*X1,-T.*X2+U.*X1,T,U,V],[1,3,2]),[],6);
W_cell = num2cell(W,[1 2]);
Cov_cell = num2cell(Cov_FLE,[1 2]);
W_Cov_Wt = blkdiag(W_cell{:})*blkdiag(Cov_cell{:})*blkdiag(W_cell{:})';
D = [[0,r(3),-r(2);-r(3),0,r(1);r(2),-r(1),0],eye(3)];
pinvC = pinv(C);
Cov_TRE = D*pinvC*W_Cov_Wt*pinvC'*D';
Cov_FRE = (eye(size(C,1))-C*pinvC)*W_Cov_Wt*(eye(size(C,1))-C*pinvC);
end