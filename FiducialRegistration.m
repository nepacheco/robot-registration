%% Fiducial Registration
clc;clear; close all;

%% Generate Ground Truth Fiducial Locations and Proper Frame Transformations
% All the transforms and locations here assume we don't have any fiducial
% localization error (FLE = 0)

N = 7; % Number of fiducials 
rng(1)
T_R_W = [rotz(120)*roty(50) [13 -20 3]'; 0 0 0 1]; % This transformation shows the proper location of the world in the robot frame
Fid_W = [(rand(3,N) - 0.5)*40]; % These are the exact locations of the fiducials in the world frame;
Fid_R = T_R_W(1:3,1:3)*Fid_W + T_R_W(1:3,4); % These are the exact locations of the fiducials in the robot frame
[R,t,FRE] = point_register(Fid_W,Fid_R);
FRE2 = FRE^2;

fprintf("Assuming the FLE = 0, we see that the resulting FRE2 = %0.4f\n",FRE2);

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
numSamples = 1000;
meanFRE2 = 0;
meanFLE2 = 0;
FidR_Vec = zeros(3,N*numSamples);
FidW_Vec = zeros(3,N*numSamples);
for i = 1:numSamples
    sdR = [0.2;0.2;0.2]; % x y and z directions in the robot frame
    sdW = [0.5;0.5;0.5]; % x y and z directions in the world frame
    Fid_W_noise = Fid_W + normrnd(0,sdW*ones(1,N));
    Fid_R_noise = Fid_R + normrnd(0,sdR*ones(1,N));
    FidR_Vec(:,(i-1)*N+1:i*N) = Fid_R_noise;
    FidW_Vec(:,(i-1)*N+1:i*N) = Fid_W_noise;
    
    nFLE2 = trace((Fid_R_noise - Fid_R)'*(Fid_R_noise-Fid_R))/N;
    nFLE2 = nFLE2 + trace((Fid_W_noise - Fid_W)'*(Fid_W_noise-Fid_W))/N;
    meanFLE2 = (meanFLE2*(i-1) + nFLE2)/i;

    [R,t,FRE] = point_register(Fid_W_noise,Fid_R_noise);
    FRE2 = FRE^2;
    
    meanFRE2 = (meanFRE2*(i-1) + FRE2)/i;
end
FLE2 = N/(N-2) * meanFRE2;

fprintf("\nAfter %d samples, with a world SD = %0.2f and a robot SD = %0.2f\n" + ...
    "the <FRE^2> = %0.3f\n",numSamples,sdW,sdR,meanFRE2);

fprintf("\nUsing <FRE^2> to compute <FLE^2> gives us a value of %0.3f\n",FLE2)
fprintf("Experimentally we find <FLE^2> to be %0.3f\n",meanFLE2);
fprintf("Theoretically the <FLE^2> = %0.3f based on our noise\n",sdR'*sdR + sdW'*sdW)

figure(1)
clf;
hold on;
FidW_Vec_Trans = T_R_W(1:3,1:3)*FidW_Vec + T_R_W(1:3,4);
scatter3(FidR_Vec(1,:),FidR_Vec(2,:),FidR_Vec(3,:),'.','DisplayName',"Robot Fiducials")
scatter3(FidW_Vec_Trans(1,:),FidW_Vec_Trans(2,:),FidW_Vec_Trans(3,:),'.','DisplayName',"World Fiducials")
plotTransforms(zeros(1,3),[1 zeros(1,3)],'FrameSize',10)
plotTransforms(T_R_W(1:3,4)',rotm2quat(T_R_W(1:3,1:3)),'FrameSize',10)
hold off;
grid on;
xlabel("X Axis")
ylabel("Y Axis")
zlabel("Z Axis")
legend('Location','best')
title("Fiducial Registration")
axis equal

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