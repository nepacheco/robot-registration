clc; clear; close all;
%% Ground Truth
rng(1)
N = 7;
TRW = [rotz(130)*roty(30) [10 -23 5]'; zeros(1,3) 1];
fidW = (rand(3,N)- 1/2)*40;
fidR = TRW(1:3,1:3)*fidW + TRW(1:3,4);

%% Example cases
stdR = 0.45075263;
stdW = 0.45075263;
numSamples = 10000;
fidWVec = zeros(3,N*numSamples);
fidRVec = zeros(3,N*numSamples);
meanFRE = 0;
meanFRE2 = 0;

for i = 1:numSamples
    fidW_ = fidW + normrnd(0,stdW,3,N);
    fidR_ = fidR + normrnd(0,stdR,3,N);
    fidWVec(:,(i-1)*N+1:(i*N)) = fidW_;
    fidRVec(:,(i-1)*N+1:(i*N)) = fidR_;
    
    [R,t,FRE] = point_register(fidW_,fidR_);
    meanFRE = (meanFRE*(i-1) + FRE)/i;
    meanFRE2 = (meanFRE2*(i-1) + FRE^2)/i;
end
FLE2 = N/(N-2) * meanFRE2;

fig = figure(1);
clf;
hold on
scatter3(fidRVec(1,:),fidRVec(2,:),fidRVec(3,:),'.',"DisplayName","Robot Fiducials")
fidWVec_trans = TRW(1:3,1:3)*fidWVec + TRW(1:3,4);
scatter3(fidWVec_trans(1,:),fidWVec_trans(2,:),fidWVec_trans(3,:),'.',"DisplayName","World Fiducials")
plotTransforms(zeros(1,3),[1 0 0 0],'FrameSize',10)
plotTransforms(TRW(1:3,4)',rotm2quat(TRW(1:3,1:3)),'FrameSize',10)
hold off
grid on;
axis equal
xlabel("X Axis")
ylabel("Y Axis")
zlabel("Z Axis");