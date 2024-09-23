clc; clear; close all;
nPoints = 10000;
Fiducials = [5.050711735552091408e-01 3.515693879346971640e-01 2.762407860116669611e-01 2.491718530947737986e-01 2.727334088413887470e-01 3.738709586037104460e-01 4.506034943387034475e-01;
-1.517067204163148331e-01 -1.282513508869667085e-01 -1.296496075800286918e-01 -2.765234465102743913e-02 1.238562332573615588e-01 1.507688123719878037e-01 1.262056640464310486e-01;
1.311176343923934340e-01 6.854201717305391517e-02 8.209621687437256010e-02 1.074948859123116274e-01 1.050382139117859992e-01 1.278565522353995454e-01 1.301764393903707406e-01];


FLE = sqrt(1.47);
centroid = mean(Fiducials,2);
rmseBreakpoints = [1,1.5,2];
rmseTargets = cell(1,length(rmseBreakpoints));
for i = 1:nPoints
    direction = rand(3,1)-0.5;
    direction = direction./norm(direction);
    magnitude = rand(1,1)*0.5;
    target = direction*magnitude;
    TRE_RMS = treapprox(Fiducials,centroid + target,FLE);
    for r = 1:length(rmseBreakpoints)
        if TRE_RMS < rmseBreakpoints(r)
            rmseTargets{r} = [rmseTargets{r} target];
            break;
        end
    end
end

Tv = treinv(Fiducials,1,FLE);
[X,Y,Z] = ellipsoid(0,0,0,norm(Tv(:,1)),norm(Tv(:,2)),norm(Tv(:,3)));
rotMatrix = [Tv(:,1)./norm(Tv(:,1)),Tv(:,2)./norm(Tv(:,2)),-Tv(:,3)./norm(Tv(:,3))];
[ax] = rotm2axang(rotMatrix);
% Tv = abs(Tv);

figure(1);
clf;
hold on;
for r = 1:length(rmseBreakpoints)
    scatter3(rmseTargets{r}(1,:),rmseTargets{r}(2,:),rmseTargets{r}(3,:),25,'filled',...
        'DisplayName',sprintf("RMSE: %0.1f",rmseBreakpoints(r)),...
        'MarkerFaceAlpha',0.4,'MarkerEdgeAlpha',1/r);
end
for i = 1:3
    x = [-Tv(1,i),Tv(1,i)];%+centroid(1);
    y = [-Tv(2,i),Tv(2,i)];%+centroid(2);
    z = [-Tv(3,i),Tv(3,i)];%+centroid(3);
    plot3(x,y,z,'LineWidth',4)
end
s = surf(X,Y,Z);
rotate(s,ax(1:3),rad2deg(ax(4)))
legend()
xlabel("X Axis (m)");
ylabel("Y Axis (m)");
zlabel("Z Axis (m)");
axis equal
grid on;