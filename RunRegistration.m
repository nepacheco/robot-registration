  clc; clear; close all;

%% Setting initial transforms
probe = [0;0;10]; % [cm] positition of the probe relative to the end effector
probeLength = norm(probe); % [cm]
fiducials = [0.5 0 0;
             0 0 0;
             -0.5 -0.5 0.5;
             -0.5 0.5 0]'*100; % [cm] location of the fiducials in the world frame
TrobotInWorld = [rotz(50)*roty(30),[50;10;-13];
                 zeros(1,3), 1]; % [cm] position of robot base in world frame

noiseMag = 0; % [cm] The amount of noise from the fiducial and probe

%% Determine probe position for each fiducial point and calculate fiducial in robot frame
numReadings = 5;
figure
triad('scale',25,'matrix',eye(4,4),'linewidth',2.0);
triad('scale',25,'matrix',TrobotInWorld,'linewidth',2.0);
hold on
posFidInRobotBase = zeros(3,size(fiducials,2));
for i = 1:size(fiducials,2)
    fid = fiducials(:,i);
    plot3(fid(1),fid(2),fid(3),'Marker','x','Markersize',15,'Color','r')

    robotTip = zeros(4,numReadings);
    for n = 1:numReadings
        theta = rand(1,1)*pi/2; % Limit slant to 0-90 degrees
        z = cos(theta)*probeLength;
        h = sin(theta)*probeLength;
        phi = rand(1,1)*2*pi; % rotation about the fiducial point
        x = cos(phi)*h;
        y = sin(phi)*h;
        
        % In reality we won't know the tip position in the world, but
        % because this is simulation we generate this value (with some
        % noise) and then convert it to the tip position in the robot base
        % frame, which we will have when performing this on a robot.
        TipInWorld = [(fid + [x;y;z] + (rand(3,1)*2-1)*noiseMag);1];
        plot3(TipInWorld(1),TipInWorld(2),TipInWorld(3),'Marker','.','Markersize',10,'Color','b')
        % Get the position of the robot tip in the base frame based on the
        % orientation of the probe, the current fiducial, and the world to
        % robot transform
        robotTip(:,n) = TrobotInWorld\TipInWorld;
    end
    % Determine the center point of the sphere that the generated tip
    % positions lie on. 
    M = robotTip';
    Y = -(robotTip(1,:)'.^2 + robotTip(2,:)'.^2 + robotTip(3,:)'.^2);
    X = (M'*M)\M'*Y;  % Moore-Penrose inverse for linear regression
    centX = -1/2*X(1);
    centY = -1/2*X(2);
    centZ = -1/2*X(3);
    radius = sqrt(centX^2 + centY^2 + centZ^2 - X(4));
    % The position of the fiducial in the robot's base frame is given by
    % that center position
    posFidInRobotBase(:,i) = [centX;centY;centZ];
end
view (45,20)
axis equal
grid on
%% Calculate the world to robot transform
% This is done by using the position of the fidiucials in the robot's frame
% and the world frame.
[R,t,FRE] = point_register(posFidInRobotBase,fiducials)

triad('scale',25,'matrix',[R,t;0 0 0 1],'linewidth',1.0,'tag','Calculated Transform');


