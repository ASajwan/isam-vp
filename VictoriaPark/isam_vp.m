%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on GTSAM Copyright 2010, Georgia Tech Research Corporation,
% Atlanta, Georgia 30332-0415
% All Rights Reserved
% Authors: Baixu Chen, Anne Gu, Jason Hoving, Ashish Sajwan
%
% See LICENSE for the license information
%
% @brief Read Robotics Institute range-only Plaza2 dataset and do iSAM
% @author Frank Dellaert
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% preliminaries
clear
nSteps = 5000;
import gtsam.*
addpath('C:\Program Files\MATLAB\R2016b\toolbox\gtsam_toolbox')
addpath('.\vicpark')
makeVideo = 0;

global AAr;
AAr = [0:360]*pi/360;

Data = load_vp_si();
M = min(nSteps, length(Data.Control.time)); % number of motion commands
K = length(Data.Laser.time); %Num of measurements
minK = 20; % first threshold for measurement update
incK = 1; % threshold for incremental measurement updates
batchInitialization = true; % used for first measurement update

%% Initalize Params
global Param

% vehicle geometry
Param.a = 3.78; % [m]
Param.b = 0.50; % [m]
Param.L = 2.83; % [m]
Param.H = 0.76; % [m]

% 2x2 process noise on control input
sigma.vc = 0.02; % [m/s]
sigma.alpha = 2*pi/180; % [rad]
Param.Qu = diag([sigma.vc, sigma.alpha].^2);

sigma.x = 0.1; % [m] 
sigma.y = 0.1; % [m]
sigma.phi = 0.5*pi/180; % [rad]
Param.Qf = diag([sigma.x, sigma.y, sigma.phi].^2);

% 2x2 observation noise
sigma.r = 0.05; % [m]
sigma.beta = 1*pi/180; % [rad]
Param.R = diag([sigma.r, sigma.beta].^2);

%% Set Noise parameters
noiseModels.prior = noiseModel.Diagonal.Sigmas([.01 .01 deg2rad(0.5)]');    % [x y theta]
noiseModels.observ = noiseModel.Diagonal.Sigmas([deg2rad(1), .05]');        % [bearing range]

%% Initialize iSAM
isam = ISAM2;

%% Add prior on first pose
pose0 = Pose2(Data.Gps.x(2),Data.Gps.y(2),36*pi/180);
newFactors = NonlinearFactorGraph;  % used for batch factor update
wholeGraph = NonlinearFactorGraph;  % used for data association

newFactors.add(PriorFactorPose2(0,pose0,noiseModels.prior));
wholeGraph.add(PriorFactorPose2(0,pose0,noiseModels.prior));   

wholeValues = Values;
wholeValues.insert(0,pose0);  
initial = Values;
initial.insert(0,pose0);
odo = Values;
odo.insert(0,pose0);
landmarkEstimates = Values;

%% configure video
if makeVideo
    try
        votype = 'avifile';
        vo = avifile('video.avi', 'fps', min(5, 1/pauseLen));
    catch
        votype = 'VideoWriter';
        vo = VideoWriter('video', 'Motion JPEG AVI');
        v.Quality = 100;
        open(vo);
    end
end

%% iSAM
% initialize counters and one-time-update conditions
k = 1;
countK = 0;
update = false;
LandMarkInit = true;
t = min(Data.Laser.time(1),Data.Control.time(1));
LandMarkCount = 0;

% initialize pose and odometry
lastPose = pose0;
odoPose = pose0;
MotionModel = createMotionModel;
JacobianModel = createMotionJacobian;

% set threshold for creating new landmarks
threshold = 6^2;

% loop through control data
for i=1:nSteps

    % get odometry measurement
    dt = Data.Control.time(i) - t;
    t = Data.Control.time(i);
    alpha = minimizedAngle(Data.Control.alpha(i));
    vc = Data.Control.ve(i) / (1- Param.H/ Param.L * tan(alpha));

    % convert dead-reckoning control to relative motion
    [delta_x, delta_y, delta_heading] = control2odom(MotionModel, lastPose, alpha, vc,dt);

    % add odometry factor
    odometry = Pose2(delta_x,delta_y,minimizedAngle(delta_heading));

    % predict pose and update odometry (for visuals)
    predictedOdo = odoPose.compose(odometry);
    odoPose = predictedOdo;
    odo.insert(i,predictedOdo);
    
    % predict pose and integrate to pose and landmark estimate
    predictedPose = lastPose.compose(odometry);
    lastPose = predictedPose;
    initial.insert(i,predictedPose);
    wholeValues.insert(i,predictedPose); 
    landmarkEstimates.insert(i,predictedPose);
  
    % compute nonlinear motion noise
    J = MotionJacobian(JacobianModel, vc,alpha ,minimizedAngle(predictedPose.theta),dt);

    % error checking
    if det(J* Param.Qu * J') < 0 && abs( det(J* Param.Qu * J') ) > det(Param.Qf)
        disp('strange')
        odometry_noise = Param.Qf;
    else
        odometry_noise = J* Param.Qu * J' + Param.Qf;
    end

    % add motion factor between previous and current pose
    wholeGraph.add(BetweenFactorPose2(i-1, i, odometry, ...
        noiseModel.Gaussian.Covariance(odometry_noise)));
    newFactors.add(BetweenFactorPose2(i-1, i, odometry, ...
        noiseModel.Gaussian.Covariance(odometry_noise)));

    % loop over measurements until current control time
    while k<=K && t>=Data.Laser.time(k)

        % obtain measurements from dataset
        measureMent = detectTreesI16(Data.Laser.ranges(k,:));
        predLandMarkBatch = zeros(size(measureMent,2),2);
        
        % skip over empty measurements
        if size(measureMent,2) == 0
            k = k+1;
            continue
        end
        
        % predict landmark location for each measurment
        for j = 1:size(measureMent,2)

            range = measureMent(1,j);
            bearing = measureMent(2,j) -pi/2;
            
            % use factor graph to predict landmark location
            shift_rotation = Pose2(0,0,minimizedAngle(bearing));
            moving_by_range = Pose2(range,0,0);
            Newpose = lastPose.compose(shift_rotation).compose(moving_by_range);
            predLandMark = Point2(Newpose.x, Newpose.y);
            predLandMarkBatch(j,:) = [predLandMark.x, predLandMark.y];
            
        end
        
        % initialize landmarks for all measurements on first time step
        if LandMarkInit
            [BatchIndex, BatchExist] = NN(wholeValues, LandMarkCount,predLandMarkBatch, threshold);
            LandMarkInit = false;
            
        % perform data association of choice
        else
        % Euclidean NN without constraints
        %[BatchIndex,BatchExist] = NNBatch_woCons(predLandMarkBatch,wholeValues, LandMarkCount, threshold);
        
        % Mahalanobis NN without constraints
        [BatchIndex,BatchExist] = MLBatch_woCons(i, lastPose,measureMent(1:2,:), wholeGraph, wholeValues, LandMarkCount, threshold);
        
        % Euclidean NN with constraints
        %[BatchIndex, BatchExist] = NNBatch(predLandMarkBatch,wholeValues, LandMarkCount, threshold);
        
        % Mahalanobis NN with constraints
        %[BatchIndex,BatchExist] = MLBatch(i, lastPose,measureMent(1:2,:), wholeGraph, wholeValues, LandMarkCount, threshold);
        end
        
        % iterate through measurement-landmark pairs
        for landIndex = 1:length(BatchIndex)
            
            exist = BatchExist(landIndex);
            nearest = BatchIndex(landIndex);
            predLandMark = Point2(predLandMarkBatch(landIndex,:)');
            
            % if the data association failed - initialize new landmark
            if ~exist
                LandMarkCount = LandMarkCount + 1;
                nearestKey =  symbol('L',LandMarkCount);
                initial.insert(nearestKey,predLandMark);
                landmarkEstimates.insert(nearestKey,predLandMark);
                wholeValues.insert(nearestKey,predLandMark);
            else
            nearestKey = symbol('L',nearest);
            end
            
            % create measurement factor between pose and landmark
            factor = BearingRangeFactor2D(i, nearestKey, Rot2(minimizedAngle(measureMent(2,landIndex)-pi/2)), measureMent(1,landIndex), noiseModels.observ);

            % add factor to batch update queue and data association graph
            newFactors.add(factor);
            wholeGraph.add(factor);
            
        end
        
        % update counters
        k = k+1;            % index for measurement data
        countK=countK+1;    % minimum # of measurements before batch update
        update = true;      % only update iSAM2 after first measurement
        
    end

    % Check whether to update iSAM 2 (batch update)
    if update && k > minK && countK > incK
        
        % only configure prior+pose the first time
        if batchInitialization
          batchOptimizer = LevenbergMarquardtOptimizer(newFactors, initial);
          initial = batchOptimizer.optimize();
          batchInitialization = false; % only once
        end

        % add new motion and measurement factors and calculate trajectory
        isam.update(newFactors, initial); 
        result = isam.calculateEstimate();
        
        % extract current pose
        lastPose = result.at(i);

        % relinearize and calculate new landmark estimates
        landmarkEstimates = Values;
        for j = 1:LandMarkCount
            key = symbol('L',j);
            landmarkEstimates.insert(key,result.at(key));
        end

        % relinearize and calculate new trajectory and map
        newFactors = NonlinearFactorGraph;
        wholeValues = Values;
        wholeValues.insert(result);
        initial = Values;
        
        % reset counter for batch update
        countK = 0;
        
    end
    
    % visuals - update graph every n steps
    if(mod(i,200) == 0)
        figure(1);clf;hold on;

        % extract and plot odometry
        XYT = utilities.extractPose2(odo);
        plot(XYT(:,1),XYT(:,2),'b-');

        % plot ground truth (until current time)
        gti = dsearchn(Data.Gps.time,Data.Control.time(i));
        plot(Data.Gps.x(1:gti),Data.Gps.y(1:gti),'g+','MarkerSize',3);

        % calculate, extract, and plot trajectory
        result = isam.calculateEstimate();
        XYT = utilities.extractPose2(result);
        plot(XYT(:,1),XYT(:,2),'r-','LineWidth',2);
        XY = utilities.extractPoint2(result);
        plot(XY(:,1),XY(:,2),'k.','MarkerSize',8);
        
        legend('odometry','ground truth','trajectory','landmarks');
        axis equal
        xlim([-180 150]);
        ylim([-200 150]);
        drawnow
        
        % for recording
        if makeVideo
            F = getframe(gcf);
            switch votype
              case 'avifile'
                vo = addframe(vo, F);
              case 'VideoWriter'
                writeVideo(vo, F);
              otherwise
                error('unrecognized votype');
            end
        end
    end
    
end

%% export video
if makeVideo
    fprintf('Writing video...');
    switch votype
      case 'avifile'
        vo = close(vo);
      case 'VideoWriter'
        close(vo);
      otherwise
        error('unrecognized votype');
    end
    fprintf('done\n');
end