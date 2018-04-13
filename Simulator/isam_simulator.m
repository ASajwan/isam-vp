%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GTSAM library - Copyright 2010, Georgia Tech Research Corporation,
% Atlanta, Georgia 30332-0415
% 
% Simulator - Copyright 2009, Ryan Eustice, University of Michigan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% preliminaries
clear
import gtsam.*
addpath('./simulator');

%% Initalize Params

global Param;
Param.initialStateMean = [180 50 0]';

% max number of landmark observations per timestep
Param.maxObs = 2;

% number of landmarks per sideline of field (minimum is 3)
Param.nLandmarksPerSide = 4;

% Motion noise (in odometry space, see p.134 in book).
Param.alphas = [0.05 0.001 0.05 0.01].^2; % std of noise proportional to alphas

% Standard deviation of Gaussian sensor noise (independent of distance)
Param.beta = [10, deg2rad(10)]; % [cm, rad]
Param.R = diag(Param.beta.^2);

sigmaInitial = 100;% draw initial landmark guess from Gaussian

%errorLimit = 20;% threshold for observation

% Step size between filter updates, can be less than 1.
Param.deltaT=0.1; % [s]
useRobust = false;
useGroundTruth = false;
numSteps = 200;
minK=Param.maxObs*30; % minimum number of range measurements to process initially
incK=Param.maxObs*5; % minimum number of range measurements to process after


%% Find and load data file
global Data;
Data = generateScript(Param.initialStateMean, numSteps, Param.maxObs, Param.alphas, Param.beta, Param.deltaT);
batchInitialization=true;


%% Set Noise parameters
noiseModels.prior = noiseModel.Diagonal.Sigmas([1 1 pi]');
noiseModels.pointPrior = noiseModel.Diagonal.Sigmas([100 100]');

if useRobust
  base = noiseModel.mEstimator.Tukey(15);
  noiseModels.range = noiseModel.Robust(base,noiseModel.Isotropic.Sigma(1, 10));
else
  noiseModels.range = noiseModel.Isotropic.Sigma(1, 10);
  %noiseModels.range = noiseModel.Diagonal.Sigmas([10 deg2rad(10)]');
  noiseModels.observ = noiseModel.Diagonal.Sigmas([deg2rad(10), 10]');
end

%% Initialize iSAM
isamModel = ISAM2;

%% Add prior on first pose
u = getControl(1);   % u = [rot1, trans, rot2], I assign the rot1 control to the last odometry control
pose0 = Pose2(Param.initialStateMean(1), Param.initialStateMean(2),Param.initialStateMean(3) + u(1));
newFactors = NonlinearFactorGraph;

newFactors.add(PriorFactorPose2(0,pose0,noiseModels.prior));


initial = Values;
initial.insert(0,pose0);
odo = Values;
odo.insert(0,pose0);
landmarkEstimates = Values;

field = getfieldinfo(Param.nLandmarksPerSide);
landmarkEverFaced = Values;

%% Loop over odometry
tic
update = false;
lastPose = pose0;
odoPose = pose0;
countK = 0;
k = 0; % observation counter
for i=1:length(Data.noisefreeControl) % M

  z = getObservations(i);
  u = getControl(i);
  rot1 = u(1); trans = u(2); rot2= u(3);
  odometry = Pose2(trans*cos(rot1),trans*sin(rot1),rot1+rot2);
  H = [-sin(rot1)*trans, cos(rot1), 0;
        cos(rot1)*trans,  sin(rot1), 0;
        1,0,1];
    M1 = Param.alphas(1) * rot1^2 + Param.alphas(2) * trans^2;
    M2 = Param.alphas(3) * trans^2 + Param.alphas(4) * rot2^2 + Param.alphas(4) * rot1^2;
    M3 = Param.alphas(1) * rot2^2 +Param.alphas(2) * trans^2;
    
    M = diag([M1,M2,M3]);
    odometry_noise = H * M * transpose(H);
    
  
  % add odometry factor
  newFactors.add(BetweenFactorPose2(i-1, i,...
      odometry,noiseModel.Gaussian.Covariance(odometry_noise)));
  
  % predict pose and update odometry
  predictedOdo = odoPose.compose(odometry);
  odoPose = predictedOdo;
  odo.insert(i,predictedOdo);
  
  % predict pose and add as initial estimate
  predictedPose = lastPose.compose(odometry);
  lastPose = predictedPose;
  initial.insert(i,predictedPose);
  landmarkEstimates.insert(i,predictedPose);
  
  % Check if there are range factors to be added
  for observe_index = 1 : size(z,2)
    bearing = z(2,observe_index);
    range = z(1,observe_index);
    id = z(3,observe_index);
    
    %% initialize points
    if ~landmarkEverFaced.exists(symbol('L',id))
        if useGroundTruth
          Lj = Point2(field.MARKER_X_POS(id),field.MARKER_Y_POS(id));
          %newFactors.add(PriorFactorPoint2(symbol('L',id),Lj,noiseModels.pointPrior));
        else
           shift_rotation = Pose2(0,0,minimizedAngle(bearing));
           moving_by_range = Pose2(range,0,0);
           pose = predictedPose.compose(shift_rotation).compose(moving_by_range);
           Lj = Point2(pose.x,pose.y);
           pose
           id
           [field.MARKER_X_POS(id),field.MARKER_Y_POS(id)]
         %newFactors.add(PriorFactorPoint2(symbol('L',id),Lj,noiseModels.pointPrior));
        end
        initial.insert(symbol('L',id),Lj);
        landmarkEstimates.insert(symbol('L',id),Lj);
        landmarkEverFaced.insert(symbol('L',id),Lj);
    end
    
    %factor = RangeFactorPosePoint2(i, symbol('L',id), range, noiseModels.range);
    factor = BearingRangeFactor2D(i, symbol('L',id), Rot2(bearing), range, noiseModels.observ);
    %factor   ange, noiseModels.observ);
    % Throw out obvious outliers based on current landmark estimates
    error=factor.unwhitenedError(landmarkEstimates);
%     if k<=minK || abs(norm(error))<errorLimit
%         newFactors.add(factor);
%     end
    newFactors.add(factor);
    k=k+1; countK=countK+1; update = true;
  end

  % Check whether to update iSAM 2
  if update && k>minK && countK>incK
    if batchInitialization % Do a full optimize for first minK ranges
      tic
      batchOptimizer = LevenbergMarquardtOptimizer(newFactors, initial);
      initial = batchOptimizer.optimize();
      toc
      batchInitialization = false; % only once
    end
    isamModel.update(newFactors, initial);
    result = isamModel.calculateEstimate();
    lastPose = result.at(i);
    % update landmark estimates
    landmarkEstimates = Values;
    for jj=1:Param.nLandmarksPerSide * 2
        key = symbol('L',jj);
        if landmarkEverFaced.exists(key)
            landmarkEstimates.insert(key,result.at(key));
        end
    end
    newFactors = NonlinearFactorGraph;
    initial = Values;
    countK = 0;
  end
      result = isamModel.calculateEstimate();
    XYT = utilities.extractPose2(result);
  plotsim(i,isamModel);
  
end
toc

%% Function definitions
%==========================================================================
function u = getControl(t)
global Data;
% noisefree control command
u = Data.noisefreeControl(:,t);  % 3x1 [drot1; dtrans; drot2]
end


%==========================================================================
function z = getObservations(t)
global Data;
% noisy observations
z = Data.realObservation(:,:,t); % 3xn [range; bearing; landmark id]
ii = find(~isnan(z(1,:)));
z = z(:,ii);
end
%==========================================================================
%%
function plotsim(t, isamModel)
global Data;
import gtsam.*

%--------------------------------------------------------------
% Graphics
%--------------------------------------------------------------

NOISEFREE_PATH_COL = 'green';
ACTUAL_PATH_COL = 'blue';

NOISEFREE_BEARING_COLOR = 'cyan';
OBSERVED_BEARING_COLOR = 'red';

GLOBAL_FIGURE = 1;

%=================================================
% data *not* available to your filter, i.e., known
% only by the simulator, useful for making error plots
%=================================================
% actual position (i.e., ground truth)
x = Data.Sim.realRobot(1,t);
y = Data.Sim.realRobot(2,t);
theta = Data.Sim.realRobot(3,t);

% real observation
observation = Data.realObservation(:,:,t);

% noisefree observation
noisefreeObservation = Data.Sim.noisefreeObservation(:,:,t);

%=================================================
% graphics
%=================================================
figure(GLOBAL_FIGURE); clf; hold on; plotfield(observation(3,:));

% draw actual path (i.e., ground truth)
plot(Data.Sim.realRobot(1,1:t), Data.Sim.realRobot(2,1:t), 'Color', ACTUAL_PATH_COL);
plotrobot( x, y, theta, 'black', 1, ACTUAL_PATH_COL);


% draw noise free motion command path
plot(Data.Sim.noisefreeRobot(1,1:t), Data.Sim.noisefreeRobot(2,1:t), 'Color', NOISEFREE_PATH_COL);
plot(Data.Sim.noisefreeRobot(1,t), Data.Sim.noisefreeRobot(2,t), '*', 'Color', NOISEFREE_PATH_COL);

for k=1:size(observation,2)
    rng = Data.Sim.noisefreeObservation(1,k,t);
    ang = Data.Sim.noisefreeObservation(2,k,t);
    noisy_rng = observation(1,k);
    noisy_ang = observation(2,k);

    % indicate observed range and angle relative to actual position
    plot([x x+cos(theta+noisy_ang)*noisy_rng], [y y+sin(theta+noisy_ang)*noisy_rng], 'Color', OBSERVED_BEARING_COLOR);

   
    % result
    result = isamModel.calculateEstimate();
    XYT = utilities.extractPose2(result);
    plot(XYT(:,1),XYT(:,2),'r-');
    XY = utilities.extractPoint2(result);
    plot(XY(:,1),XY(:,2),'r*');
     
    % indicate ideal noise-free range and angle relative to actual position
    plot([x x+cos(theta+ang)*rng], [y y+sin(theta+ang)*rng], 'Color', NOISEFREE_BEARING_COLOR);
end
end
