function varargout = runsim_2(stepsOrData, pauseLen)

global Param;
global Data;
global State;

if ~exist('pauseLen','var')
    pauseLen = 0.3; % seconds
end

% Initalize Params
%===================================================
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

% Step size between filter updates, can be less than 1.
Param.deltaT=0.1; % [s]

if isscalar(stepsOrData)
    % Generate a data set of motion and sensor info consistent with
    % noise models.
    numSteps = stepsOrData;
    Data = generateScript(Param.initialStateMean, numSteps, Param.maxObs, Param.alphas, Param.beta, Param.deltaT);
else
    % use a user supplied data set from a previous run
    Data = stepsOrData;
    numSteps = size(Data, 1);
end
%===================================================

% Initialize State
%===================================================
State.Ekf.mu = Param.initialStateMean;
State.Ekf.Sigma = eye(3);

Mu = Param.initialStateMean;
Sigma = eye(3);
traj = zeros(numSteps+1,2);
traj(1,:) = [180,50];

for t = 1:numSteps
    plotsim(t, traj);

    %=================================================
    % data available to your filter at this time step
    %=================================================
    u = getControl(t);
    z = getObservations(t);


    %=================================================
    %TODO: update your filter here based upon the
    %      motionCommand and observation
    %=================================================
    
    [Mu, Sigma]=ekfupdate(z,Mu, Sigma,u);
    %=================================================
    %TODO: plot and evaluate filter results here
    %=================================================
     
     mean_x = Mu(1); mean_y = Mu(2);
     mean = [mean_x;mean_y];
     plot_Sigma = Sigma(1:2,1:2);
     draw_ellipse(mean,plot_Sigma,3);
     traj(t+1,:) = mean';

     
     
     
     
     if length(Mu) > 3
         for i = 1: (length(Mu) - 3) /2
             mean = [Mu(3+2*(i)-1); Mu(3+2*i)];
             plot_Sigma = Sigma(3+2*(i)-1:3+2*i,3+2*(i)-1:3+2*i);
             draw_ellipse(mean,plot_Sigma,3);
         end
     end

    drawnow;
    if pauseLen > 0
        pause(pauseLen);
    end
end

if nargout >= 1
    varargout{1} = Data;
end

%==========================================================================
function u = getControl(t)
global Data;
% noisefree control command
u = Data.noisefreeControl(:,t);  % 3x1 [drot1; dtrans; drot2]


%==========================================================================
function z = getObservations(t)
global Data;
% noisy observations
z = Data.realObservation(:,:,t); % 3xn [range; bearing; landmark id]
ii = find(~isnan(z(1,:)));
z = z(:,ii);

%==========================================================================
function plotsim(t, traj)
global Data;

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

     %draw predict path
     plot(traj(1:t,1), traj(1:t,2), 'Color', [0.7 0.7 0.7]);
     
    % indicate ideal noise-free range and angle relative to actual position
    plot([x x+cos(theta+ang)*rng], [y y+sin(theta+ang)*rng], 'Color', NOISEFREE_BEARING_COLOR);
end
