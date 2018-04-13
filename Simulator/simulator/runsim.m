function varargout = runsim(stepsOrData, pauseLen)

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
Param.maxObs = 4;

% number of landmarks per sideline of field (minimum is 3)
Param.nLandmarksPerSide = 6;

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

traj_cov = [];
makeVideo = 0;

if makeVideo
    try
        votype = 'avifile';
        vo = avifile('video.avi', 'fps', min(5, 1/pauseLen));
    catch
        votype = 'VideoWriter';
        vo = VideoWriter('video', 'Motion JPEG AVI');
        v.Quality = 100;
        %set(vo, 'FrameRate', min(5, 1/pauseLen));
        open(vo);
    end
end

dA_correct = [];
dA_actual = [];

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
    
    [Mu, Sigma, dataAssociation]=ekfupdate(z,Mu, Sigma,u);
    dA_actual = [dA_actual;dataAssociation'];
    dA_correct = [dA_correct;z'];
    
    %=================================================
    %TODO: plot and evaluate filter results here
    %=================================================
     
     mean_x = Mu(1); mean_y = Mu(2);
     mean = [mean_x;mean_y];
     plot_Sigma = Sigma(1:2,1:2);
     draw_ellipse(mean,plot_Sigma,3);
     traj(t+1,:) = mean';
     
     [~, cov_size] = size(Sigma);
     
     if cov_size == 3 + 2 *Param.nLandmarksPerSide*2
         traj_cov = [traj_cov; Sigma];
     end


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
plot_DA(dA_actual,dA_correct);
plotHW(Sigma, traj_cov);

if nargout >= 1
    varargout{1} = Data;
end


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

function plotHW(Sigma, traj_cov)
global Param
    [~, cov_size] = size(Sigma);
    sigma_each = zeros(cov_size);
    for i = 1 : cov_size
        sigma_each(i) = sqrt(Sigma(i,i));
    end
    for i = 1: cov_size
        for j = 1:cov_size
            Sigma(i,j) = Sigma(i,j) / (sigma_each(i) * sigma_each(j));
        end
    end
    Sigma = 1 - Sigma;
    I = mat2gray(Sigma);
    I = imresize(I, 20, 'nearest');
    figure;
    imshow(I);
    
    
    cov_landmarks = zeros(Param.nLandmarksPerSide*2,Param.nLandmarksPerSide*2);
    Sigma = 1 - Sigma;
    for i = 1 : Param.nLandmarksPerSide*2
        for j = 1: Param.nLandmarksPerSide*2
            cov_landmarks(i,j) = Sigma(3+2*i-1, 3+2*j-1);
        end
    end
    figure;
    plot(cov_landmarks)
    legend('landmark1','landmark2', 'landmark3','landmark4',...
           'landmark5','landmark6','landmark7','landmark8');
    xlabel('landmark id')
    ylabel('\rho')
    title('correlation cofficients \rho between landmarks')
    
    [row, ~] = size(traj_cov);
    n = 3 + 2 * 2 * Param.nLandmarksPerSide;
    time = row / n;
    
    det_plot = zeros(time, Param.nLandmarksPerSide*2);
    for t = 1: time
        sigma_that_time = traj_cov((t-1)*n + 1:t*n,:);
        for i = 1:Param.nLandmarksPerSide*2
            submatrix = sigma_that_time(3+2*i-1:3+2*i, 3+2*i-1:3+2*i);
            det_plot(t,i) = log(det(submatrix));
        end
    end
    figure;
    plot(det_plot)
%     legend('landmark1','landmark2', 'landmark3','landmark4',...
%            'landmark5','landmark6','landmark7','landmark8');
    xlabel('time t')
    ylabel('log det_{\Sigma_{ii}}')
    title('determinant of each landmarks covariance over time')
    
           
function plot_DA(dA_actual, dA_correct)
    global Param
    landmark = zeros(Param.nLandmarksPerSide*2,1);
    time = length(dA_correct);
    correct_da_hist = zeros(time,1);
    assign_da_hist = zeros(time,1);
    for t = 1:time
        correct_da = dA_correct(t,3);
        assign_da = dA_actual(t,3);
        if landmark(correct_da,1) == 0
            landmark(correct_da,1) = assign_da;
            assign_da_hist(t,1) = correct_da;
        else
            if assign_da == landmark(correct_da,1)
                assign_da_hist(t,1) = correct_da;
            else
                assign_da_hist(t,1) = -1;
            end
        end
        correct_da_hist(t,1) = correct_da;
    end
    figure;
    scatter(1:t, assign_da_hist);
    hold on
    scatter(1:t, correct_da_hist);
    legend('actual association','correct association');
    xlabel('time t')
    ylabel('landmark')
    title('Data Association over Time')
    
    





