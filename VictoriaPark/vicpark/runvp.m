function runvp(nSteps,pauseLen)

global Param;
global State;
global Data;

if ~exist('nSteps','var') || isempty(nSteps)
    nSteps = inf;
end

if ~exist('pauseLen','var')
    pauseLen = 0; % seconds
end

if Param.makeVideo
    try
        votype = 'avifile';
        vo = avifile('video.avi', 'fps', min(5, 1/pauseLen));
    catch
        votype = 'VideoWriter';
        vo = VideoWriter('video', 'MPEG-4');
        set(vo, 'FrameRate', min(5, 1/pauseLen));
        open(vo);
    end
end

Data = load_vp_si();

% Initalize Params
%===================================================
% vehicle geometry
Param.a = 3.78; % [m]
Param.b = 0.50; % [m]
Param.L = 2.83; % [m]
Param.H = 0.76; % [m]

% 2x2 process noise on control input
sigma.vc = 0.02; % [m/s]
sigma.alpha = 2*pi/180; % [rad]
Param.Qu = diag([sigma.vc, sigma.alpha].^2);

% 3x3 process noise on model error
sigma.x = 0.1; % [m]
sigma.y = 0.1; % [m]
sigma.phi = 0.5*pi/180; % [rad]
Param.Qf = diag([sigma.x, sigma.y, sigma.phi].^2);

% 2x2 observation noise
sigma.r = 0.05; % [m]
sigma.beta = 1*pi/180; % [rad]
Param.Q = diag([sigma.r, sigma.beta].^2);

Param.R = .99;
%===================================================

% Initialize State
%===================================================
State.Ekf.mu = [Data.Gps.x(2), Data.Gps.y(2), 36*pi/180]';
State.Ekf.Sigma = zeros(3);
predMu = State.Ekf.mu;
predSigma = State.Ekf.Sigma;

global AAr;
AAr = [0:360]*pi/360;


figure(1); clf;
axis equal;

ci = 1; % control index
t = min(Data.Laser.time(1), Data.Control.time(1));
for k=1:min(nSteps, length(Data.Laser.time))
    
    while (Data.Control.time(ci) < Data.Laser.time(k))
       % control available
       dt = Data.Control.time(ci) - t;
       t = Data.Control.time(ci);
       u = [Data.Control.ve(ci), Data.Control.alpha(ci)]';

       tstart = tic;
       [predMu, predSigma] = ekfpredict_vp(u, dt);
       Data.stopwatch{1}(ci) = Data.Control.time(ci);
       Data.stopwatch{2}(ci) = toc(tstart);
       
       ci = ci+1;
    end
    
    % observation available
    dt = Data.Laser.time(k) - t;
    t = Data.Laser.time(k);
    z = detectTreesI16(Data.Laser.ranges(k,:));
    
    tstart = tic;
    ekfupdate(z,predMu,predSigma)
    Data.stopwatch{3}(k) = Data.Laser.time(k);
    Data.stopwatch{4}(k) = toc(tstart);
    Data.stopwatch{5}(k) = length(State.Ekf.mu);
    Data.trajectory.x(k) = State.Ekf.mu(1);
    Data.trajectory.y(k) = State.Ekf.mu(2);
    
    doGraphics(z);
    drawnow;
    if pauseLen > 0
        pause(pauseLen);
    end
    
    if Param.makeVideo
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

if Param.makeVideo
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
function doGraphics(z)
% Put whatever graphics you want here for visualization
%
% WARNING: this slows down your process time, so use sparingly when trying
% to crunch the whole data set!

global Param;
global State;

% plot the robot and 3-sigma covariance ellipsoid
plotbot(State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.mu(3), 'black', 1, 'blue', 1);
hold on;

plotcov2d( State.Ekf.mu(1), State.Ekf.mu(2), State.Ekf.Sigma, 'blue', 0, 'blue', 0, 3);

for i = 4:2:length(State.Ekf.mu)
    plotcov2d(State.Ekf.mu(i),State.Ekf.mu(i+1),State.Ekf.Sigma(i:i+1,i:i+1),'r',0,'r',0,9);
%                 plotcov2d(State.Ekf.mu(i),State.Ekf.mu(i+1),[.1 0;0 .1],'r',0,'r',0,3);
end

% restrict view to a bounding box around the current pose
BB=20;
axis([[-BB,BB]+State.Ekf.mu(1), [-BB,BB]+State.Ekf.mu(2)]);
% axis([-100 0 -100 0])

% project raw sensor detections in global frame using estimate pose
xr = State.Ekf.mu(1);
yr = State.Ekf.mu(2);
tr = State.Ekf.mu(3);
for k=1:size(z,2)
    r = z(1,k);
    b = z(2,k);
    xl = xr + r*cos( minimizedAngle(b+tr-pi/2) );
    yl = yr + r*sin( minimizedAngle(b+tr-pi/2) );
    plot([xr; xl], [yr; yl],'r',xl,yl,'r*');
end


hold off;

