function [UpdatedMu, UpdatedSigma, z]=ekfupdate(z,mu, sigma,u)
% EKF-SLAM update step for both simulator and Victoria Park data set

global Param;
global State;

% --------------------------------------------
% Prediction step
% --------------------------------------------

% EKF prediction of mean and covariance
x = mu(1); y = mu(2); theta = minimizedAngle(mu(3));
rot1 = minimizedAngle(u(1)); tran = u(2); rot2 = minimizedAngle(u(3));

M1 = Param.alphas(1)*rot1^2 + Param.alphas(2) * tran^2;
M2 = Param.alphas(3) * tran^2 + Param.alphas(4) * rot2^2 + Param.alphas(4) * rot1^2;
M3 = Param.alphas(1) * rot2^2 +Param.alphas(2) * tran^2;
M = [M1,0,0;
    0, M2,0;
    0,0, M3];

angle_sum = minimizedAngle(theta+rot1);


% Gdela is the jacobian with respect to error
% F = zeros(3, 3+2*(State.Ekf.nL));
% F(1:3,1:3) = eye(3);

V= [-tran * sin(angle_sum), cos(angle_sum), 0;
      tran * cos(angle_sum), sin(angle_sum), 0;
      1,0,1];
  
G =  [1,0, -tran * sin(angle_sum);
      0,1, tran * cos(angle_sum);
      0,0,1];
          
R = V * M * transpose(V);
sigma(1:3,1:3) = G * sigma(1:3,1:3) *  transpose(G) + R;

sigma(1:3,4:end) = G *sigma(1:3,4:end);
sigma(4:end, 1:3) = sigma(1:3,4:end)';

newRob = predictionNew(mu(1:3), u);
mu(1:3) = newRob;
          
% G = eye(3 + 2*(State.Ekf.nL)) + transpose(F)*G_local *(F);
% R_local = V * M * transpose(V);
% R = transpose(F)*R_local *(F);
% newRob = predictionNew(mu(1:3), u);
% mu(1:3) = newRob;
%  
% sigma = G * sigma * transpose(G) + R;


% returns state vector indices pairing observations with landmarks
switch lower(Param.dataAssociation)
    case 'known'
        Li = da_known(z(3,:));
    case 'nn'
        [Li,z] = da_nn(z(1:2,:), mu, sigma);
    case 'jcbb'
        Li = da_jcbb(z(1:2,:), Param.R);
    otherwise
        error('unrecognized data association method: "%s"', Param.dataAssociation);
end


[temp, n_observation] = size(z);

%--------------------------------------------------------------
% Correction step
%--------------------------------------------------------------
UpdatedMu = mu;
UpdatedSigma = sigma;
debug  = 0;

landmarks = getfieldinfo(4);

if strcmp(Param.updateMethod,'batch')
    H_stack =[];
    diff_stack = [];
    count = 0;
    
    for i = 1:n_observation
        index = findIndex(z(3,i));
        if index == -1
            continue
        end
        count = count + 1;
        
%         landmark_x = landmarks.MARKER_X_POS(z(3,i));
%         landmark_y = landmarks.MARKER_Y_POS(z(3,i));
        
        landmark_x = UpdatedMu(3+2*(index)-1);
        landmark_y = UpdatedMu(3+2*(index));

        % Compute expected observation and Jacobian
        rob_x = UpdatedMu(1); rob_y = UpdatedMu(2);  
        rob_theta = minimizedAngle(UpdatedMu(3));

        F = zeros(5, 3+2*(State.Ekf.nL));
        F(1:3,1:3) = eye(3);
        F(4:5,3+2*index - 1: 3+2*index) = eye(2);
           

         delta_x = landmark_x - rob_x;
         delta_y = landmark_y - rob_y;
         q = delta_x^2 +  delta_y^2;

         H = [-sqrt(q)*delta_x, -sqrt(q)*delta_y, 0,sqrt(q)*delta_x,sqrt(q)*delta_y;
                delta_y, -delta_x, -q,  -delta_y, delta_x]/q * F;

%          H = [-sqrt(q)*delta_x, -sqrt(q)*delta_y, 0
%                 delta_y, -delta_x, -q]/q;
            
         theta_diff = minimizedAngle(z(2,i)-(atan2(delta_y, delta_x) - rob_theta));
         distance_diff = z(1,i) - sqrt(q);
         diff = [distance_diff; theta_diff];
         
         H_stack = [H_stack;H];
         diff_stack = [diff_stack; diff];
    end
    if count > 0
        % Correction
            
        R = zeros(2*count,2*count);
        for i = 1:count
            R(2*(i-1)+1:2*i, 2*(i-1)+1:2*i) = Param.R;
        end
        % Innovation / residual covariance
%          S = H * UpdatedSigma * transpose(H) + Param.R;

         S_stack = H_stack * UpdatedSigma * transpose(H_stack) + R;

        % Kalman gain
%          K = UpdatedSigma * transpose(H) / S;  % (3+2n)*2;
           K_stack = UpdatedSigma * transpose(H_stack) / S_stack;  % (3+2n)*2;
           UpdatedMu = UpdatedMu + K_stack * diff_stack;
           UpdatedSigma = (eye(3 + 2*(State.Ekf.nL)) - K_stack * H_stack) * UpdatedSigma;

    end
else
    for i = 1:n_observation
        index = findIndex(z(3,i));
        if index == -1
            continue
        end

        landmark_x = UpdatedMu(3+2*(index)-1);
        landmark_y = UpdatedMu(3+2*(index));
        % Compute expected observation and Jacobian
        rob_x = UpdatedMu(1); rob_y = UpdatedMu(2);  
        rob_theta = minimizedAngle(UpdatedMu(3));

        F = zeros(5, 3+2*(State.Ekf.nL));
        F(1:3,1:3) = eye(3);
        F(4:5,3+2*index - 1: 3+2*index) = eye(2);

         delta_x = landmark_x - rob_x;
         delta_y = landmark_y - rob_y;
         q = delta_x^2 +  delta_y^2;

         H = [-sqrt(q)*delta_x, -sqrt(q)*delta_y, 0,sqrt(q)*delta_x,sqrt(q)*delta_y;
                delta_y, -delta_x, -q,  -delta_y, delta_x]/q * F;




        % Innovation / residual covariance
         S = H * UpdatedSigma * transpose(H) + Param.R;

        % Kalman gain
         K = UpdatedSigma * transpose(H) / S;  % (3+2n)*2;

        % Correction
         theta_diff = minimizedAngle(z(2,i)-(atan2(delta_y, delta_x) - rob_theta));
         distance_diff = z(1,i) - sqrt(q);
         diff = [distance_diff; theta_diff];

         UpdatedMu = UpdatedMu + K * diff;
         UpdatedSigma = (eye(3+2*(State.Ekf.nL)) - K * H ) * UpdatedSigma;

    end
end


% --------------------------------------------
% Augument step
% --------------------------------------------
[UpdatedMu,UpdatedSigma] = ekfpredict_sim(z,UpdatedMu, UpdatedSigma);

 
 
end



function NewRob = predictionNew(state, motion)
NewRob(3)=minimizedAngle(state(3)+motion(1));
NewRob(1)=state(1)+motion(2)*cos(state(3));
NewRob(2)=state(2)+motion(2)*sin(state(3));
NewRob(3)=NewRob(3)+motion(3);
NewRob(3)=minimizedAngle(NewRob(3));
end

function index = findIndex(landID)
    global State;
    index = -1;
    for i = 1:State.Ekf.nL
        if State.Ekf.sL(i) == landID
            index = State.Ekf.iL{i};
            break;
        end
    end
 
end