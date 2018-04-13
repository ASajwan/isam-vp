function [Li,indexed_z] = da_nn(z, Mu, Sigma)
% perform nearest-neighbor data association
global Param;
global State;

[~, n_observation] = size(z);
Li = ones(1,n_observation) * -1;
index_start = State.Ekf.nL + 1;

for i = 1: n_observation
    Distance_History = [];
    n_landmark = (length(Mu) - 3 ) /2;
    threshHold = 15;   %15 works for the task simulator    %200 for task vp
    min_distance = 1000;
    min_index = -1;
    for j = 1: n_landmark
        landmark_x = Mu(3+2*(j)-1);
        landmark_y = Mu(3+2*(j));
    

        % Compute expected observation and Jacobian
        rob_x = Mu(1); rob_y = Mu(2);  
        rob_theta = minimizedAngle(Mu(3));

        F = zeros(5, 3+2*(State.Ekf.nL));
        F(1:3,1:3) = eye(3);
        F(4:5,3+2*j - 1: 3+2*j) = eye(2);

         delta_x = landmark_x - rob_x;
         delta_y = landmark_y - rob_y;
         q = delta_x^2 +  delta_y^2;

         H = [-sqrt(q)*delta_x, -sqrt(q)*delta_y, 0,sqrt(q)*delta_x,sqrt(q)*delta_y;
                delta_y, -delta_x, -q,  -delta_y, delta_x]/q * F;
            
            
         S = H * Sigma * transpose(H) + Param.R;
            
         theta_diff = minimizedAngle(z(2,i)-(atan2(delta_y, delta_x) - rob_theta));
         distance_diff = z(1,i) - sqrt(q);
         diff = [distance_diff; theta_diff];
         
         distance = diff' / S * diff;
         Distance_History = [Distance_History,distance];
         
         if distance < threshHold && distance < min_distance
             min_distance = distance;
             min_index = j;
         end
    end
    if min_index ~= -1
        Li(i) = min_index;
    else 
        %hope the efkpredict_sim will find this is new landmark and
        %augument the state
        Li(i) = index_start;
        index_start = index_start + 1;
    end
end
indexed_z = [z;Li];


% Innovation / residual covariance
end


function [x,y] = predictLandMark(z,predMu)
    xRob = predMu(1); yRob = predMu(2); thetaRob = predMu(3);
    angle = minimizedAngle(minimizedAngle(z(2)) + minimizedAngle(thetaRob));
    x = xRob + z(1) * cos(angle);
    y = yRob + z(1) * sin(angle);
    
end
