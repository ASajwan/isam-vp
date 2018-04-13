function [augumentedMu,augumentedSigma] = ekfpredict_sim(z,mu, sigma)

% EKF-SLAM prediction for simulator process model

global Param;
global State;



% Augument the state if needed
augumentedSigma = sigma;
augumentedMu = mu;
[~, n_observation] = size(z);
for i = 1:n_observation
    if contains(z(3,i)) == 1  %continue, omit
        continue;
    else
        State.Ekf.sL = [State.Ekf.sL, z(3,i)];
        State.Ekf.iL    = [State.Ekf.iL, State.Ekf.nL + 1];
        State.Ekf.nL = State.Ekf.nL + 1;
        [landmarkX, landmarkY] = predictLandMark([z(1,i), z(2,i)],mu);

        augumentedMu = zeros(3+2*(State.Ekf.nL),1);
        augumentedMu(4:3+2*(State.Ekf.nL-1)) = mu(4:end);
        augumentedMu(1:3) = mu(1:3);
        augumentedMu(3+2*(State.Ekf.nL)-1:end) = [landmarkX,landmarkY];
        augumentedMu(3) = minimizedAngle(augumentedMu(3));
        
        %augument the state covariance, Gr is the jacobian matrix with repesct to
        %rob position, (2*3)
        angle_sum = minimizedAngle(minimizedAngle(mu(3)) +minimizedAngle(z(2,i)));
        r = z(1,i);
        Gr= [1,0,-r*sin(angle_sum);
            0,1,r*cos(angle_sum)];
        %Gdelta is the jacobian matrix with respect to error (2*2)
        Gdelta = [cos(angle_sum), -r*sin(angle_sum);
                   sin(angle_sum), r*cos(angle_sum)];
        
        sigmaNewRob = Gr*sigma(1:3,1:3);   %2*3
        sigmaNewLand = Gr*sigma(1:3, 4:end); %2*3, 3*2n = 2*2n
        sigmaNewNew = Gr*sigma(1:3,1:3)*transpose(Gr) + ...
                                Gdelta*Param.R*transpose(Gdelta);
                            
        newLength = 3+2*State.Ekf.nL;
        augumentedSigma = zeros(newLength, newLength);
        augumentedSigma(1:newLength-2, 1:newLength-2) = sigma;
        augumentedSigma(newLength -1:end, 1:3) = sigmaNewRob;
        augumentedSigma(1:3, newLength -1:end)= transpose(sigmaNewRob);
        augumentedSigma(newLength -1:end, 4:newLength-2) = sigmaNewLand;
        augumentedSigma(4:newLength-2, newLength -1:end) = transpose(sigmaNewLand);
        augumentedSigma(newLength-1:end, newLength -1:end) = sigmaNewNew;
        
    end
    mu = augumentedMu;
    sigma = augumentedSigma;
end
end

function exist = contains(landmark)
    global State;
    exist = -1;
    for i = 1: length(State.Ekf.sL)
        if State.Ekf.sL(i) == landmark
            exist = 1;
            break;
        end
    end
end



% z = [distance, angle]
function [x,y] = predictLandMark(z,predMU)
    xRob = predMU(1); yRob = predMU(2); thetaRob = predMU(3);
    angle = minimizedAngle(minimizedAngle(z(2)) + minimizedAngle(thetaRob));
    x = xRob + z(1) * cos(angle);
    y = yRob + z(1) * sin(angle);
    
end


