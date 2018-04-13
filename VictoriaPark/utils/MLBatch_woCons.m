function  [index, exist] = MLBatch_woCons(poseID, pose, z, graph, nodeEstimate, LandMarkCount, threshold)
%
% Description: This function provides the index of the landmark
% correspoding to a measurement based on maximum likehood using Mahalanobis
% distance, given the landmarks' coordinates, robot's pose and measurement
%
% Arguments:
% Inputs-
%       poseID - ID of the pose from where the measurement was taken
%       pose   - Pose of the robot when the measurement was taken
%       z      - Measurement [bearing; range]
%       graph  - Factor graph till now
%       nodeEstimate - Numerical estimates of the nodes
%       LandMarkCount - No. of landmarks seen
%       threshold - Threshold sq. distance for declaring a new landmark
%
% Outputs-
%       index - index of the associated landmark
%       exist - boolean array indicating if landmark exists for the
%       corresponding measurement
%
marginal = gtsam.Marginals(graph, nodeEstimate);
KeyVector = gtsam.KeyVector;
KeyVector.push_back(poseID);
for k = 1: LandMarkCount
    key = gtsam.symbol('L',double(k));
    KeyVector.push_back(key);
end
fullMatrix = marginal.jointMarginalCovariance(KeyVector).fullMatrix;

global Param;
D = ones(size(z,2), LandMarkCount + ...
        max(1,floor(size(z,2)/2)));
D = D * threshold;
for k = 1: size(z,2)
    for i = 1:LandMarkCount
        key = gtsam.symbol('L',double(i));
        row = [1:3,3+2*i-1:3+2*i];
        Cov = fullMatrix(row,row);

        delta = [nodeEstimate.at(key).x, nodeEstimate.at(key).y] - [pose.x,pose.y] ;
        delta_x = delta(1);   delta_y = delta(2);
        q = delta * delta';    
        H = 1/q * [-sqrt(q)*delta_x, -sqrt(q)*delta_y, 0, sqrt(q)*delta_x, sqrt(q)*delta_y
                   delta_y, -delta_x, -q, -delta_y, delta_x];
        diff = z(:,k) - [sqrt(q); minimizedAngle(atan2(delta_y, delta_x) - pose.theta + pi/2)];
        %diff(2) = minimizedAngle(diff(2));
        diff(2) = minimizedAngle(diff(2));
%         diff = rot90(rot90(diff));

        Q = H * Cov * H' + Param.R;
        distance = diff' * (Q\diff);
        D(k,i) = distance;
    end
end
index = ones(size(z,2),1);
for i = 1: size(z,2)
    smallestIndex = LandMarkCount + 1;
    distance = threshold;
    for j = 1:  LandMarkCount
        if D(i,j) < distance
            distance = D(i,j);
            smallestIndex = j;
        end
    end
    index(i) = smallestIndex;
end
exist = index <= LandMarkCount;
end

