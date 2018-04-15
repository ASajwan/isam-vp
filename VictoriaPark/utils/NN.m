function [index, exist]= NN(factor, LandMarkCount,predLandMarkBatch, threshold)
%
% Description: This function provides the index of the landmark
% correspoding to a measurement based on nearest neighbour algorithm, given
% the landmarks' coordinates, robot's pose and measurement
%
% Arguments:
% Inputs-
%       poseID - ID of the pose from where the measurement was taken
%       pose   - Pose of the robot when the measurement was taken
%       z      - Measurement [bearing; range]
%       graph  - Factor graph till now
%       LandMarkCount - No. of landmarks seen
%       predLandMarkBatch - Predicted Landmark coordinates
%       threshold - Threshold sq. distance for declaring a new landmark
%
% Outputs-
%       index - index of the associated landmark
%       exist - boolean array indicating if landmark exists for the
%       corresponding measurement
%

index = zeros(length(predLandMarkBatch),1);
exist = zeros(length(predLandMarkBatch),1);

for j = 1 : length(predLandMarkBatch)

    nearest = -1;
    nearestDist = threshold;
    ext = false;
    predLandMark = predLandMarkBatch(j);

    for i = 1:LandMarkCount
        key = gtsam.symbol('L',double(i));
        Lx = factor.at(key).x;
        Ly = factor.at(key).y;
        
        % Find distance between ldmk predicted using measurement and pose,
        % and the actual ldmk
        dist =  (Lx-predLandMark.x)^2 + (Ly-predLandMark.y)^2 ;
        
        % Choose the ldmk with the smallest distance yet
        if dist < nearestDist
            nearest = i;
            ext = true;
            nearestDist = dist;
        end
    end
    index(j) = nearest;
    exist(j) = ext;

end

end

