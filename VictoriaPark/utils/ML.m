function [nearest, exist]= ML(poseID, pose, z, graph, nodeEstimate, LandMarkCount,predLandMark, threshold,fullMatrix)
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
%       predLandMark - Predicted landmarks' coordinates
%       threshold - Threshold sq. distance for declaring a new landmark
%       fullMatrix - Full covariance matrix
%
% Outputs-
%       index - index of the associated landmark
%       exist - boolean array indicating if landmark exists for the
%       corresponding measurement
%
global Param;


nearest = gtsam.symbol('L',0);
exist = false;
nearestDist = threshold;


    
for i = 1:LandMarkCount
    key = gtsam.symbol('L',double(i));
    if nodeEstimate.exists(key)  %insanity check
        row = [1:3,3+2*i-1:3+2*i];
        Cov = fullMatrix(row,row);
        %% strange, it currently is 3 x 2, it is not symmetric
% % % %         %Attention !!!! Cov will be like [bearing, bearing_range
% % % %                                           %bearing_range, range];
                                          

%         Cov_shift = Cov;
%         Cov_shift(1,:) = Cov(2,:);
%         Cov_shift(2,:) = Cov(1,:);
%         Cov_shift(:,1) = Cov(:,2);
%         Cov_shift(:,2) = Cov(:,1);


        delta = [nodeEstimate.at(key).x, nodeEstimate.at(key).y] - [pose.x,pose.y] ;
        delta_x = delta(1);   delta_y = delta(2);
        q = delta * delta';    
        H = 1/q * [delta_y, -delta_x, -q, -delta_y, delta_x;
                   -sqrt(q)*delta_x, -sqrt(q)*delta_y, 0, sqrt(q)*delta_x, sqrt(q)*delta_y];

        diff = z - [sqrt(q); atan2(delta_y, delta_x) - pose.theta + pi/2];
        diff(2) = minimizedAngle(diff(2));
        diff = rot90(rot90(diff));

        Q = H * Cov * H' + rot90(rot90(Param.R));
        distance = diff' * (Q\diff);

        if distance < nearestDist
            nearest = key;
            nearestDist = distance;
            exist = true;
        end
    end
end


end
