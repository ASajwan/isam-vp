function [BatchIndex,BatchExist] = NNBatch_woCons(predLandMarkBatch,graph, LandMarkCount, threshold)
%
% Description: This function provides the indices of landmarks
% correspoding to measurements based on nearest neighbour algorithm, given
% the landmarks' coordinates, robot's pose and measurement
%
% Arguments:
% Inputs-
%       graph  - Factor graph till now
%       LandMarkCount - No. of landmarks seen
%       predLandMarkBatch - Predicted Landmark coordinates
%       threshold - Threshold sq. distance for declaring a new landmark
%
% Outputs-
%       BatchIndex - indices of the associated landmarks
%       BatchExist - boolean array indicating if landmarks exist for the
%       corresponding measurements
%
import gtsam.*                              
AllLandMarks = zeros(LandMarkCount,2);
for j = 1:LandMarkCount
  key = symbol('L',j);
  AllLandMarks(j,:) = [graph.at(key).x, graph.at(key).y];
end
D = ones(size(predLandMarkBatch,1), LandMarkCount + ...
        max(1,floor(size(predLandMarkBatch,1)/3)));
D = D * threshold;
X_diff = predLandMarkBatch(:,1) - AllLandMarks(:,1)';
Y_diff = predLandMarkBatch(:,2) - AllLandMarks(:,2)';
D(1:size(predLandMarkBatch,1), 1: LandMarkCount) = X_diff.*X_diff + Y_diff.*Y_diff;


BatchIndex = ones(size(predLandMarkBatch,1),1);
for i = 1: size(predLandMarkBatch,1)
    smallestIndex = LandMarkCount + 1;
    distance = threshold;
    for j = 1:  LandMarkCount
        if D(i,j) < distance
            distance = D(i,j);
            smallestIndex = j;
        end
    end
    BatchIndex(i) = smallestIndex;
end
BatchExist = BatchIndex <= LandMarkCount;
end


