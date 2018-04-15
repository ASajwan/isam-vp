function [BatchIndex,BatchExist] = NNBatch(predLandMarkBatch,graph, LandMarkCount, threshold)
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

% calculate the nearest neighboor distance cost
% Dij = || x_pred_i - x_landmark_j ||^2

%Suppose there are k measurement at each time while there are n landmarks
%For ith row, there are the distances between current predicted landmark,
%based on the rob pose and ith measurement and all the landmarks
%following with D_max in indicate new landmarks

% D = [D_11, D_12, D_13....   D_1n, D_max... D_max
%      D_21, D_22, D_23....   D_2n, D_max... D_max
%      ...     ...   ...       ...   ...       ...
%      D_k1, D_k2, D_k3....   D_kn, D_max... D_max]


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
try
    BatchIndex = lapjv(D);
    BatchExist = BatchIndex <= LandMarkCount;
catch
    [BatchIndex, BatchExist]= NNSingle(graph, LandMarkCount,predLandMarkBatch, threshold);
end
end


