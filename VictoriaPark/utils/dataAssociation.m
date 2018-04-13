function [nearest,exist] = dataAssociation(type,poseID,pose,z, graph, factor, LandMarkCount,predLandMark,threshold,fullMatrix)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
switch type
    case 'ML'
        [nearest, exist]= ML(poseID,pose,z, graph, factor, LandMarkCount,predLandMark,threshold,fullMatrix);
    case 'NN'
        [nearest, exist]= NN(poseID,pose,z, graph, factor, LandMarkCount,predLandMark,threshold,fullMatrix);
end
        
end

