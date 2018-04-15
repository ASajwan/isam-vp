# isam-VicPark
Decription:
This dataset comes from http://www-personal.acfr.usyd.edu.au/nebot/victoria_park.htm
The video are produced based on one-step update frequency and Constrained Nearnest Neighboor data association method
Along with the vehicle moving, at each observation time, it will take at most 2-12 observations, with bearing and range. It will incorporate observations and odometry input and update the its positions and landmarks' positions after "incK" steps.

Define graph procedure:
To be more precise, the initial position and initial covariance are set by value.insert(initPose), and graph.add(PriorFactorPose2(noiseModels.priorNoise)), which are provided by the gtsam library. They should be defined in global coordinate.
The odometry are inserted into the graph by graph.add(BetweenFactorPose2(lastPoseIndex, thisPoseIndex,odometry,noiseModel))
The observation are inserted into the graph by graph.add(BearingRangeFactor2D(robPoseIndex, landmarkIndex, bearing, range, noiseModels))
To solve the least square problem, by firstly define the Optimizer like LevenbergMarquardtOptimizer, then solve it by calling optimizer optimize function. After solving each batch, we will collect another batch of newFactors and initial values of nodes, so that we will 
reinitialize the graphs and values correspondingly 
There is a handon tutorial: https://borg.cc.gatech.edu/sites/edu.borg/files/downloads/gtsam.pdf

Data Association methods:
1.nearest neighboor without constraint. the cost will be defined as the Euclidean distance between the virtual predicted landmark based on the rob pose and measurement with the each actual ever seen landmarks. If the closest distance falls below the threshold, we will assign the closest landmark to that measurement. Otherwise, we will initiate a landmark for this measurement.

2.Maximum Likelihood(Mahalanobis nearest neighboor) without constraint. We will first project the landmark uncertainties and rob pose uncertainties into the measurement space, and calculate the Mahalanobis distance between the actual measurement with the predicted measurement of actual landmark. Same, If the closest distance falls below the threshold, we will assign the closest landmark to that measurement. Otherwise, we will initiate a landmark for this measurement.

3.Constraint data association:
We enforce that at most one third measurement can initiate a new landmark to enhance the probability to get loop closure. We also enforce the mutual exclusive assignments to seek an joint compatible data association solutions. It is an linear programming problem with Euclidean cost or Mahalanobis cost, with the minimum total cost as the objective. We adopt the online library to solve this problem
https://www.mathworks.com/matlabcentral/fileexchange/26836-lapjv-jonker-volgenant-algorithm-for-linear-assignment-problem-v3-0


Notes
In order to run the demo, please change the path to gtsam_toolbox appropriately.




 
