function [delta_x, delta_y, delta_theta] = control2odom(model, lastPose, alpha, vc,T)
global Param

%u = [vc, Data.Control.alpha(i)]';
delta=model(Param.a,minimizedAngle(alpha),Param.b,Param.L,minimizedAngle(lastPose.theta),T,vc);
delta_x = delta(1);
delta_y = delta(2);
delta_theta = minimizedAngle(delta(3));
end