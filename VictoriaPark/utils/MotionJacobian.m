function J=MotionJacobian(JacobianModel, vc_,alpha_,phi_,t_)
global Param
J=JacobianModel(Param.a,minimizedAngle(alpha_),Param.b,Param.L,minimizedAngle(phi_),t_,vc_);
end
