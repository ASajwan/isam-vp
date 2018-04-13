function JacobianModel = createMotionJacobian
syms vc l alp phi a b t 

%rot=[cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
del=[t*(vc*cos(phi)-(vc/l)*tan(alp)*(a*sin(phi)+b*cos(phi)));t*(vc*sin(phi)+(vc/l)*tan(alp)*(a*cos(phi)-b*sin(phi)));t*(vc/l)*tan(alp)];
%measurement=rot*del;
measurement= del;

H=jacobian(measurement,[vc, alp]);
JacobianModel = matlabFunction(H); %Inputs  g(a,alpha,b,l,phi,t,vc)
end

