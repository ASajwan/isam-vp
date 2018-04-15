function [g] = createMotionModel
syms vc l alp phi a b t 

%This motion mode refers to Vicotrial DataSet motion model
%http://www-personal.acfr.usyd.edu.au/nebot/experimental_data/modeling_info/Ute_modeling_info.htm

rot=[cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1];
del=[t*(vc*cos(phi)-(vc/l)*tan(alp)*(a*sin(phi)+b*cos(phi)));t*(vc*sin(phi)+(vc/l)*tan(alp)*(a*cos(phi)-b*sin(phi)));t*(vc/l)*tan(alp)];
measurement=rot*del;
g = matlabFunction(measurement); %Inputs  g(a,alpha,b,l,phi,t,vc)
end

