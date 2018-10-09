% find initial config 
% Note: the initial values for the params are from the CAD model (i.e., angles and dimensions)
% 12/22/2017 Writen by Alirza Ramezani, Pasadena, CA
clear all; close all; clc;

               
%param & and initial guess for the angles (from CAD model) 
func_params;

% x_sol:
% alpha1
% alpha2
% alpha3
% alpha4
% alpha5
% alpha6
[x_sol]=func_sol_mechanism(x0,q1,q2,R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4);
alpha1 = x_sol(1);
alpha2 = x_sol(2);
alpha3 = x_sol(3);
alpha4 = x_sol(4);
alpha5 = x_sol(5);
alpha6 = x_sol(6);



disp('computed unkowns are ...');
disp(['alpha1 = ',num2str(alpha1*180/pi)]);
disp(['alpha2 = ',num2str(alpha2*180/pi)]);
disp(['alpha3 = ',num2str(alpha3*180/pi)]);
disp(['alpha4 = ',num2str(alpha4*180/pi)]);
disp(['alpha5 = ',num2str(alpha5*180/pi)]);
disp(['alpha6 = ',num2str(alpha6*180/pi)]);
% 
% computing joint positions ...
disp('computing fwd kinematics terms ...');
% computing fwd kinematics terms ...
% Note: these are symbolically generated...
% x:
% alpha1
% alpha2
% alpha3
% alpha4
% alpha5
% alpha6
% motor angles:
% q1
% q2
% params:
% R1 R2 R3 Rc1 Rc2 Rc3 e1 e2 e3 e4

% output:
% p1, p2 ,p3, ... p9
% v1, v2 ,v3, ... v9
[p1,p2,p3,p4,p5,p6,p7,p8,p9,v1,v2,v3,v4,v5,v6,v7,v8,v9] = ...
    func_leg_fwd_kin(x_sol, q1, q2, R1, R2, R3, Rc1, Rc2, Rc3, e1, e2, e3, e4);

% coordinates of the virtual leg
p_h = (p1+p2)./2;
p_e = p9;

%   
% plot leg mechanism
scrsz = get(groot,'ScreenSize');
fh = figure('Name','LEO leg mechanism (side view)',...
                'Renderer','opengl',...
                    'GraphicsSmoothing','on');
ah = axes('Box','on',...
            'XGrid','off',...
                'YGrid','off',...
                        'DataAspectRatio',[1,1,1],...
                              'PlotBoxAspectRatio',[1,1,1],...
                                    'Parent',fh);
                                xlabel(ah,'[m]');
                                ylabel(ah,'[m]');
                                hold(ah,'on');

% plot joints ...                                
plot(ah,p1(1),p1(2),'ko');
plot(ah,p2(1),p2(2),'ko');
plot(ah,p3(1),p3(2),'ko');
plot(ah,p4(1),p4(2),'ko');
plot(ah,p5(1),p5(2),'ko');
plot(ah,p6(1),p6(2),'ko');
plot(ah,p7(1),p7(2),'ko');
plot(ah,p8(1),p8(2),'ko');
plot(ah,p9(1),p9(2),'ko');

% plot links ...
% first loop
line(ah,[p1(1),p2(1)],[p1(2),p2(2)],'Color','b');  
line(ah,[p2(1),p3(1)],[p2(2),p3(2)],'Color','b');  
line(ah,[p3(1),p4(1)],[p3(2),p4(2)],'Color','b');  
line(ah,[p4(1),p5(1)],[p4(2),p5(2)],'Color','b');  
line(ah,[p5(1),p1(1)],[p5(2),p1(2)],'Color','b');  
line(ah,[p1(1),p3(1)],[p1(2),p3(2)],'Color','r','LineStyle','--');  
line(ah,[p3(1),p5(1)],[p3(2),p5(2)],'Color','r','LineStyle','--');  

% second loop
line(ah,[p6(1),p4(1)],[p6(2),p4(2)],'Color','b');  
line(ah,[p7(1),p6(1)],[p7(2),p6(2)],'Color','b');  
line(ah,[p8(1),p7(1)],[p8(2),p7(2)],'Color','b');
line(ah,[p4(1),p8(1)],[p4(2),p8(2)],'Color','b');  
line(ah,[p9(1),p7(1)],[p9(2),p7(2)],'Color','b');  
line(ah,[p6(1),p8(1)],[p6(2),p8(2)],'Color','r','LineStyle','--');  

% virtual leg
line(ah,[p_h(1),p_e(1)],[p_h(2),p_e(2)],'Color','g','LineStyle','--');  



          