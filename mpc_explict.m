clc;
clf;
close all;
format long;
A=[-4.59 -0.94;1.52 -4.44];
Bu=[2.29;-0.76];
Bd=[2.30;10.67];
C=[1 0;0 1];
T=0.02;
[Ga,Hu]=c2d(A,Bu,T);
[Gb,Hd]=c2d(A,Bd,T);
A=Ga;
B=[Hu Hd];
C=eye(2,2);
D = zeros(2,2);
model = LTISystem('A', A, 'B', B, 'C', C, 'D', D);

model.x.min = [-1; -0.85];
model.x.max = [1; 0.85];
P=[1,0;0,1];
model.x.penalty = QuadFunction(P);
model.u.max = [inf;0.1];
model.u.min = [-inf;0.1];
model.u.penalty = QuadFunction(5.0*eye(model.nu));

horizon = 10;
ctrl = MPCController(model, horizon);

disp('Generating the explicit solution:');
expctrl = ctrl.toExplicit()

% Compare optimal solutions
x0 = [0; 0];
disp('Optimal control input obtained by evaluating the explicit solution:');
Uexp = expctrl.evaluate(x0)

disp('Optimal control input obtained by solving the optimization problem on-line:');
Uonl = ctrl.evaluate(x0)

close all
% plot the explicit optimizer


% plot the value function
figure
expctrl.cost.fplot();
title('Explicit cost function');

% plot the regions
figure
expctrl.partition.plot()
title('Regions of the polyhedral partition');

loop = ClosedLoop(ctrl, model);
T=0.02;
x0 = [0.0;0.0];
N_sim = 100;
data = loop.simulate(x0, N_sim);
figure(3)
plot((1:1:N_sim)*T,data.U(1,:),'-','LineWidth',3);
xlabel('t/s');
title('Control Inpout of u');
legend('u')
hh=findobj('tag','legend') %|
set(hh,'fontsize',10) %| 设置legend字号大小
set(hh,'LineWidth',1.0); %| 设置图形线宽
set(gca,'linewidth',1.5) %| 设置图形外边框的线宽1.5
set(gca,'box','off') %| 去图形外筐
figure(4)
plot((1:1:N_sim)*T,data.Y(1,:), '-','LineWidth',3);
title('\beta');
legend('\beta');
xlabel('t/s');
hh=findobj('tag','legend') %|
set(hh,'fontsize',10) %| 设置legend字号大小
set(hh,'LineWidth',1.0); %| 设置图形线宽
set(gca,'linewidth',1.5) %| 设置图形外边框的线宽1.5
set(gca,'box','off') %| 去图形外筐
figure(5)
plot((1:1:N_sim)*T,data.Y(2,:), '-','LineWidth',3);
title('r');
legend('r');
xlabel('t/s');
hh=findobj('tag','legend') %|
set(hh,'fontsize',10) %| 设置legend字号大小
set(hh,'LineWidth',1.0); %| 设置图形线宽
set(gca,'linewidth',1.5) %| 设置图形外边框的线宽1.5
set(gca,'box','off') %| 去图形外筐
