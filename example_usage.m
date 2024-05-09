clc;clear;close all

% This file is created for purpose on showing an example usage of unconOPT 
% nonlinear unconstrained optimization solver usage. After running this
% cell you can run any method cell.

% --------------------------------------------
% The solver can be use below search method
% 'steepest-descend'
% 'classical-conjugate'
% 'hestenes-stiefel-conjugate'
% 'fletcher-reeves-conjugate'
% 'polak-ribiere-conjugate'
% 'modified-newton'
% 'davidon-fletcher-powell'
% 'broyden-fletcher-goldfarb-shanno'
% --------------------------------------------


% First define the problem f, and x initial
% There are many f function is written below. If you want to solve any
% below cost function you can uncomment it and comment others.

f = @(x) 0.1*x(1)^2 + x(2)^2 -10;
%f = @(x) (x(1)^2 + 20*x(2)^2 - 6*x(1)*x(2));
%f = @(x) x(1)^2 + x(2)^2 -2*x(1)*x(2);
%f = @(x) x(1)^2 + 2*x(2)^2 + 2*x(3)^2 + 2*x(1)*x(2) + 2*x(2)*x(3);
%f = @(x) 3*x(1)^2 + 2*x(1)*x(2) + 2*x(2)^2 + 7
%f = @(x) 10*x(1)^4 - 20*x(1)^2*x(2) + 10*x(2)^2 + x(1)^2 - 2*x(1) + 5
%f = @(x) 5*x(1)^2 + 2*x(1)*x(2) + x(2)^2 + 7;

% Specify x initial value to search starting point
x_init = [1 ; 1];

% Each cell below use different method. You can run each cell seperetaly.
%%
method = 'steepest-descend';
solver = unconOPT(x_init,f,method);
% solver.i
% solver.x_opt(end,:)
% solver.f_opt(end)

figure
subplot(1,2,1)
plot([1:solver.i],solver.f_opt,'--r')
hold on
plot([1:solver.i],solver.f_opt,'.b','MarkerSize',20)
xlabel('iteration')
ylabel('function evaluation')
title(method)

subplot(1,2,2)
for i=1:length(solver.x_opt(1,:))
    plot([1:solver.i],solver.x_opt(:,i),'.','MarkerSize',20)
    hold on
    plot([1:solver.i],solver.x_opt(:,i),'--')
    hold on
    
end
xlabel('iteration')
ylabel('design variable value')

legend_list = [];

for i=1:length(solver.x_opt(1,:))
    legend_list = [legend_list, ['x_',num2str(i)] ,"-"];
end
legend(legend_list)
title(method)


disp('-------------------------------')
disp([method, ' Method Solution'])
fprintf('Optimal point is [%s]\n',num2str(solver.x_opt(end,:)))
fprintf('Optimal cost function value is [%s]\n',num2str(solver.f_opt(end)))
disp('-------------------------------')

%%
method = 'classical-conjugate';
solver = unconOPT(x_init,f,method);

figure
subplot(1,2,1)
plot([1:solver.i],solver.f_opt,'--r')
hold on
plot([1:solver.i],solver.f_opt,'.b','MarkerSize',20)
xlabel('iteration')
ylabel('function evaluation')
title(method)

subplot(1,2,2)
for i=1:length(solver.x_opt(1,:))
    plot([1:solver.i],solver.x_opt(:,i),'.','MarkerSize',20)
    hold on
    plot([1:solver.i],solver.x_opt(:,i),'--')
    hold on
    
end
xlabel('iteration')
ylabel('design variable value')

legend_list = [];

for i=1:length(solver.x_opt(1,:))
    legend_list = [legend_list, ['x_',num2str(i)] ,"-"];
end
legend(legend_list)
title(method)

disp('-------------------------------')
disp([method, ' Method Solution'])
fprintf('Optimal point is [%s]\n',num2str(solver.x_opt(end,:)))
fprintf('Optimal cost function value is [%s]\n',num2str(solver.f_opt(end)))
disp('-------------------------------')


%%
method = 'hestenes-stiefel-conjugate';
solver = unconOPT(x_init,f,method);

figure
subplot(1,2,1)
plot([1:solver.i],solver.f_opt,'--r')
hold on
plot([1:solver.i],solver.f_opt,'.b','MarkerSize',20)
xlabel('iteration')
ylabel('function evaluation')
title(method)

subplot(1,2,2)
for i=1:length(solver.x_opt(1,:))
    plot([1:solver.i],solver.x_opt(:,i),'.','MarkerSize',20)
    hold on
    plot([1:solver.i],solver.x_opt(:,i),'--')
    hold on
    
end
xlabel('iteration')
ylabel('design variable value')

legend_list = [];

for i=1:length(solver.x_opt(1,:))
    legend_list = [legend_list, ['x_',num2str(i)] ,"-"];
end
legend(legend_list)
title(method)

disp('-------------------------------')
disp([method, ' Method Solution'])
fprintf('Optimal point is [%s]\n',num2str(solver.x_opt(end,:)))
fprintf('Optimal cost function value is [%s]\n',num2str(solver.f_opt(end)))
disp('-------------------------------')

%%
method = 'fletcher-reeves-conjugate';
solver = unconOPT(x_init,f,method);

figure
subplot(1,2,1)
plot([1:solver.i],solver.f_opt,'--r')
hold on
plot([1:solver.i],solver.f_opt,'.b','MarkerSize',20)
xlabel('iteration')
ylabel('function evaluation')
title(method)

subplot(1,2,2)
for i=1:length(solver.x_opt(1,:))
    plot([1:solver.i],solver.x_opt(:,i),'.','MarkerSize',20)
    hold on
    plot([1:solver.i],solver.x_opt(:,i),'--')
    hold on
    
end
xlabel('iteration')
ylabel('design variable value')

legend_list = [];

for i=1:length(solver.x_opt(1,:))
    legend_list = [legend_list, ['x_',num2str(i)] ,"-"];
end
legend(legend_list)
title(method)

disp('-------------------------------')
disp([method, ' Method Solution'])
fprintf('Optimal point is [%s]\n',num2str(solver.x_opt(end,:)))
fprintf('Optimal cost function value is [%s]\n',num2str(solver.f_opt(end)))
disp('-------------------------------')

%%
method = 'polak-ribiere-conjugate';
solver = unconOPT(x_init,f,method);

figure
subplot(1,2,1)
plot([1:solver.i],solver.f_opt,'--r')
hold on
plot([1:solver.i],solver.f_opt,'.b','MarkerSize',20)
xlabel('iteration')
ylabel('function evaluation')
title(method)

subplot(1,2,2)
for i=1:length(solver.x_opt(1,:))
    plot([1:solver.i],solver.x_opt(:,i),'.','MarkerSize',20)
    hold on
    plot([1:solver.i],solver.x_opt(:,i),'--')
    hold on
    
end
xlabel('iteration')
ylabel('design variable value')

legend_list = [];

for i=1:length(solver.x_opt(1,:))
    legend_list = [legend_list, ['x_',num2str(i)] ,"-"];
end
legend(legend_list)
title(method)

disp('-------------------------------')
disp([method, ' Method Solution'])
fprintf('Optimal point is [%s]\n',num2str(solver.x_opt(end,:)))
fprintf('Optimal cost function value is [%s]\n',num2str(solver.f_opt(end)))
disp('-------------------------------')

%%
method = 'modified-newton';
solver = unconOPT(x_init,f,method);

figure
subplot(1,2,1)
plot([1:solver.i],solver.f_opt,'--r')
hold on
plot([1:solver.i],solver.f_opt,'.b','MarkerSize',20)
xlabel('iteration')
ylabel('function evaluation')
title(method)

subplot(1,2,2)
for i=1:length(solver.x_opt(1,:))
    plot([1:solver.i],solver.x_opt(:,i),'.','MarkerSize',20)
    hold on
    plot([1:solver.i],solver.x_opt(:,i),'--')
    hold on
    
end
xlabel('iteration')
ylabel('design variable value')

legend_list = [];

for i=1:length(solver.x_opt(1,:))
    legend_list = [legend_list, ['x_',num2str(i)] ,"-"];
end
legend(legend_list)
title(method)

disp('-------------------------------')
disp([method, ' Method Solution'])
fprintf('Optimal point is [%s]\n',num2str(solver.x_opt(end,:)))
fprintf('Optimal cost function value is [%s]\n',num2str(solver.f_opt(end)))
disp('-------------------------------')

%%
method = 'davidon-fletcher-powell';
solver = unconOPT(x_init,f,method);

figure
subplot(1,2,1)
plot([1:solver.i],solver.f_opt,'--r')
hold on
plot([1:solver.i],solver.f_opt,'.b','MarkerSize',20)
xlabel('iteration')
ylabel('function evaluation')
title(method)

subplot(1,2,2)
for i=1:length(solver.x_opt(1,:))
    plot([1:solver.i],solver.x_opt(:,i),'.','MarkerSize',20)
    hold on
    plot([1:solver.i],solver.x_opt(:,i),'--')
    hold on
    
end
xlabel('iteration')
ylabel('design variable value')

legend_list = [];

for i=1:length(solver.x_opt(1,:))
    legend_list = [legend_list, ['x_',num2str(i)] ,"-"];
end
legend(legend_list)
title(method)

disp('-------------------------------')
disp([method, ' Method Solution'])
fprintf('Optimal point is [%s]\n',num2str(solver.x_opt(end,:)))
fprintf('Optimal cost function value is [%s]\n',num2str(solver.f_opt(end)))
disp('-------------------------------')

%%
method = 'broyden-fletcher-goldfarb-shanno';
solver = unconOPT(x_init,f,method);

figure
subplot(1,2,1)
plot([1:solver.i],solver.f_opt,'--r')
hold on
plot([1:solver.i],solver.f_opt,'.b','MarkerSize',20)
xlabel('iteration')
ylabel('function evaluation')
title(method)

subplot(1,2,2)
for i=1:length(solver.x_opt(1,:))
    plot([1:solver.i],solver.x_opt(:,i),'.','MarkerSize',20)
    hold on
    plot([1:solver.i],solver.x_opt(:,i),'--')
    hold on
    
end
xlabel('iteration')
ylabel('design variable value')

legend_list = [];

for i=1:length(solver.x_opt(1,:))
    legend_list = [legend_list, ['x_',num2str(i)] ,"-"];
end
legend(legend_list)
title(method)

disp('-------------------------------')
disp([method, ' Method Solution'])
fprintf('Optimal point is [%s]\n',num2str(solver.x_opt(end,:)))
fprintf('Optimal cost function value is [%s]\n',num2str(solver.f_opt(end)))
disp('-------------------------------')
