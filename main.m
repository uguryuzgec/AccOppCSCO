clear 
close all
clc
% mex cec14_func.cpp -DWINDOWS
% 17-22 hibrid func ve 29-30...
% 1-3 unimodal func.
% 4-16 multimodal func.
% 23-28 composition func.
func_num=25; % function number
runs=1; % number of runs
D=2; % dimension
Xmin=-100;
Xmax=100;
pop_size=1; % it is 1 for single candidate optimizer
iter_max=5000;
fhd=str2func('cec14_func');
empty_solution.cost=[];
empty_solution.position=[];
empty_solution.t=[];
solution=repmat(empty_solution,func_num,runs);
X_suru=Xmin+(Xmax-Xmin).*rand(pop_size,D);

plot_function
drawnow

% Optimization with SCO...
[gbest,gbestval,FES,t]= SCO_func(fhd,D,pop_size,iter_max,Xmin,Xmax,X_suru,func_num);
solution(func_num,1).position = gbest;
solution(func_num,1).cost = gbestval-func_num*100;
solution(func_num,1).t = t;
fprintf('Func no: %d -> %d. run : best error = %1.2e\n',func_num,runs,solution(func_num,1).cost);


fonk_numara = func_num;
FESindex = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]*FES;

% convergence curve for selected problem...
figure (2)
semilogy(FESindex,solution(fonk_numara,1).t,'-db','LineWidth',2);
xlabel('Function Evaluations');
ylabel('Error Value');
str = sprintf('Convergence Analysis of FN%d',fonk_numara);
title(str);

% Optimization with AccOppCSCO...

[i_gbest,i_gbestval,i_FES,i_t]= AccOppCSCO_func(fhd,D,pop_size,iter_max,Xmin,Xmax,X_suru,func_num);
solution1(func_num,1).position = i_gbest;
solution1(func_num,1).cost = i_gbestval-func_num*100;
solution1(func_num,1).t = i_t;
fprintf('\n---------------------------------------------------------------\n');
fprintf('Optimization with AccOppCSCO\n ');
fprintf('Func no: %d -> %d. run : best error = %1.2e\n',func_num,runs,solution1(func_num,1).cost);

% convergence curve for selected problem...
figure (2)
hold on
semilogy(FESindex,solution1(fonk_numara,1).t,'-sr','LineWidth',2);
legend('SCO','AccOppCSCO')