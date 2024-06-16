% Single Candidate Optimizer...

function [gbest,gbestval,fitcount,t_c]=SCO_func(fhd,Dimension,ps,Max_iter,VRmin,VRmax,X_suru,varargin)
rand('state',sum(100*clock));
D=Dimension;
gbest_position=zeros(Max_iter,D);


if length(VRmin)==1
    VRmin=repmat(VRmin,1,D); 
    VRmax=repmat(VRmax,1,D);
end

Range=VRmax-VRmin;

% % % S=VRmin+(VRmax-VRmin).*rand(ps,D);
% % % if D==2
% % %     S=X_suru;
% % % end

S=X_suru;

gbestval=feval(fhd,S',varargin{:});
gbest=S;
gbest_position(1,:)=gbest;
fitcount=ps;

POO=0; %% Intial counter to count unsuccessful fitness improvements
m=5; %% number of unsuccessful attempts to improve the fitness
alpha=round(Max_iter/3); % number of function evaluations in the First phase
b=2.4;
P=0;

% Distance = abs(pbest(1,1)-pbest(:,1)); % Analiz icin eklendi...
% DistanceSum(1) = sum(Distance)/(ps-1);

% 2. Search History Analysis
if D==2 % 2 boyut icin yapilacak...
	refresh = 100;
	figure(11)

% initialize
        x1 = linspace(VRmin(1), VRmax(1), 101);
        x2 = linspace(VRmin(2), VRmax(2), 101);
        x3 = zeros(length(x1), length(x2));
	for i = 1:length(x1)
		for j = 1:length(x2)
			x3(i, j) = feval(fhd,[x1(i);x2(j)],varargin{:});
		end
    end
    
    
% titles, labels, legend
	str = sprintf('Search History of FN%d',varargin{:});
	xlabel('x_1'); ylabel('x_2');  title(str);
	contour(x1', x2', x3'); hold on;
	plot(S(:,1),S(:,2),'bs','MarkerSize',8,'MarkerFaceColor','k');
	drawnow
	% plot(gbest(1),gbest(2),'k*','MarkerSize',8);
    figure(13)
    contour(x1', x2', x3'); hold on; 
    plot(gbest(1),gbest(2),'bs','MarkerSize',8,'MarkerFaceColor','k');
    str = sprintf('Trajectory of FN%d',varargin{:});
	xlabel('x_1'); ylabel('x_2');  title(str);
end
t_c=[];t=1;
    while fitcount<Max_iter %% && abs(gbestval-varargin{:}*100) > 1e-8
		w(t) = exp(-(b*t/Max_iter)^b); %% Equation (3) in the paper 
		if t>alpha
			if sum(P)==0       %% Counting the number of unsuccessful fitness improvements
				POO=1+POO;      %% Counter to count unsuccessful fitness improvements
			end
		end
		K=rand;
		for j = 1:D
			EE= w(t)*K*Range(j);
			if t<alpha 
				if rand<0.5
					x(j) = S(j)+(w(t)*abs(S(j)));
                else                        %% Equation (2) in the paper 
					x(j) = S(j)-(w(t)*abs(S(j)));
                end 
            else
				if POO==m
					POO=0;      %% Reset counter
					if rand<0.5
						x(j) = S(j)+rand*Range(j); 
					else                           %% Equation (5) in the paper
						x(j) = S(j)-rand*Range(j);
					end 
				else
					if rand<0.5
						x(j) = S(j)+EE; 
					else                           %% Equation (4) in the paper
						x(j) = S(j)-EE;
					end   
				end
			end
    
    %% Check if a dimension of the candidate solution goes out of boundaries
       if x(j)>VRmax(j)
           x(j)=S(j);
       end              %% Equation (6) in the paper
        if x(j)<VRmin(j)
           x(j)=S(j);
        end
    end
    
 %% Evaluate the fitness of the newly generated candidate solution
     F(t)=feval(fhd,x',varargin{:});
	 fitcount=fitcount+ps;
     if F(t)<gbestval
        gbestval=F(t);
		S=x;
        P=1;
    else
        P=0;
    end
	gbest=S;
	
		if fitcount==.01*Max_iter||fitcount==.02*Max_iter||fitcount==.03*Max_iter||fitcount==.05*Max_iter....
			||fitcount==.1*Max_iter||fitcount==.2*Max_iter||fitcount==.3*Max_iter||fitcount==.4*Max_iter||fitcount==.5*Max_iter....
			||fitcount==.6*Max_iter||fitcount==.7*Max_iter||fitcount==.8*Max_iter||fitcount==.9*Max_iter||fitcount==Max_iter....
				t_c=[t_c;abs(gbestval-varargin{:}*100)];              
        end	
	
	t=t+1;
	gbest_position(t,:)=gbest;
	% 2. Search History Analysis
		if D==2 && (rem(fitcount/ps,refresh) == 0)% 2 boyut icin yapilacak...
			figure(11)
			plot(x(:,1),x(:,2),'bs','MarkerSize',8,'MarkerFaceColor','k');
            xlabel('x_1'); ylabel('x_2'); title(str);
			drawnow
        end       	
    end
		
		if D==2
            figure(11)
            plot(gbest(1),gbest(2),'bs','MarkerSize',14,'MarkerFaceColor','b');
			% 3. Trajectory Analysis
			figure(12) 
% % % 			str = sprintf('Trajectory of Elite FN%d',varargin{:});
            subplot(211)
            plot(ps:ps:fitcount,gbest_position(:,1),'b','LineWidth',2);
			xlabel('Iteration'); ylabel('x_1 position');
            subplot(212)
            plot(ps:ps:fitcount,gbest_position(:,2),'b','LineWidth',2);
			xlabel('Iteration'); ylabel('x_2 position');
            
% % %             subplot(2,2,[1 3]);
            figure(13)
% % %             contour(x1', x2', x3'); 
            hold on; 
            plot(gbest_position(:,1)',gbest_position(:,2)','LineWidth',2,'Color','b');
            plot(gbest_position(end,1),gbest_position(end,2),'bd','MarkerSize',8,'MarkerFaceColor','k');

% % % 			plot(gbest_position(:,1), gbest_position(:,2), 'kd', 'MarkerSize', 7);
% % %             hold on 
% % % 			plot(gbest_position(end,1), gbest_position(end,2), 'rd', 'MarkerSize', 7);
% % % 			str = sprintf('Trajectory of Elite FN%d',varargin{:});
% % % 			xlabel('x_1'); ylabel('x_2');  title(str);
% % % 			subplot(2,2,2);  
% % % 			plot(ps:ps:fitcount,gbest_position(:,1),'b');
% % % 			xlabel('Iteration'); ylabel('x_1 position');
% % % 			subplot(2,2,4);  
% % % 			plot(ps:ps:fitcount,gbest_position(:,2),'b');
% % % 			xlabel('Iteration'); ylabel('x_2 position');
            drawnow
% % % 			% 4. Average Distance Analysis
% % % 			figure(13); 
% % % 			plot(DistanceSum,'b-');
% % % 			xlabel('iteration');
% % % 			ylabel('Average distance');
% % % 			str = sprintf('Average Distance of FN%d',varargin{:});
% % %             title(str);
% % %             drawnow
% % % 		end
end