% Accelerated Opposition Learning based Chaotic Single Candidate Optimization Algorithm (AccOppCSCO) 
% This code was written by Ugur Yuzgec on January 2024.
% ---------------------------------------------------------------------------------------------------------------

function [gbest,gbestval,fitcount,t_c]=AccOppCSCO_func(fhd,Dimension,ps,Max_iter,VRmin,VRmax,X_suru,varargin)
rand('state',sum(100*clock));
D=Dimension;
gbest_position=zeros(Max_iter,D);

%-------------------------------------------------
% Opposition parameter...
jump_rate = 0.5; % jumping rate constant.
jump_rate_sayac = 0;
e_max=1;
e_min=1e-5;
%-------------------------------------------------

if length(VRmin)==1
    VRmin=repmat(VRmin,1,D); 
    VRmax=repmat(VRmax,1,D);
end

Range = VRmax-VRmin;
lu = [VRmin;VRmax];

% % % S=VRmin+(VRmax-VRmin).*rand(ps,D);
% % % if D==2
% % %     S=X_suru;
% % % end

S=X_suru; % First candidate is taken from outside of the function as parameter

% Opposition based learning initiliazing...
OppS = (VRmin + VRmax) - S;

gbestval=feval(fhd,S',varargin{:});
gbestval_Opp=feval(fhd,OppS',varargin{:});

if gbestval_Opp < gbestval
   S = OppS;
   gbestval = gbestval_Opp;
end

gbest=S;
gbest_position(1,:)=gbest;
fitcount=ps;

POO=0; %% Intial counter to count unsuccessful fitness improvements
m=5; %% number of unsuccessful attempts to improve the fitness
alpha=round(Max_iter/3); % number of function evaluations in the First phase
b=2.4;
P=0;

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
	plot(S(:,1),S(:,2),'ro','MarkerSize',8,'MarkerFaceColor','k');
	drawnow
	% plot(gbest(1),gbest(2),'k*','MarkerSize',8);
    figure(13)
% % %     contour(x1', x2', x3'); 
    hold on; 
    plot(gbest(1),gbest(2),'ro','MarkerSize',8,'MarkerFaceColor','k');
end
t_c=[];t=1;
x=VRmin+(VRmax-VRmin).*rand(ps,D);
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
                    if rand<0.5
                        x(j) = S(j)+(w(t)*sin((2*pi)*rand)*abs(S(j))); %% updated Equation (2) in the paper
                    else                         
                        x(j) = S(j)+(w(t)*cos((2*pi)*rand)*abs(S(j))); 
                    end 
                else
                  	if rand<0.5
                    	x(j) = S(j)+(w(t)*abs(S(j)));
                    else                                %% Equation (2) in the paper 
                        x(j) = S(j)-(w(t)*abs(S(j)));
                    end 
                end
            else
				if POO==m
					POO=0;      %% Reset counter
                    r1=rand;
					% mutation operator based on chaotic functions
                    if r1>0.7
                        x(j)=cauchy(S(j)); % cauchy mutation...
                    elseif r1<0.7 && r1>0.3
                        x(j)=gauss(S(j)); % gauss mutation... 
                    else
                        x(j)=levy(S(j)); % levy mutation...
                    end
				else
					if rand<0.5
                        if rand<0.5
                            x(j) = S(j)+EE*sin((2*pi)*rand); 
                        else                           %% updated Equation (4) in the paper
                            x(j) = S(j)+EE*cos((2*pi)*rand);
                        end
                    else
                        if rand<0.5
                            x(j) = S(j)+EE; 
                        else                           %% Equation (4) in the paper
                            x(j) = S(j)-EE;
                        end
                    end
				end
			end
		end % for j...
       
		 %% Check if a dimension of the candidate solution goes out of boundaries
		%------------------------------------------------------------------------
		x = boundConstraint(x,S,lu); % new boundary mechanism...
		%------------------------------------------------------------------------ 
 
		e_ac=e_max-t*(e_max-e_min)/Max_iter; % Accelerating mechanism of opposition learning
	 
	 %% Evaluate the fitness of the newly generated candidate solution
     F(t)=feval(fhd,x',varargin{:});
	 fitcount=fitcount+ps;
	 % opposition based generation jumping...........
	 %---------------------------------------------%
	 if (rand() < jump_rate),
	 % fprintf(1,'\n opposition based generation jumping......');
		jump_rate_sayac = jump_rate_sayac + 1;
		% Oppx = rand(1,dim)*(lb + ub-2*x) - x;
		Oppx = e_ac.*(rand(1,D).*VRmin + rand(1,D).*VRmax) - x; % NEW FORM...
        %% Check if a dimension of the candidate solution goes out of boundaries
		%------------------------------------------------------------------------
		Oppx = boundConstraint(Oppx,S,lu); % new boundary mechanism...
		%------------------------------------------------------------------------ 
		F_opp = feval(fhd,Oppx',varargin{:});
		if F_opp < F(t)
			x = Oppx;
			F(t) = F_opp;
		end

	end
	%---------------------------------------------%
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
		if D==2 && (rem(fitcount/ps,refresh) == 0)% only be done for two dimension...
			figure(11)
			plot(x(:,1),x(:,2),'ro','MarkerSize',8,'MarkerFaceColor','k');
			xlabel('x_1'); ylabel('x_2'); title(str);
			drawnow
        end         
    end
		
		if D==2
            figure(11)
            plot(gbest(1),gbest(2),'ro','MarkerSize',10,'MarkerFaceColor','r');
			% 3. Trajectory Analysis
            figure(12) 
            str = sprintf('Trajectory of Elite FN%d',varargin{:});
            subplot(211)
            hold on
            plot(ps:ps:fitcount,gbest_position(:,1),'r','LineWidth',2);
			xlabel('Iteration'); ylabel('x_1 position');title(str);
            legend('SCO','AccOppCSCO')
            subplot(212)
            hold on
            plot(ps:ps:fitcount,gbest_position(:,2),'r','LineWidth',2);
			xlabel('Iteration'); ylabel('x_2 position');
            legend('SCO','AccOppCSCO')
            
            figure(13)
            hold on; 
            plot(gbest_position(:,1)',gbest_position(:,2)','LineWidth',2,'Color','r');
            plot(gbest_position(end,1),gbest_position(end,2),'rd','MarkerSize',8,'MarkerFaceColor','k');

            drawnow

end