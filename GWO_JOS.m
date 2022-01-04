function [Alpha_score,Alpha_pos,GWOJOS]=GWO_JOS(N,T,maxRun,maxFE,BFid,nD,fhd,Jr)

if nargin ~= 8
    N=30;               % Population size
    maxRun = 10;        % Maximum Run
    BFid = 1;           % Number id of benchmark function
    nD = 10;            % Number of dimensions
    maxFE = 10000*nD;   % Number of function evaluations
    Jr=0.25;            % Jumping Rate
    T=ceil(maxFE/N);    % Maximum number of iterations
    fhd=str2func('cec17_func');
end

lb = -100*ones(1,nD); ub = 100*ones(1,nD);
dim = nD;

disp('GWO-JOS is now tackling your problem')

GWOJOS=zeros(maxRun+1,T);

for run=1:maxRun
    tic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %initializing boundary for opposition
    Boundary_no= size(ub,2); % numnber of boundaries
    
    % If the boundaries of all variables are equal and user enter a signle
    % number for both ub and lb
    if Boundary_no==1
        for i=1:dim
            upper(1,i)=ub;
            lower(1,i)=lb;
        end
        % If each variable has a different lb and ub
    else
        for i=1:dim
            upper(1,i)=ub(i);
            lower(1,i)=lb(i);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % initialize alpha, beta, and delta_pos
    Alpha_pos=zeros(1,dim);
    Alpha_score=inf; %change this to -inf for maximization problems
    
    Beta_pos=zeros(1,dim);
    Beta_score=inf; %change this to -inf for maximization problems
    
    Delta_pos=zeros(1,dim);
    Delta_score=inf; %change this to -inf for maximization problems
    
    %Initialize the positions of search agents
    X=initialization(N,dim,ub,lb);
    
    for i=1:N
        OP(i,:)=((ub-lb).*rand(size(lb)))+lb-X(i,:);
        RO(i,:) = rand*OP(i,:);
        DO(i,:) = X(i,:) + rand*(RO(i,:)-X(i,:));
        
        for j=1:dim
            if DO(i,j)<lb(1,j)
                DO(i,j)=lb(1,j);
            end
            if DO(i,j)>ub(1,j)
                DO(i,j)=ub(1,j);
            end
        end
    end
    X=DO;
    
    t=0;% Loop counter
    nFE=0;
    Row_Id=1;
    while nFE<maxFE
        for i=1:size(X,1)
            
            % Return back the search agents that go beyond the boundaries of the search space
            Flag4ub=X(i,:)>ub;
            Flag4lb=X(i,:)<lb;
            X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            
            % Calculate objective function for each search agent
            fitness(i)=feval(fhd,X(i,:)',BFid);
            nFE = nFE + 1;
            if nFE > maxFE; break; end
            
            % Update Alpha, Beta, and Delta
            if fitness(i)<Alpha_score
                Alpha_score=fitness(i); % Update alpha
                Alpha_pos=X(i,:);
                Row_Id=i;
            end
            
            if fitness(i)>Alpha_score && fitness(i)<Beta_score
                Beta_score=fitness(i); % Update beta
                Beta_pos=X(i,:);
                Row_Id=i;
            end
            
            if fitness(i)>Alpha_score && fitness(i)>Beta_score && fitness(i)<Delta_score
                Delta_score=fitness(i); % Update delta
                Delta_pos=X(i,:);
                Row_Id=i;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %updating boundary for opposition after every iteration
            for x=1:size(X,1)
                for y=1:size(X,2)
                    if upper(1,y)<X(x,y)
                        upper(1,y)=X(x,y);
                    end
                    if lower(1,y)>X(x,y)
                        lower(1,y)=X(x,y);
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        t=t+1;
        GWOJOS(1,t) = nFE; %nFE
        GWOJOS(run+1,t) = Alpha_score;
        
        a=2-t*((2)/T); % a decreases linearly fron 2 to 0
        threshold=a;
        X=corOppose2(X,ub,lb,upper,lower,dim,threshold,Row_Id);
        
        % Update the Position of search agents including omegas
        for i=1:size(X,1)
            for j=1:size(X,2)
                
                r1=rand(); % r1 is a random number in [0,1]
                r2=rand(); % r2 is a random number in [0,1]
                
                A1=2*a*r1-a; % Equation (3.3)
                C1=2*r2; % Equation (3.4)
                
                D_alpha=abs(C1*Alpha_pos(j)-X(i,j)); % Equation (3.5)-part 1
                X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
                
                r1=rand();
                r2=rand();
                
                A2=2*a*r1-a; % Equation (3.3)
                C2=2*r2; % Equation (3.4)
                
                D_beta=abs(C2*Beta_pos(j)-X(i,j)); % Equation (3.5)-part 2
                X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2
                
                r1=rand();
                r2=rand();
                
                A3=2*a*r1-a; % Equation (3.3)
                C3=2*r2; % Equation (3.4)
                
                D_delta=abs(C3*Delta_pos(j)-X(i,j)); % Equation (3.5)-part 3
                X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3
                
                X(i,j)=(X1+X2+X3)/3;% Equation (3.7)
                
            end
        end
        
        if rand < Jr && nFE+N < maxFE
            
            for i=1:N
                OP(i,:)=((ub-lb).*rand(size(lb)))+lb-X(i,:);
                RO(i,:) = rand*OP(i,:);
                DO(i,:) = X(i,:) + rand*(RO(i,:)-X(i,:));
                for j=1:dim
                    if DO(i,j)<lb(1,j)
                        DO(i,j)=lb(1,j);
                    end
                    if DO(i,j)>ub(1,j)
                        DO(i,j)=ub(1,j);
                    end
                end
            end
            X=DO;
        end
    end
    toc
end
% display(['The best solution obtained by GWO is : ', num2str(Best_pos)]);
% display(['The best optimal value of the objective funciton found by GWO is : ', num2str(Best_score)]);
end

function Positions=initialization(N,dim,ub,lb)
Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    Positions=rand(N,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(N,1).*(ub_i-lb_i)+lb_i;
    end
end
end

function [Positions]=corOppose2(Positions,ub,lb,upper,lower,dim,threshold,Row_Id)

[n,b] = size(Positions);
for i=1:n(1) 
    
    if i ~= Row_Id
        
        sum=0;
        greater=[];
        less=[];
        x=1;z=1;y=1;
        
        for j=1:b
            d(x)=abs(Positions(Row_Id,j)-Positions(i,j));
            if d(x)<threshold
                greater(y)=j;
                y=y+1;
            else
                less(z)=j;
                z=z+1;
            end
            sum=sum+d(x)*d(x);
            x=x+1;
        end
        %     sum
        src=1-(double(6*sum))/(double(n(1)*(n(1)*n(1)-1)));
        %     src
        if src<=0
            if size(greater)<size(less)
                %             for j=1:size(less)
                %                 dim=less(j);
                %                 Positions(i,dim)=ub(dim)+lb(dim)-Positions(i,dim);
                %             end
            else
                for j=1:size(greater)
                    dim=greater(j);
                    Positions(i,dim)=(upper(1,dim)+lower(1,dim)-Positions(i,dim));
                end
            end
        end
    end
end
end




