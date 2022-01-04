function [Leader_score,Leader_pos,WOACNVG]=WOA_JOS(N,T,maxRun,maxFE,BFid,nD,fhd,Jr)

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

lb=-100*ones(1,nD); ub=100*ones(1,nD);
dim = nD;

disp('WOA-JOS is now tackling your problem')

WOACNVG=zeros(maxRun+1,T);

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
    
    % initialize position vector and score for the leader
    Leader_pos=zeros(1,dim);
    Leader_score=inf; %change this to -inf for maximization problems
    
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
            %fitness=fobj(X(i,:));
            
            fitness=feval(fhd,X(i,:)',BFid);
            nFE = nFE + 1;
            if nFE > maxFE; break; end
            
            % Update the leader
            if fitness<Leader_score % Change this to > for maximization problem
                Leader_score=fitness; % Update alpha
                Leader_pos=X(i,:);
                Row_Id = i;
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
        WOACNVG(1,t) = nFE; %nFE
        WOACNVG(run+1,t) = Leader_score;
        
        a=2-t*((2)/T); % a decreases linearly fron 2 to 0
        % Oppose the least fitness elements
        threshold=a;
        
        X=corOppose2(X,ub,lb,upper,lower,dim,threshold,Row_Id);
        % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12) in WOA paper

        a2=-1+t*((-1)/T);
        
        % Update the Position of search agents
        for i=1:size(X,1)
            r1=rand(); % r1 is a random number in [0,1] in WOA paper
            r2=rand(); % r2 is a random number in [0,1] in WOA paper
            
            A=2*a*r1-a;  % Eq. (2.3) in the paper of WOA
            C=2*r2;      % Eq. (2.4) in the paper of WOA
            
            
            b=1;               %  parameters in Eq. (2.5) in WOA paper
            l=(a2-1)*rand+1;   %  parameters in Eq. (2.5) in WOA paper
            
            p = rand();        % p in Eq. (2.6) in WOA paper
            
            for j=1:size(X,2)
                
                if p<0.5
                    if abs(A)>=1
                        rand_leader_index = floor(N*rand()+1);
                        X_rand = X(rand_leader_index, :);
                        D_X_rand=abs(C*X_rand(j)-X(i,j)); % Eq. (2.7) in WOA paper
                        X(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8) in WOA paper
                        
                    elseif abs(A)<1
                        D_Leader=abs(C*Leader_pos(j)-X(i,j)); % Eq. (2.1) in WOA paper
                        X(i,j)=Leader_pos(j)-A*D_Leader;      % Eq. (2.2) in WOA paper
                    end
                    
                elseif p>=0.5
                    
                    distance2Leader=abs(Leader_pos(j)-X(i,j));
                    % Eq. (2.5) in WOA paper
                    X(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
                end
            end
        end
        
        if rand < Jr && nFE+size(X,1) < maxFE
            
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



