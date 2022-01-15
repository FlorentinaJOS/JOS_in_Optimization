function [Best_flame_score,Best_flame_pos,MFOJOS]=MFO_JOS(N,T,maxRun,maxFE,BFid,nD,fhd,Jr)

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

MFOJOS=zeros(maxRun+1,T);

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
    %Initialize the positions of moths
    Moth_pos=initialization(N,dim,ub,lb);
    
    % Modified from the original Dynamic Opposite (DO)
    for i=1:N
        OP(i,:)=((ub-lb).*rand(size(lb)))+lb-Moth_pos(i,:);
        RO(i,:) = rand*OP(i,:);
        DO(i,:) = Moth_pos(i,:) + rand*(RO(i,:)-Moth_pos(i,:));
        
        for j=1:dim
            if DO(i,j)<lb(1,j)
                DO(i,j)=lb(1,j);
            end
            if DO(i,j)>ub(1,j)
                DO(i,j)=ub(1,j);
            end
        end
    end
    Moth_pos=DO;
    
    nFE=0;
    Iter=0;
    Row_Id=1;
    while nFE<maxFE
        % Number of flames Eq. (3.14) in the paper MFO
        Flame_no=round(N-Iter*((N-1)/T));
        for i=1:size(Moth_pos,1)
            % Check if moths go out of the search spaceand bring it back
            Flag4ub=Moth_pos(i,:)>ub;
            Flag4lb=Moth_pos(i,:)<lb;
            Moth_pos(i,:)=(Moth_pos(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            % Calculate the fitness of moths
            Moth_fitness(1,i)=feval(fhd,Moth_pos(i,:)',BFid);
            nFE=nFE+1;
         
            if Iter==0
                % Sort the first population of moths
                [fitness_sorted I]=sort(Moth_fitness);
                sorted_population=Moth_pos(I,:);
                % Update the flames
                best_flames=sorted_population;
                best_flame_fitness=fitness_sorted;
                Row_Id=i;
            else
                % Sort the moths
                double_population=[previous_population;best_flames];
                double_fitness=[previous_fitness best_flame_fitness];
                
                [double_fitness_sorted, I]=sort(double_fitness);
                double_sorted_population=double_population(I,:);
                
                fitness_sorted=double_fitness_sorted(1:N);
                sorted_population=double_sorted_population(1:N,:);
                
                % Update the flames
                best_flames=sorted_population;
                best_flame_fitness=fitness_sorted;
                Row_Id=i;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %updating boundary for opposition after every iteration
            for x=1:size(Moth_pos,1)
                for y=1:size(Moth_pos,2)
                    if upper(1,y)<Moth_pos(x,y)
                        upper(1,y)=Moth_pos(x,y);
                    end
                    if lower(1,y)>Moth_pos(x,y)
                        lower(1,y)=Moth_pos(x,y);
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        % Update the position best flame obtained so far
        Best_flame_score=fitness_sorted(1);
        Best_flame_pos=sorted_population(1,:);
        
        previous_population=Moth_pos;
        previous_fitness=Moth_fitness;
        
        Iter=Iter+1;
        MFOJOS(1,Iter)=nFE;
        MFOJOS(run+1,Iter)=Best_flame_score;
        
        % a linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
        a=-1+Iter*((-1)/T);
        threshold=a;
        
        % Selective Leading Opposition (SLO) modified from the original Selective Opposition (SO)
        Moth_pos=corOppose2(Moth_pos,ub,lb,upper,lower,dim,threshold,Row_Id);
        
        for i=1:size(Moth_pos,1)
            for j=1:size(Moth_pos,2)
                if i<=Flame_no % Update the position of the moth with respect to its corresponsing flame
                    % D in Eq. (3.13) in MFO paper
                    distance_to_flame=abs(sorted_population(i,j)-Moth_pos(i,j));
                    b=1;
                    t=(a-1)*rand+1;
                    % Eq. (3.12) in MFO paper
                    Moth_pos(i,j)=distance_to_flame*exp(b.*t).*cos(t.*2*pi)+sorted_population(i,j);
                end
                
                if i>Flame_no % Upaate the position of the moth with respct to one flame
                    % Eq. (3.13) in MFO paper
                    distance_to_flame=abs(sorted_population(i,j)-Moth_pos(i,j));
                    b=1;
                    t=(a-1)*rand+1;
                    % Eq. (3.12) in MFO paper
                    Moth_pos(i,j)=distance_to_flame*exp(b.*t).*cos(t.*2*pi)+sorted_population(Flame_no,j);
                end
            end
        end
        
        % Modified from the original Dynamic Opposite (DO)
        if rand < Jr && nFE+N < maxFE
            
            for i=1:N
                OP(i,:)=((ub-lb).*rand(size(lb)))+lb-Moth_pos(i,:);
                RO(i,:) = rand*OP(i,:);
                DO(i,:) = Moth_pos(i,:) + rand*(RO(i,:)-Moth_pos(i,:));
                for j=1:dim
                    if DO(i,j)<lb(1,j)
                        DO(i,j)=lb(1,j);
                    end
                    if DO(i,j)>ub(1,j)
                        DO(i,j)=ub(1,j);
                    end
                end
            end
            Moth_pos=DO;
        end
    end
    toc
end
% display(['The best solution obtained by MFO is : ', num2str(Best_flame_score)]);
% display(['The best optimal value of the objective funciton found by MFO is : ', num2str(Best_flame_pos)]);
end

function X=initialization(SearchAgents_no,dim,ub,lb)
Boundary_no= size(ub,2); % numnber of boundaries
% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    X=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end
% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
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



