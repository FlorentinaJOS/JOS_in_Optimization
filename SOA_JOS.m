function [Score,Position,SOAJOS]=SOA_JOS(N,T,maxRun,maxFE,BFid,nD,fhd,Jr)

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

disp('SOA-JOS is now tackling your problem')

SOAJOS=zeros(maxRun+1,T);

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
    
    
    Position=zeros(1,dim);
    Score=inf;
    
    Positions=init(N,dim,ub,lb);
    
    for i=1:N
        OP(i,:)=((ub-lb).*rand(size(lb)))+lb-Positions(i,:);
        RO(i,:) = rand*OP(i,:);
        DO(i,:) = Positions(i,:) + rand*(RO(i,:)-Positions(i,:));
        
        for j=1:dim
            if DO(i,j)<lb(1,j)
                DO(i,j)=lb(1,j);
            end
            if DO(i,j)>ub(1,j)
                DO(i,j)=ub(1,j);
            end
        end
    end
    Positions=DO;
    
    t=0;
    nFE=0;
    Row_Id=1;
    
    while nFE<maxFE
        
        for i=1:size(Positions,1)
            
            Flag4ub=Positions(i,:)>ub;
            Flag4lb=Positions(i,:)<lb;
            Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            
            fitness=feval(fhd,Positions(i,:)',BFid);
            nFE = nFE + 1;
            if nFE > maxFE; break; end
            if fitness<Score
                Score=fitness;
                Position=Positions(i,:);
                Row_Id=i;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %updating boundary for opposition after every iteration
            for x=1:size(Positions,1)
                for y=1:size(Positions,2)
                    if upper(1,y)<Positions(x,y)
                        upper(1,y)=Positions(x,y);
                    end
                    if lower(1,y)>Positions(x,y)
                        lower(1,y)=Positions(x,y);
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        t=t+1;
        SOAJOS(1,t) = nFE; %nFE
        SOAJOS(run+1,t) = Score;
        
        Fc=2-t*((2)/T);
        
        % Oppose the least fitness elements
        threshold=Fc;
        
        Positions=corOppose2(Positions,ub,lb,upper,lower,dim,threshold,Row_Id);
        
        for i=1:size(Positions,1)
            for j=1:size(Positions,2)
                
                r1=rand();
                r2=rand();
                
                A1=2*Fc*r1-Fc;
                C1=2*r2;
                b=1;
                ll=(Fc-1)*rand()+1;
                
                D_alphs=Fc*Positions(i,j)+A1*((Position(j)-Positions(i,j)));
                X1=D_alphs*exp(b.*ll).*cos(ll.*2*pi)+Position(j);
                Positions(i,j)=X1;
                
            end
        end
        if rand < Jr && nFE+size(Positions,1) < maxFE
            
            for i=1:N
                OP(i,:)=((ub-lb).*rand(size(lb)))+lb-Positions(i,:);
                RO(i,:) = rand*OP(i,:);
                DO(i,:) = Positions(i,:) + rand*(RO(i,:)-Positions(i,:));
                for j=1:dim
                    if DO(i,j)<lb(1,j)
                        DO(i,j)=lb(1,j);
                    end
                    if DO(i,j)>ub(1,j)
                        DO(i,j)=ub(1,j);
                    end
                end
            end
            Positions=DO;
        end
    end
    % display(['The best solution obtained by SOA is : ', num2str(Position)]);
    % display(['The best optimal value of the objective funciton found by SOA is : ', num2str(Score)]);
    toc
end
end
function Pos=init(SearchAgents,dim,ub,lb)

Boundary= size(ub,2);
if Boundary==1
    Pos=rand(SearchAgents,dim).*(ub-lb)+lb;
end

if Boundary>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Pos(:,i)=rand(SearchAgents,1).*(ub_i-lb_i)+lb_i;
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



