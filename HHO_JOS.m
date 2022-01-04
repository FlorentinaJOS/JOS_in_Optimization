function [Rabbit_Energy,Rabbit_Location,HHOJOS]=HHO_JOS(N,T,maxRun,maxFE,BFid,nD,fhd,Jr)

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

disp('HHO-JOS is now tackling your problem')

HHOJOS=zeros(maxRun+1,T);

for run=1:maxRun
    tic
    % initialize the location and Energy of the rabbit
    Rabbit_Location=zeros(1,dim);
    Rabbit_Energy=inf;
    
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
    
    %Initialize the locations of Harris' hawks
    X=initialization(N,dim,ub,lb);
    %X=initialization(N,dim,up,down);
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
    
    t=0; % Loop counter
    nFE=0;
    h=size(X,1);
    fitness = zeros(1,h); %row

    while nFE<maxFE
 
        for i=1:h
            % Check boundries
            FU=X(i,:)>ub;
            FL=X(i,:)<lb;
            X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;

            fitness(i)=feval(fhd,X(i,:)',BFid);
            nFE = nFE + 1;
            if nFE > maxFE; break; end        
            
            % Update the location of Rabbit
            if fitness(i)<Rabbit_Energy
                Rabbit_Energy=fitness(i);
                Rabbit_Location=X(i,:);
                Rabbit_Row_Id = i;  %%
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
        
        %nFE=nFE+h;
        t=t+1;
        HHOJOS(1,t) = nFE; %nFE
        HHOJOS(run+1,t) = Rabbit_Energy;
        
        E1=2*(1-(t/T)); % factor to show the decreaing energy of rabbit
        
        % Oppose the least fitness elements
        threshold=E1;
      
        X=corOppose2(X,fitness,ub,lb,upper,lower,dim,threshold,Rabbit_Row_Id);
        
        % Update the location of Harris' hawks
        for i=1:h
            E0=2*rand()-1; %-1<E0<1
            Escaping_Energy=E1*(E0);  % escaping energy of rabbit
            
            if abs(Escaping_Energy)>=1 %15.4 percent
                %% Exploration:
                % Harris' hawks perch randomly based on 2 strategy:
                
                q=rand();
                rand_Hawk_index = floor(N*rand()+1);
                X_rand = X(rand_Hawk_index, :);
                if q<0.5
                    % perch based on other family members
                    X(i,:)=X_rand-rand()*abs(X_rand-2*rand()*X(i,:)); %TLBO
                elseif q>=0.5
                    % perch on a random tall tree (random site inside group's home range)
                    X(i,:)=(Rabbit_Location(1,:)-mean(X))-rand()*((ub-lb)*rand+lb); %QPSO
                end
                
            elseif abs(Escaping_Energy)<1
                %% Exploitation: 84.6 percent
                % Attacking the rabbit using 4 strategies regarding the behavior of the rabbit
                
                %% phase 1: surprise pounce (seven kills)
                % surprise pounce (seven kills): multiple, short rapid dives by different hawks
                
                r=rand(); % probablity of each event
                
                if r>=0.5 && abs(Escaping_Energy)<0.5 % Hard besiege  %%21.149 percent
                    X(i,:)=(Rabbit_Location)-Escaping_Energy*abs(Rabbit_Location-X(i,:));
                end
                
                if r>=0.5 && abs(Escaping_Energy)>=0.5  % Soft besiege %% 21.149
                    Jump_strength=2*(1-rand()); % random jump strength of the rabbit
                    X(i,:)=(Rabbit_Location-X(i,:))-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
                end
                
                %% phase 2: performing team rapid dives (leapfrog movements)
                if r<0.5 && abs(Escaping_Energy)>=0.5, % Soft besiege % rabbit try to escape by many zigzag deceptive motions
                    
                    Jump_strength=2*(1-rand()); % 0< Jump_strength <2
                    X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:));
                    
                    if feval(fhd,X1',BFid) < feval(fhd,X(i,:)',BFid)
                        X(i,:)=X1;
                        
                    else % hawks perform levy-based short rapid dives around the rabbit
                        X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-X(i,:))+rand(1,dim).*Levy(dim);
                        if feval(fhd,X2',BFid) < feval(fhd,X(i,:)',BFid)
                            X(i,:)=X2;
                        end
                    end
                end
                
                if r<0.5 && abs(Escaping_Energy)<0.5, % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                    % hawks try to decrease their average location with the rabbit
                    Jump_strength=2*(1-rand());
                    X1=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X));
                    
                    if feval(fhd,X1',BFid) < feval(fhd,X(i,:)',BFid)
                        X(i,:)=X1;
                    else % Perform levy-based short rapid dives around the rabbit
                        X2=Rabbit_Location-Escaping_Energy*abs(Jump_strength*Rabbit_Location-mean(X))+rand(1,dim).*Levy(dim);
                        if feval(fhd,X2',BFid) < feval(fhd,X(i,:)',BFid)
                            X(i,:)=X2;
                        end
                    end
                end
                %%
            end
        end
        
        
        if rand < Jr && nFE+h < maxFE
            
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
end

% ___________________________________
function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);
step=u./abs(v).^(1/beta);
o=step;
end

function [X]=initialization(N,dim,up,down)

if size(up,1)==1
    X=rand(N,dim).*repmat((up-down)+down,N,1);
end
if size(up,1)>1
    for i=1:dim
        high=up(i);low=down(i);
        X(:,i)=rand(1,N).*(high-low)+low;
    end
end
end

function [Positions]=corOppose2(Positions,fitness,ub,lb,upper,lower,dim,threshold,Rabbit_Row_Id)

[n,b] = size(Positions);
for i=1:n(1) 
    
    if i ~= Rabbit_Row_Id
        
        sum=0;
        greater=[];
        less=[];
        x=1;z=1;y=1;
        
        for j=1:b
            d(x)=abs(Positions(Rabbit_Row_Id,j)-Positions(i,j));
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



