clc
clear
close all
format shortG

%% Parameters Setting

nvar=5;                % Number Of Variable
lb=-10*ones(1,nvar);    % Lower Bound
ub= 10*ones(1,nvar);    % Upper Bound

maxiter=1000;     % Maximum Of Iteration

npop=20;          % Harmony Memory Size

HMCR=0.8;         % Harmony Memory Consideration Rate

PAR=0.9;          % Piterch Adjustment Rate

BW=1;           % Bandwidth

BW_RF=0.95;      % BandWidth Reduction Factor 

%% Initial Population
tic
emp.x=[];
emp.fit=[];

pop=repmat(emp,npop,1);

for i=1:npop
    pop(i).x=unifrnd(lb,ub);
    pop(i).fit=fitness(pop(i).x);
end

% Sorting
[~, ind]=sort([pop.fit]);
pop=pop(ind);

gpop=pop(1);  % Global Solutuion



%%  Main Loop

BEST=zeros(maxiter,1);
MEAN=zeros(maxiter,1);

for iter=1:maxiter
    
    newpop=pop;

    for i=1:npop
        
        % Consideration
        J=rand(1,nvar)<HMCR;
        J=find(J==1);
                
        for j=J
             k=randi([1 npop]);
             newpop(i).x(j)=pop(k).x(j);
             
%              if rand<PAR
%                  d=BW*unifrnd(-1,1);
%                  newpop(i).x(j)=newpop(i).x(j)+d;
%              end    

        end
        
        
        % Piterch Adjustment
         J=rand(1,nvar)<PAR;
         J=find(J==1);     
         d=BW*unifrnd(-1,1,size(J));
         newpop(i).x(J)=newpop(i).x(J)+d; 
         
        % Check Bound 
        newpop(i).x=CB(newpop(i).x,lb,ub);          
         
        newpop(i).fit=fitness(newpop(i).x);
        
    end
    
    % Merge
    [pop]=[pop;newpop];
    
    % Sorting
    [~, ind]=sort([pop.fit]);
    pop=pop(ind);
    pop=pop(1:npop);
    
    gpop=pop(1); 

    BW=BW*BW_RF;

    
    BEST(iter)=gpop.fit;
    MEAN(iter)=mean([pop.fit]);

    disp(['Iter ' num2str(iter) ' BEST = ' num2str(BEST(iter))]);
    
    
end

%% Results

disp(' ')
disp([ ' BEST solution = '  num2str(gpop.x)]);
disp([ ' BEST fitness = '  num2str(gpop.fit)]);
disp([ ' Time = '  num2str(toc)]);

figure
semilogy(BEST,'r');
hold on
semilogy(MEAN,'b');
xlabel('Iteration');
ylabel('Fitness');
legend('BEST','MEAN');
