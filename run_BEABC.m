clear all
close all
clc

global sequence
sequence = [0 1 0 1 0 0 1 0 1 0 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 0 1 0 1 0 0 1 0 0 1 0 1 0 0 1 0 1 0 0 1 0 0 1 0 1 0 0 1 0]; % D55
D = size(sequence,2)*2 - 5;
maxCycle = 5000; 
runtime = 5;

alpha = 0.9;
pops = 0;
objfun = 'libai';

NP = 40;
FoodNumber=NP/2; 
GlobalMins=zeros(1,runtime);
ObjVal = zeros(FoodNumber,1);
Fitness = ObjVal;
trial=zeros(1,FoodNumber);
ObjValSol=ObjVal;
FitnessSol = ObjVal;
beabc=zeros(runtime,maxCycle);
qq=zeros(1,maxCycle);


for r=1:runtime

ub = ones(1,D)*180;
lb = ones(1,D).*-180;
Range = repmat((ub-lb),[FoodNumber 1]);
Lower = repmat(lb, [FoodNumber 1]);
Foods = rand(FoodNumber,D) .* Range + Lower;


for ky = 1:FoodNumber
    ObjVal(ky,1) = libai(Foods(ky,:));
    Fitness(ky,1) = calculateFitness(ObjVal(ky,1));
end

trial = ones(1,FoodNumber);

BestInd=find(ObjVal==min(ObjVal));
BestInd=BestInd(end);
GlobalMin=ObjVal(BestInd);
GlobalParams=Foods(BestInd,:);

iter=1;
while ((iter <= maxCycle)),
    sequenca = randperm(D);
    for i=1:(FoodNumber)

                                sequenca = randperm(D);
                                list = sequenca(1,1:trial(i));
                                sol = Foods(i,:);
                                for j = 1:trial(i)
                                    kiko = list(1,j);
                                    neighbour1 = fix(rand*(FoodNumber)) + 1;
                                    neighbour2 = fix(rand*(FoodNumber)) + 1;
                                    while(neighbour2 == i)
                                        neighbour2 = fix(rand*(FoodNumber)) + 1;
                                    end;
                                    vv = trial(i)./(trial(i) + trial(neighbour2));
                                    sol(kiko) = Foods(neighbour1,kiko) + (Foods(i,kiko) - Foods(neighbour2,kiko)).*(2.*rand-1).*vv;
                                end


                               %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
                                ind=find(sol<lb);
                                sol(ind)=lb(ind);
                                ind=find(sol>ub);
                                sol(ind)=ub(ind);

                                %evaluate new solution
                                ObjValSol=feval(objfun,sol);
                                FitnessSol=calculateFitness(ObjValSol);

                               % /*a greedy selection is applied between the current solution i and its mutant*/
                               if (FitnessSol>Fitness(i)) %/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
                                    Foods(i,:)=sol;
                                    Fitness(i)=FitnessSol;
                                    ObjVal(i)=ObjValSol;
                                    trial(i) = 1;
                               end;   
    end;


    prob = Fitness./max(Fitness);
    i=1;
    t=0;
    for t = 1:FoodNumber
        if(rand < prob(i))
                                        %/*The parameter to be changed is determined randomly*/
                                        Param2Change = fix(rand*D)+1;
                                        neighbour = fix(rand*(FoodNumber))+1;      
                                            while(neighbour==i)
                                                neighbour=fix(rand*(FoodNumber))+1;
                                            end;
                                       sol=Foods(i,:);
                                       vv = trial(i)./(trial(i) + trial(neighbour));
                                       sol(Param2Change)=Foods(i,Param2Change)+(Foods(i,Param2Change)-Foods(neighbour,Param2Change))*(rand-0.5)*2.*vv;

                                       %  /*if generated parameter value is out of boundaries, it is shifted onto the boundaries*/
                                        ind=find(sol<lb);
                                        sol(ind)=lb(ind);
                                        ind=find(sol>ub);
                                        sol(ind)=ub(ind);

                                        %evaluate new solution
                                        ObjValSol=feval(objfun,sol);
                                        FitnessSol=calculateFitness(ObjValSol);

                                       % /*a greedy selection is applied between the current solution i and its mutant*/
                                       if (FitnessSol>Fitness(i)) %/*If the mutant solution is better than the current solution i, replace the solution with the mutant and reset the trial counter of solution i*/
                                            Foods(i,:)=sol;
                                            Fitness(i)=FitnessSol;
                                            ObjVal(i)=ObjValSol;
                                            trial(i) = 1;
                                       else
                                            trial(i) = trial(i) + 1; %/*if the solution i can not be improved, increase its trial counter*/
                                       end;
        end;
        i=i+1;
        if (i==(FoodNumber)+1) 
            i=1;
        end;  
    end;
    ind=find(ObjVal == min(ObjVal));
    ind=ind(end);
    if (ObjVal(ind)<GlobalMin)
        GlobalMin=ObjVal(ind);
        GlobalParams=Foods(ind,:);
    end;
    
    ind = find(trial > D);
    trial(ind) = D;
    if (mean(trial) > (D*alpha))
                pops = pops + 1;
                sequence = randperm(NP);
                num = fix(NP*alpha) + 1;
                list3 = sequence(1:num);
                for ii = 1:num 
                    sol = (ub-lb).*rand(1,D) + lb;
                    ObjValSol = feval(objfun,sol);
                    FitnessSol = calculateFitness(ObjValSol);
                    Foods(list3(ii),:) = sol;
                    Fitness(list3(ii)) = FitnessSol;
                    ObjVal(list3(ii)) = ObjValSol;
                end
    end;            

            
            
    if (mod(iter-1,100)==0)
        fprintf('iter=%d ObjVal=%f\n',iter,GlobalMin);
    end
    be(r,iter) = GlobalMin;
    iter=iter+1;

    end % End of ABC
    GlobalMins(r)=GlobalMin;
end

save ('be_d55.mat', 'be')
plot(mean(be))

