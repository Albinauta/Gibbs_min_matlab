function X = FFA2(funmin,limiti,UB,LB,d,u0)%%%algoritmo generico firefly che usa gli oggetti SpecieChimica
rng("shuffle")
%%%%INIZIALIZZAZIONE VALORI DI CONTROLLO DELLA GENERAZIONE POP E DEL LORO
%%%%MOVIMENTO
alfa = 0.2;
%alfa = 0.01*sum(UB-LB);
%beta0 = 0.8;
beta0 = 1;
gamma = 1;
%gamma = 0.01;
delta = 0.97;
N = 100;
POPS = {};
risultati = [];
bestfit = {};
solution = [];
MAXGENERATIONS = 100;
%%%%%%GENERAZIONE POPOLAZIONE%%%%%%%%%%%%%%%
for i = 1:N/2
    pop = LB +(u0-LB).*rand(1,d);
    POPS{end +1} = pop;
end
for i = 1:N/2%%%i pop generati sono incentrati su u0
    pop = u0 +(UB-u0).*rand(1,d);
    POPS{end +1} = pop;
end
%%%%%%%CONFRONTO DELLE LUCCIOLE%%%%%%%%%%
for k = 1:MAXGENERATIONS
    risultati = [];
   for i = 1: N
        L1 = POPS{i};
        F1 = funmin(L1) + limiti(L1);%valore attrazione
        risultati(end+1) = F1;
   end
   [value,index] = sort(risultati);
   bestfirefly = POPS{index(1)};
   bestfirefly_value = value(1);
   for j = 1:N
       if j ~= index(1)
           L2 = POPS{j};
           %F2 = risultati(j);
           D = sqrt(sum((L2-bestfirefly).^2));
           %I1 = bestfirefly_value./D;% luminosità
           %I2 = F2./D;
           L2_new = L2 + beta0*exp(-gamma*D^2)*(L2-bestfirefly) + alfa.*(rand -1/2)*(UB-LB);
           %%%controllo vincoli di bordo%%%%%
           matchUB = find(L2_new > UB);
           matchLB = find(L2_new < LB);
           if numel(matchUB) > 0
               for h = 1:numel(matchUB)
                   L2_new(matchUB(h)) = UB(h);
               end
           elseif numel(matchLB) > 0
               for h =1:numel(matchLB)
                   L2_new(matchLB(h)) = LB(h);
               end
           end
           POPS{j} = L2_new;
       end
   end
   if numel(bestfit) == 0
       bestfit{1} = bestfirefly;
       solution(end+1) = bestfirefly_value;
       %iter(end+1) = numel(solution);
   else
       Fbest = solution(end);
       if bestfirefly_value < Fbest
           bestfit{1} = bestfirefly;
           solution(end+1) = bestfirefly_value;
           %iter(end+1) = numel(solution);
       end
   end
alfa = alfa*delta; %riduce la randomness ad ogni iterazione
end
X = {bestfit{1},solution};
end
