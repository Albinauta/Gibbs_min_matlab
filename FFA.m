function X = FFA(funmin,limiti,UB,LB,d,u0)%%%algoritmo generico firefly che usa gli oggetti SpecieChimica
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
maxiter = 20;
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
    for i = 1: N
        L1 = POPS{i};
        F1 = funmin(L1) + limiti(L1);%valore attrazione
        risultati(end+1) = F1;
    end
    for i = 1:N
        L1 = POPS{i};
        F1 = risultati(i);
        for j = 1:maxiter
            %j = randi(N,1);
            if i ~= j
                L2 = POPS {j};
                F2 = risultati(j);
                D = sqrt(sum((L2-L1).^2));
                I1 = F1./D;% luminosità
                I2 = F2./D;
                %L1_new = zeros(1,d);
                if I1 > I2
                    %L1_new = L1 + beta0*exp(-gamma*D^2)*(L1-L2) + alfa*(rand +1/2)*(UB-LB);%%%movimento della lucciola
                    L1_new = L1 + beta0*exp(-gamma*D^2)*(L1-L2) + alfa.*(rand -1/2)*(UB-LB);
                    %%%controllo vincoli di bordo%%%%%
                    matchUB = find(L1_new > UB);
                    matchLB = find(L1_new < LB);
                    if numel(matchUB) > 0
                        for h = 1:numel(matchUB)
                            L1_new(matchUB(h)) = UB(h);
                        end
                    elseif numel(matchLB) > 0
                        for h =1:numel(matchLB)
                            L1_new(matchLB(h)) = LB(h);
                        end
                    end
                    F1_new = funmin(L1_new) + limiti(L1_new);
                    D_new = sqrt(sum((L2-L1_new).^2));
                    I1_new = F1_new./D_new;
                    if I1_new < I1%confronto attrazione tra nuova posizione e vecchia
                        L1 = L1_new;
                        F1 = F1_new;
                        POPS{i} = L1;
                        risultati(i) = F1;
                    end
                end
            end
        end
        %iter = [];
        %%%%AGGIORNAMENTO DEL BESTFIT%%%%%%%%%%
        if numel(bestfit) == 0
            bestfit{1} = L1;
            solution(end+1) = F1;
            %iter(end+1) = numel(solution);
        else
            Fbest = solution(end);
            if F1 < Fbest
                bestfit{1} = L1;
                solution(end+1) = F1;
                %iter(end+1) = numel(solution);
            end
        end
    end
    alfa = alfa*delta; %riduce la randomness ad ogni iterazione
end
X = {bestfit{1},solution};
end
