function x = minG2(database,P,varargin)% funzione di minimizzazione dell'energia libera di gibbs di una miscela di composti, in condizioni isotermiche
if numel(varargin) < 2
    if numel(varargin{1}) > 1%%La funzione accetta gruppi di singoli oggetti o un cell array che li contiene
        varargin = varargin{1};
    else
        fprintf("ERROR: There must be at least 2 species as input in order to find equilibrium's composition")
        return
    end
end
species = {};
for i = 1 : numel(varargin)
    %z = varargin{1,i};
    species{i} = varargin{i}.Name;
end
%NOMI = [species{1,:}];
R = 0.008314;  %[KJ/molK]
T = varargin{i}.Temp; %K
%%%ESTRAZIONE GIBS FORMAZIONE A T%%%%%%%%%%%
GJ0 = {};
for i = 1:numel(varargin)
    if varargin{i}.DG0 ~= 0
        if varargin{i}.DG ~= 0
            GJ0{i} = varargin{i}.DG;
        else
            varargin{i}.calcDG(database)
            GJ0{i} = varargin{i}.DG;
        end
    else
        varargin{i}.calcDG0(database)
        if varargin{i}.DG ~= 0
            GJ0{i} = varargin{i}.DG;
        else
            varargin{i}.calcDG(database)
            GJ0{i} = varargin{i}.DG;
        end
    end
end
GJ0 = [GJ0{1,:}];
%%%DEFINIZIONE FUNZIONE DI GIBBS DA MINIMIZZARE%%%%
    function G = func(nj)
        Enj = sum(nj);
        G = sum(nj.*(GJ0/R/T + log(nj/Enj*P)));
    end
Aeq = zeros(7,numel(varargin));
for i = 1 : numel(varargin)
    Aeq(:,i) = varargin{i}.Comp;%matrice delle composizioni di ogni specie che compare nell'equilibrio
end
beq = zeros(7,1);%vettore moli iniziali; la funzione minimizza con il vincolo Aeq*ns = beq.
Mtot = 0; %calcolo le moli totali entranti per definire una distribuzione iniziale dell'equilibrio
for i = 1:numel(varargin)
    comp = varargin{i}.Comp;
    Mtot = Mtot + varargin{i}.MassFlux;
    beq(1) = beq(1) + (varargin{i}.MassFlux *comp(1));
    beq(2) = beq(2) + (varargin{i}.MassFlux *comp(2));
    beq(3) = beq(3) + (varargin{i}.MassFlux *comp(3));
    beq(4) = beq(4) + (varargin{i}.MassFlux *comp(4));
    beq(5) = beq(5) + (varargin{i}.MassFlux *comp(5));
    beq(6) = beq(6) + (varargin{i}.MassFlux *comp(6));
    beq(7) = beq(7) + (varargin{i}.MassFlux *comp(7));
end
LB = zeros(1,numel(varargin)); % limite inferiori per le moli
x0 = zeros(1,numel(varargin)); % inzializzo il vettore della config iniziale
rng(14,'twister');
for i = 1 : numel(varargin)
    x0(i) = Mtot/numel(varargin);%la config iniziale assume una distribuzione media delle moli all'equilibrio
end
%%%%MINIMIZZAZIONE%%%%%
%options = optimset('Algorithm','sqp');
opts = optimoptions(@fmincon,'Algorithm','sqp');
problem = createOptimProblem('fmincon','x0',x0,'objective',@func,'Aeq',Aeq,'beq',beq,'lb',LB,'options',opts);%%E' RICHIESTA LA Global Optimization Toolbox PER GIRARE IL PROGRAMMA
gs = GlobalSearch('Display','iter');
x = run(gs,problem);
%x = fmincon(@func,x0,[],[],Aeq,beq,LB,[],[],options);
%for i=1:numel(x)
%    fprintf('%5d%10s%10.3g\n',i,species{i},x(i))
%end
end