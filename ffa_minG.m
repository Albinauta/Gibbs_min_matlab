function y = ffa_minG(database,P,varargin)%%%funzione euristica, usa il firefly minimization algorithm per trovare una possibile composizione, da usare se findcomp da problemi
if numel(varargin) < 2
    if numel(varargin{1}) > 1
        varargin = varargin{1};
    else
        fprintf("ERROR: There must be at least 2 species as input in order to find equilibrium's composition")
        return
    end
end
%%%%%PARAMETRI%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%COSTRUZIONE VETTORI REAGENTI%%%%%
R = 0.008314;  %[KJ/molK]
n0 = [];%%MOLI REAGENTI
Re = {};
numeratore = 0;%numeratore della media pesata per il calcolo della temperatura di equilibrio
denominatore = 0; % denominatore media pesata
active_species = {};
inert = {};
inert_prod= {};
nInert = [];
n_prodotti = 0;
Matrix_const = zeros(numel(varargin),8);%%matrice constanti termodinamiche di tutte le specie
DHR = zeros(numel(varargin,1));
for i = 1 : numel(varargin)
    if varargin{i}.MassFlux > 0
        n0(end+1) = varargin{i}.MassFlux;%%%include anche gli inerti iniziali!!!
        Re{end+1} = varargin{i};%%%include anche gli inerti iniziali!
        cp = varargin{i}.getcp(varargin{i}.Temp);
        numeratore = numeratore + (varargin{i}.MassFlux * cp * varargin{i}.Temp);
        denominatore = denominatore + (varargin{i}.MassFlux * cp);
        if varargin{i}.Inert == true
            inert{end+1} = varargin{i};
            inert_prod{end+1} = copy(varargin{i}.copy);
            nInert(end+1)= varargin{i}.MassFlux;

        end
    else
        n_prodotti = n_prodotti + 1;
    end
    if varargin{i}.Inert == false
        active_species{end+1} = copy(varargin{i}.copy);
    end
    Matrix_const(i,:) = varargin{i}.getConst;
    DHR(i) = str2num(varargin{i}.Properties(12));
end
Matrix_const_act = zeros(numel(active_species),8);%%% solo specie attive
COMPeq = [active_species, inert_prod];%%array oggetti equilibrio
Teq = numeratore / denominatore;
DHin = [];
Aeq = zeros(7,numel(active_species));
beq = zeros(7,1);%vettore moli iniziali; la funzione minimizza con il vincolo Aeq*ns = beq.
for i = 1 : numel(active_species)
    Aeq(:,i) = active_species{i}.Comp;%matrice delle composizioni di ogni specie che compare nell'equilibrio(eccetto inerti)
    comp = varargin{i}.Comp;
    Matrix_const_act(i,:) = active_species{i}.getConst;
    beq(1) = beq(1) + (active_species{i}.MassFlux *comp(1));
    beq(2) = beq(2) + (active_species{i}.MassFlux *comp(2));
    beq(3) = beq(3) + (active_species{i}.MassFlux *comp(3));
    beq(4) = beq(4) + (active_species{i}.MassFlux *comp(4));
    beq(5) = beq(5) + (active_species{i}.MassFlux *comp(5));
    beq(6) = beq(6) + (active_species{i}.MassFlux *comp(6));
    beq(7) = beq(7) + (active_species{i}.MassFlux *comp(7));
end
for i = 1 : numel(Re)
    DHin(end+1) = Re{i}.intgH(database,Teq);%%%tiene conto anche degli inerti iniziali!!!!!
end
%%%%%%DEFINIZIONE DELLE FUNZIONI PROBLEMA E VINCOLI%%%%%%%%
    function G = cost(nj) %%%%FUNZIONE PROBLEMA: E' LA FUNZIONE CHE CERCO DI MINIMIZZARE%%%%
        T = nj(end);
        T = T/1000;
        nmass = nj(1,1:end-1);
        T_vecH = [T;T^2 /2; T^3 /3; T^4 /4; -1/T;1;0;-1];
        t_vecS = [log(T);T;T^2 /2;T^3 / 3; 1/(2*T^2);0;1;0];%%%non ho usato gli oggetti perchè voglio un cp indipendente dalla T(programma + veloce)
        DH = Matrix_const_act * T_vecH;
        DS = Matrix_const_act * t_vecS;
        G0 = DH - T * DS;
        Enj = sum(nj(1,1:end-1));
        G = sum(nmass.*(G0/R/T + log(nmass/Enj*P)));
        G = sum(G);
    end
    function V = vincoli(nj)%%chiamata per calcolare la "penalità" di ogni composizione generata
        lambda = 10^5;
        v1 = sum(abs(Aeq*nj(1,1:end-1)' -beq));
        T = nj(end)/1000;
        T_vecH = [T;T^2 /2; T^3 /3; T^4 /4; -1/T;1;0;-1];
        DHf = DHR' + Matrix_const *T_vecH ;
        ntot = [nj(1,1:end-1),nInert];
        v2 = sum(abs(n0*(DHin') - ntot*(DHf)));
        V = lambda*(v1+v2);%%da rivedere
    end
nRe = [];
M = 0;
for i = 1:numel(Re)
    if Re{i}.Inert == false
        m = Re{i}.MassFlux;
        nRe(end+1) = m;
        if m > M
            M = m;
        end
    end
end
%%%%COMPOSIZIONI LIMITI SUPERIORI ED INFERIORI%%%%%
nre = sum(n0) - sum(nInert);
neq =nre.*ones(1, n_prodotti);
neq = 5.*ones(1,n_prodotti);
nEQ = [nRe,neq];
LB = [zeros(1,numel(active_species)),298];
UB = [nEQ,5*Teq];
d = numel(active_species) +1;
u0 = findcomp(database,P,COMPeq);%%%prima stima della comp
%u0 = [findcomp(database,P,active_species),Teq];
%u0 = LB + (UB-LB).*rand(1,d);
z = FFA2(@cost,@vincoli,UB,LB,d,u0);
plot(z{2},'b--o')
y = z{1};
end