function y = findcomp(database,P,varargin)%%funzione per la minimizzazione di gibbs in condizioni adiabatiche, viene applicato il vincolo di entalpia
if numel(varargin) < 2
    if numel(varargin{1}) > 1
        varargin = varargin{1};% la funzione accetta oggetti singoli oppure un cell array di oggetti
    else
        fprintf("ERROR: There must be at least 2 species as input in order to find equilibrium's composition")
        return
    end
end
%%%%%%%%%COSTRUZIONE VETTORI REAGENTI%%%%%
n0 = [];%%MOLI REAGENTI
Re = {};
numeratore = 0;%numeratore della media pesata per il calcolo della temperatura di equilibrio
denominatore = 0; % denominatore media pesata
active_species = {};
inert = {};
inert_prod= {};
nInert = [];
for i = 1: numel(varargin)
    if varargin{i}.MassFlux > 0
        n0(end+1) = varargin{i}.MassFlux;%%include anche le moli degli inerti iniziali!!
        Re{end+1} = varargin{i};%%%include gli inerti!!!
        cp = varargin{i}.getcp(varargin{i}.Temp);
        numeratore = numeratore + (varargin{i}.MassFlux * cp * varargin{i}.Temp);
        denominatore = denominatore + (varargin{i}.MassFlux * cp);
        if varargin{i}.Inert == true
            inert{end+1} = varargin{i};
            inert_prod{end+1} = copy(varargin{i}.copy);
            nInert(end+1)= varargin{i}.MassFlux;

        end
    end
    %%%%%%%%%%FILTRO SPECIE INERTI%%%%%%%%%%%%%%%%%%%
    if varargin{i}.Inert == false
        active_species{end+1} = copy(varargin{i}.copy);
    end
end
COMPeq = [active_species, inert_prod];%%array oggetti equilibrio
%%%%CALCOLO ENTALPIA REAGENTI%%%%%%
Teq = numeratore / denominatore;
DHin = [];
for i = 1 : numel(Re)
    DHin(end+1) = Re{i}.intgH(database,Teq);
end
%%%%%VALORI DI INIZIALIZZAZIONE DEL CICLO WHILE%%%%%%%
T_0 = Teq;
DT = 10000000;%%tolleranza errore
Niter = 0;
tr = 298.2/1000;
%TR = [tr;tr^2 /2; tr^3 / 3; tr^4 / 4; -1/tr];
%%%% VINCOLO DI ADIABATICITA'%%%%%
function DH = findT(T)
        T = T/1000;
        T_vec = [T;T^2 /2; T^3 /3; T^4 /4; -1/T;1;0;-1];
        %test = Matrix_const * TR;
        j = DHR' + Matrix_const * T_vec;
        DH = n0*(DHin')-ntot*j;
end
vettore_temp = [];
vettore_iter = [];
%%%%%%CICLO WHILE%%%%%
while DT > 0,001;
        if Niter == 20
            fprintf('MAX NUMBER OF ITERATIONS REACHED')
            vettore_temp(end+1) = T_0;
            vettore_iter(end+1) = Niter;
            plot(vettore_iter,vettore_temp,'b--o')
            y = [Neq,T_0];
            for i=1:numel(Neq)
                fprintf('%5d%10s%10.3g\n',i,COMPeq{i}.Name,Neq)
            end
        return
        end
        if Niter >= 10%%%nel caso di oscillazioni della soluzione, viene restituita la composizione alla T media tra i limiti di oscillazione
            if T_0 - vettore_temp(end) > 20
                if T_0 - vettore_temp(end-1) < 10 && vettore_temp(end) - vettore_temp...
                        (end-2) < 10
                    T_0 = (T_0+vettore_temp(end))/2;
                     for i = 1: numel(active_species)
                         prop = database.getVec(active_species{i}.Name,T_0);
                         active_species{i}.Properties = prop;
                         active_species{i}.Temp = T_0;
                         active_species{i}.calcDG(database);
                     end
                     Neq = minG2(database,P,active_species);
                     vettore_iter(end+1) = Niter;
                     vettore_temp(end+1) = T_0;
                     y = [Neq,T_0];
                     for i=1:numel(Neq)
                         fprintf('%5d%10s%10.3g\n',i,COMPeq{i}.Name,Neq(i))
                     end
                     plot(vettore_iter,vettore_temp,'b--o')
                     return

                end
            end
        end
    %%%AGGIORNAMENTO PROPRIETA' MISCELA%%%%
    vettore_temp(end+1) = T_0;
    for i = 1: numel(active_species)
        prop = database.getVec(active_species{i}.Name,T_0);
        active_species{i}.Properties = prop;
        active_species{i}.Temp = T_0;
        active_species{i}.calcDG(database);
    end
    if numel(inert_prod) > 0
        for i = 1 : numel(inert_prod)
            prop = database.getVec(inert_prod{i}.Name,T_0);
            inert_prod{i}.Properties = prop;
            inert_prod{i}.Temp = T_0;
            inert_prod{i}.calcDG(database);
        end
    end
    Neq = minG2(database,P,active_species); % composizione all'equilibrio a T_0
    ntot = [Neq,nInert];
    Matrix_const = zeros(numel(COMPeq),8);%%matrice constanti termodinamiche
    DHR = zeros(numel(COMPeq,1));%vettore entalpie standard di formazione
    for i = 1:numel(COMPeq)
        Matrix_const(i,:) = COMPeq{i}.getConst;
        DHR(i) = str2num(COMPeq{i}.Properties(12));
    end
    x = lsqnonlin(@findT,T_0,298.2);%%%temperatura che risolve il vincolo di adiabacità per la data composizione all'equilibrio
    %x = fsolve(@findT,T_0);
    DT = abs(T_0 - x);
    T_0 = x;
    Niter = Niter + 1;
    vettore_iter(end+1) = Niter;
end

y = [Neq,T_0];
for i=1:numel(Neq)
    fprintf('%5d%10s%10.3g\n',i,COMPeq{i}.Name,Neq(i))
end
plot(vettore_iter,vettore_temp,'b--o')
end


