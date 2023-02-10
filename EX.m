function y = EX(database,set1,compEQ,varargin)%%%funzione di calcolo exergia, accetta un cell array di oggetti e un array di numeri ( ottenuti dalle funzioni adiabatiche)
if numel(varargin) == 0
    varargin = zeros(1,numel(set1));
else
    varargin = varargin{1};
end
R = 0.008314;%Kj/molK
numeratore = 0;
denominatore = 0;
%%%calcolo exergia
T = compEQ(end);
exfisR = [];
exfisP = [];
exchimR = [];
exchimP = [];
nreg = 0;
%DStot = 0;
%%%%calcolo massa totale reagenti e temperatura media%%%%
for i = 1 :numel(set1)
    if set1{i}.MassFlux~= 0
        %DH = set1{i}.intgH(database,T);
        %DS = set1{i}.intgS(database,T);
        %DStot = DStot + DS;
        %exfisR(end+1) = DH + 298.2 * DS;
        nreg = nreg + set1{i}.MassFlux;
        cp = set1{i}.getcp(set1{i}.Temp);
        numeratore = numeratore + (set1{i}.MassFlux * cp * set1{i}.Temp);
        denominatore = denominatore + (set1{i}.MassFlux * cp);
    end
end
Tr = numeratore/denominatore;
%%%%%calcolo exergia chimica e fisica dei reagenti e prodotti%%%%%
for i = 1 :numel(set1)
    if set1{i}.MassFlux~= 0
        DHr = set1{i}.intgH(database,Tr);%Kj
        DSr = set1{i}.intgS(database,Tr);%j
        %mass = set1{i}.MassFlux;
        %c = DHr-str2num(set1{i}.Properties(12));
        %d = (DSr-str2num(set1{i}.Properties(13)))*10^-3;
        xk = (set1{i}.MassFlux) / nreg;
        exfisR(end+1) = xk*(DHr-str2num(set1{i}.Properties(12))+varargin(i) + 298.2 * (DSr-str2num(set1{i}.Properties(13)))*10^-3);%tiene conto di eventuali contributi ulteriori di entalpia
        exchimR(end+1) = xk*str2num(set1{i}.Properties(14)) + R*298.2*xk*log(xk);
        exfisP(end+1) = 0;
        exchimP(end+1) = 0;
    else
        DH = set1{i}.intgH(database,T);
        DS = set1{i}.intgS(database,T);
        xkP = (compEQ(i))/sum(compEQ(1:end-1));
        exfisP(end+1) = xkP*(DH-str2num(set1{i}.Properties(12)) + 298.2*(DS-str2num(set1{i}.Properties(12)))*10^-3);
        exchimP(end+1) = xkP*str2num(set1{i}.Properties(14)) + R*298.2*xkP*log(xkP);
    end   
end
ExfisR = sum(exfisR);
ExfisP = sum(exfisP);
ExchimR = sum(exchimR);
ExchimP = sum(exchimP);%% valori di debugging
ExR = sum(exfisR + exchimR);
ExP = sum(exfisP + exchimP);
%%%%calcolo exergia distrutta%%%%
%Sgen = DStot + loss/T;
%ExD = 298.2*Sgen;
%s=sum(Ex);
y = ExP/ExR;

end
