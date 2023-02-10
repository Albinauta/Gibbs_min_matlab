function y = HEAT(database,Spchim,Tin)
Tfin = Spchim.Temp;
DHin = Spchim.intgH(database,Tin);
DHtot = Spchim.intgH(database,Tfin);
y = DHtot-DHin;% il segno segue la convenzione del calore: positivo se deve essere fornito, negativo se deve essere sottratto(ossia viene emesso/ devo raffreddare)
end