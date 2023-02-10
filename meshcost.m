function y = meshcost(flowrate,CO2perc)
H2Operc = 12.54;
O2perc = 6.98;
Arperc = 0.7;
MaxPerc = 100-(H2Operc+O2perc+Arperc);
x = linspace(CO2perc(1),CO2perc(2),1000);
q = linspace(flowrate(1),flowrate(2),1000);
A = 40.01093;
mesh = [];
for i = 1:numel(x)
    for j = 1:numel(q)
        pm = 28*(MaxPerc-x(i))/100 + 44*x(i)/100+ 32*O2perc/100 + 18*H2Operc/100 + 40*Arperc/100;
        Totmol = q(j)*1000/pm;
        CO2mol = Totmol*x(i)/100;
        A = (25.9106 + O2perc/x(i) *17.4846 + (MaxPerc-x(i))/x(i) *16.4599 +H2Operc/x(i) *60.3696 + 1/1.0797 *30.4728...
            +1/0.2082 + 1/(1.0797*7.6286) *(11.3696+72.7773) + 6*1/(1.0797*7.6286) *(-14.2089) + 1/(1.0797*7.6286*(-0.026592))...
            + 6/(1.0797*7.6286*1.6875*(-0.008673)) + 1/(1.0797*7.6286*3.096*(-0.01910)))*3600*8760*10^-9 *1/10.11 *10^-3 *1.2;
        a = A/(1.2*365) * 0.1 * 10^6;%Tonnellate/day di biogas necessario; utilizzato per dimensionamento ASU.
        Maqua = ((CO2mol/0.0252)*0.9*10^3)/(35.22*490);
        Wmec = Maqua*59.9238*0.90*0.60;%potenza meccanica turbina: Macqua*Lavoro*efficienzaturb*effcompressore
        b = 4125*(Wmec)^0.7 *(1+(0.05/(1-0.4))^3)*(1+exp((813-866)/10.42))*0.1*10^-6;
        B = 6/(1.0797*7.6286*0.0079)*1/0.80*1500*0.101*10^-12;
        C = 2*pi*3*8*994*0.1*10^-6;%densità acciaio 8000kg/m3;semplificato in 8 t/m3.
        D = 13140*(CO2mol/0.2082)^0.67 *10^-6 * 0.1;
        mesh(i,j) = (A+B) * CO2mol + 22*((a*CO2mol)/432)^0.6 + b + C*sqrt(q(j)/(pi*3))+D-(CO2mol/(1.0797*7.6286) *3600*56*8760*25*10^-12);
    end
end

y= mesh;
surf(x,q,mesh,LineStyle=":")
end