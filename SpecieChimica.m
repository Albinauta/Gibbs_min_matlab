classdef SpecieChimica < matlab.mixin.SetGet & dynamicprops %%%L'oggetto è una handle class, ciò mi permette di modificarne i valori senza creare copie indipendenti di esso
    

properties
Name
MassFlux
Temp
Inert% definisce se l'oggetto si comporta come un inerte o meno. importante per la funzione di iterazione in condizioni adiabatiche
Properties
DG0 = 0;
DG = 0;
end
properties(Dependent)
    %H
    %S
    Comp
end

methods
%%%%%%%%%%%%%METODO COSTRUTTORE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function thisSpecieChimica = SpecieChimica(Name,Massflux,Temp,Inert,p)
        thisSpecieChimica.Name = Name;
        thisSpecieChimica.MassFlux = Massflux;
        thisSpecieChimica.Temp = Temp;
        thisSpecieChimica.Inert = Inert;
        thisSpecieChimica.Properties = p;
    end
    %%%METODO DI COPIA OGGETTO%%%%%% crea un oggetto indipendente da quello
    %%%fornito, ai fini di iterazione e modifica
    function new_obj = copy(obj)
        nome = obj.Name;
        mass = obj.MassFlux;
        temp = obj.Temp;
        inert = obj.Inert;
        p = obj.Properties;
        new_obj = SpecieChimica(nome,mass,temp,inert,p);
    end
 %%%%%%% METODO DI CALCOLO DEL DELTA G IN CONDIZIONI STANDARD A T=298.2 K %%%
%%%%obsoleto
    function calcDG0(obj,database)
        elements = ["C","O2","H2","N2","S","Caα","Feγ"];
        fstech = obj.Comp;
        fstech(2) = fstech(2) * 0.5;
        fstech(3) = fstech(3) * 0.5;
        fstech(4) = fstech(4) * 0.5;
        DSr = 0;
        for i = 1 : length(obj.Comp)
            v = database.getVec(elements(i),298.2);
            DSr = DSr +(fstech(i)*str2num(v(13)));
        end
        DS = str2num(obj.Properties(13)) - DSr;
        DS = DS/1000;
        DG01 = str2num(obj.Properties(12)) - 298.2*DS;
        obj.DG0 = DG01;   
    end
 %%%% METODO CALCOLO DG CONDIZIONI STANDARD A T QUALSIASI%%%%%%%%
    function calcDG(obj,database)
        DHt = obj.intgH(database,obj.Temp);
        DSt = obj.intgS(database,obj.Temp);
        t = obj.Temp/1000;
        DG1 = DHt -(t*DSt);
        obj.DG = DG1;
    end
%%% METODO CALCOLO DH_formazione A T GENERICA%%%% la T di partenza è a 25 C
    function y = intgH(obj,database,Tf)
        v = database.getVec(obj.Name,Tf);
        tf = Tf/1000;
        DHt = str2num(v(12)) + str2num(v(4))*tf + str2num(v(5))*(tf^2)/2 + str2num(v(6))*(tf^3)/3 + ...
            str2num(v(7))*(tf^4)/4 - str2num(v(8))/tf + str2num(v(9)) - str2num(v(11));
        y = DHt; %escono in kj
    end
%%%METODO CALCOLO DS_formazione A T GENERICA%%%
    function y = intgS(obj,database,Tf)
        v = database.getVec(obj.Name,Tf);
        tf = Tf/1000;
        y = str2num(v(4))*log(tf) + str2num(v(5))*tf + ((str2num(v(6))*tf^2)/2) + ((str2num(v(7))*tf^4)/4) - (str2num(v(8))/(2*tf^2)) + str2num(v(10)); %escono in joule
    end
%%%METODO PER CALCOLARE IL CALORE SPECIFICO, AI FINI DEL CALCOLO DELLA Teq
%%%DI UNA MISCELA
    function y = getcp(obj,t)
        t = t/1000;
        costanti = obj.Properties(4:8);
        costanti = str2double(costanti);
        T = [1;t;t^2;t^3;1/t^2];
        y = costanti*T;
    end
%%%%METODO CHE FORNISCE LE COSTANTI A-H PER CREARE MATRICE COSTANTI
    function const = getConst(obj)
        %consttemp = double(obj.Properties(4:9));
        %const = [consttemp,obj.Properties(11)];
        const = double(obj.Properties(4:11));
    end
    %function H = get.H(theSpecieChimica)
        %A = str2num(theSpecieChimica.Properties(4));
        %B = str2num(theSpecieChimica.Properties(5));
        %C = str2num(theSpecieChimica.Properties(6));
        %D = str2num(theSpecieChimica.Properties(7));
        %E = str2num(theSpecieChimica.Properties(8));
        %F = str2num(theSpecieChimica.Properties(9));
        %G = str2num(theSpecieChimica.Properties(10));
        %h = str2num(theSpecieChimica.Properties(11));
        %t = theSpecieChimica.Temp/1000;
        %H = A*t + B*(t^2)/2 + C*(t^3)/3 + D*(t^4)/4 - E/t + F -h;
        %return 
    %end
    %function S = get.S(theSpecieChimica)
        %A = str2num(theSpecieChimica.Properties(4));
        %B = str2num(theSpecieChimica.Properties(5));
        %C = str2num(theSpecieChimica.Properties(6));
        %D = str2num(theSpecieChimica.Properties(7));
        %E = str2num(theSpecieChimica.Properties(8));
        %F = str2num(theSpecieChimica.Properties(9));
        %G = str2num(theSpecieChimica.Properties(10));
        %h = str2num(theSpecieChimica.Properties(11));
        %t = theSpecieChimica.Temp/1000;
        %S = A * log(t) + B*t + C*(t^2)/2 + D*(t^3)/3 -E/(2*t^2) + G;
        %return
    %end
%%%% METODO DI CALCOLO COMPOSIZIONE DELLA SPECIE CHIMICA %%%%%
    function Comp = get.Comp(theSpecieChimica)
        C = str2num(theSpecieChimica.Properties(15));
        O = str2num(theSpecieChimica.Properties(16));
        H = str2num(theSpecieChimica.Properties(17));
        N = str2num(theSpecieChimica.Properties(18));
        S = str2num(theSpecieChimica.Properties(19));
        Ca = str2num(theSpecieChimica.Properties(20));
        Fe = str2num(theSpecieChimica.Properties(21));
        Comp = [C,O,H,N,S,Ca,Fe]';
        return
    end
end


end 