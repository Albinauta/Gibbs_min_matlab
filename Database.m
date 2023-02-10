classdef Database

properties

filename
sheetname
range
end
properties(Dependent)
    Table
end
methods
    %%%%METODO COSTRUTTORE%%%
    function thisDatabase = Database(filename,sheetname,range)
        thisDatabase.filename = filename;
        thisDatabase.sheetname = sheetname;
        thisDatabase.range = range;
    end
    %%%%FUNZIONE ESTRAI PROPRIETA' %%%%%%%%
    %data una formula fornisce le proprietà termochimiche associate,
    %estraendole dalla tabella
    function vec = getVec(thisDatabase,species_name,T)
        Names = thisDatabase.Table.Species;
        matchnumbers = find(Names == species_name);
        if numel(matchnumbers) == 0
            vec = 0;
            return
        end
        vec = [];
        for n = 1 : length(matchnumbers)%prima ricerca, itera sugli indici trovati in base alla stringa fornita
            indice = matchnumbers(n);
            if T >= thisDatabase.Table{indice,2} && T <= thisDatabase.Table{indice,3}
                %if species_name == Names(n)
                vec = thisDatabase.Table{indice,:};
                return
            end
        end
        if numel(vec) == 0%seconda ricerca, itera dal composto con i limiti dell'intervallo di T più basso verso quello con i limiti più alti(che per costruzione del file excel si trova "dopo"..
            if T > thisDatabase.Table{matchnumbers(end),2}
                n0 = matchnumbers(end) + 1;
                while n0 <= numel(Names)
                    if extract(Names(n0),1) == extract(species_name,1) && strlength(Names(n0)) == strlength(species_name)
                        if T >= thisDatabase.Table{n0,2} && T <= thisDatabase.Table{n0,3}
                            vec = thisDatabase.Table{n0,:};
                            %n0 = numel(Names) +1;
                            return
                        else
                            n0 = n0 +1;
                        end
                    else
                        n0 = n0 +1;
                    end
                end
                vec = thisDatabase.Table{matchnumbers(end),:};
                return
            else % terza ricerca, questa volta itera verso l'alto, l'effetto finale è che verrà restituito sempre un composto, che al limite avrà come intervalli di T quelli più estremi(con i limiti più piccoli o più grandi)
                n0 = matchnumbers(1) - 1;
                while n0 >= 1
                    if extract(Names(n0),1) == extract(species_name,1) && strlength(Names(n0)) == strlength(species_name)
                        if T >= thisDatabase.Table{n0,2} && T <= thisDatabase.Table{n0,3}
                            vec = thisDatabase.Table{n0,:};
                            %n0 = numel(Names) +1;
                            return
                        else
                            n0 = n0 - 1;
                        end
                    else
                        n0 = n0 - 1;
                    end
                end
                vec = thisDatabase.Table{matchnumbers(1),:};
                return
            end
        end
        end 

%%%%% METODO COSTRUZIONE DELLA TABELLA%%%%%%%
    function Table = get.Table(theDatabase)

        if nargin == 1 || isempty(string(theDatabase.sheetname))
            sheetName = 1;
        end

        % If row start and end points are not specified, define defaults
        if length(theDatabase.range) <=2
            Range = [6, 57];%  IN CASO DI MODIFICA DEL FILE EXCEL, BISOGNA VARIARE QUESTI LIMITI
        else
            Range = theDatabase.range;
        end

        %% Set up the Import Options and import the data
        opts = spreadsheetImportOptions("NumVariables", 21);

        % Specify sheet and range
        opts.Sheet = sheetName;
        opts.DataRange = "C" + Range(1, 1) + ":W" + Range(1, 2);

        % Specify column names and types
        opts.VariableNames = ["Species", "T1", "T2", "A", "B", "(C)", "D", "E", "F", "G", "(H)", "Hf", "Sf", "Exf", "C", "O", "H", "N", "S", "Ca", "Fe"];
        opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
        opts.VariableNamingRule = ["preserve"];
        % Import the data
        Table = readtable(string(theDatabase.filename), opts, "UseExcel", false);
    end
end

end