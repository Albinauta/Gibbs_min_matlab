function obj = createEl(database,species_name,temp,inert_status,varargin)% Funzione che automatizza la costruzione dell'oggetto specie chimica, rendendo il processo più user friendly
if nargin < 5
    mass_flow = 0;
else
    mass_flow = varargin{1};
end
%%% ESTRAZIONE PROPRIETA' CHIMICHE DEL COMPOSTO DAL DATABASE
prop = database.getVec(species_name,temp);
obj = SpecieChimica(species_name,mass_flow,temp,inert_status,prop);
return
end