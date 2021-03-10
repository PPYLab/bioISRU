%Script to add 2,3-BDO production into E. coli
%Base model = 
%Pathway based on Xu et al (2014) -  Systematic metabolic engineering of 
%E. coli for high-yield production of fuel bio-chemical 2,3-butanediol

function [model] = add23BDO(input)
%add necessary metabolites
model = addMetabolite(input,'acetoin_c','R-acetoin','C4H8O2');
model = addMetabolite(model,'23bdo_c','2,3-butanediol','C4H10O2');
model = addMetabolite(model,'23bdo_e','2,3-butanediol extracellular','C4H10O2');

%add necessary reactions, starts from ?-acetolactone (alac__S_c)
model = addReaction(model,'ALDC',{'alac__S_c','acetoin_c','co2_c'},...
    [-1 1 1],false); %?-acetolactate decarboxylase
model = addReaction(model,'BDD',{'acetoin_c', 'nadh_c', '23bdo_c', 'nad_c'},...
    [-1 -1 1 1], true); %2,3-butanediol dehydrogenase
model = addReaction(model, '23BDOtex',{'23bdo_c', '23bdo_e'},...
    [-1 1], true); %2,3-BDO transport
model = addExchangeRxn(model, {'23bdo_e'}); %2,3-BDO exchange rxn
end
