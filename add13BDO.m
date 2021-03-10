%Script to add 1,3-BDO production into E. coli
%Base model = iML1515
%Pathway based on Kataoka et al (2013) - Improvement of (R)-1,3-butanediol 
%production by engineered Escherichia coli

function [model] = add13BDO(input)
%add necessary metabolites
model = addMetabolite(input,'3hbald_c','3-hydroxybutyraldehyde','C4H8O2');
model = addMetabolite(model,'13bdo_c','1,3-butanediol','C4H10O2');
model = addMetabolite(model,'13bdo_e','1,3-butanediol extracellular','C4H10O2');

%add necessary reactions, starts from 3-hydroxybutyryl-CoA (3hbcoa_c)
model = addReaction(model,'AdhE',{'3hbcoa_c', 'nadh_c','3hbald_c','nad_c', 'coa_c'},...
    [-1 -1 1 1 1],true); %dehydrogenase
model = addReaction(model,'AdhE_bdo',{'3hbald_c', 'nadh_c', '13bdo_c', 'nad_c'},...
    [-1 -1 1 1], true); %dehydrogenase
model = addReaction(model, '13BDOtex',{'13bdo_c', '13bdo_e'},...
    [-1 1], true); %1,3-BDO transport
model = addExchangeRxn(model, {'13bdo_e'}); %1,3-BDO exchange rxn
end
