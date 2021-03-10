%Script to add 1,2-PDO production into E. coli
%Base model = iML1515
%Pathway based on Niu et al (2019) - Metabolic engineering of Escherichia coli 
%for the de novo stereospecific biosynthesis of 1,2-propanediol through lactic acid

%most of the compounds are endogenous, add bypass step to remove
%methylglyoxyl buildup to avoid cytotoxicity

function [model] = add12PDO(input)
%add necessary metabolites
model = addMetabolite(input,'lac_coa_c','lactoyl-coa','C24H40N7O18P3S');

%add necessary reactions, starts from D-lactate (lac__D_c) - just need
%to add connection of lactic acid to lactaldehyde (lald__D_c) through new
%laconyl-coa
model = addReaction(model,'lacCoA_trans',{'lac__D_c', 'accoa_c','lac_coa_c','ac_c'},...
    [-1 -1 1 1],true); %coa transferase
model = addReaction(model,'lactald_dehyd',{'lac_coa_c', 'nadh_c', 'lald__D_c', 'nad_c', 'coa_c'},...
    [-1 -1 1 1 1], true); %dehydrogenase

%remove reactions to fully bypass methylglyoxal production
model = removeRxns(model, {'MGSA','AACTOOR'})
end
