%determines max theoretical yields from glucose for 1,3-propanediol,
%1,3-butanediol and 2,3-butanediol

%uses genome scale model for E. coli - iML1515, Monk et al (2017) - iML1515, 
%a knowledgebase that computes Escherichia coli traits

model = readCbModel('iML1515.mat');

gluc_up = -10;
model = changeRxnBounds(model,{'EX_glc__D_e' 'EX_o2_e'},[gluc_up -15],'l');

model_23bdo = add23BDO(model);
model_23bdo = changeObjective(model_23bdo,'EX_23bdo_e');
solution_23bdo = optimizeCbModel(model_23bdo);

model_13bdo = add13BDO(model);
model_13bdo = changeObjective(model_13bdo, 'EX_13bdo_e');
solution_13bdo = optimizeCbModel(model_13bdo)

model_12pdo_nat = changeObjective(model, 'EX_12ppd__R_e');
solution_12pdo_nat = optimizeCbModel(model_12pdo_nat)

model_12pdo_eng = add12PDO(model);
model_12pdo_eng = changeObjective(model_12pdo_eng, 'EX_12ppd__R_e');
solution_12pdo_eng = optimizeCbModel(model_12pdo_eng);

yield =     [solution_23bdo.f/abs(gluc_up);
            solution_13bdo.f/abs(gluc_up);
            solution_12pdo_nat.f/abs(gluc_up);
            solution_12pdo_eng.f/abs(gluc_up)]
