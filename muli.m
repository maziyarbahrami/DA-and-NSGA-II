function z=muli(x)
x1 = abs(x(1));
x2 = abs(x(2));
x3 = abs(x(3));
x4 = abs(x(4));
x5 = abs(x(5));
load('baseflowfin','baseflowfin')
load('totflowtr','totflowtr')
load('totflow','totflow')
totflow(124)=44;
load('fulltimerunfin','fulltimerunfin')
load('totflowms','totflowms')
load('totflowcm','totflowcm')
load('totflowper','totflowper')

Train_Ac = fulltimerunfin(1:800);
Test_Ac = fulltimerunfin(801:1200);
cmorph_Train = totflowcm(1:800);
cmorph_Test = totflowcm(801:1200);
ms_Train = totflowms(1:800);
ms_Test = totflowms(801:1200);
per_Train = totflowper(1:800);
per_Test = totflowper(801:1200);
tr_Train = totflowtr(1:800);
tr_Test = totflowtr(801:1200);
gr_Train = totflow(1:800);
gr_Test = totflow(801:1200);
Train_Pr =x1*gr_Train + x2*cmorph_Train + x3*ms_Train + x4*per_Train + x5*tr_Train;
Test_Pr = x1*gr_Test + x2*cmorph_Test + x3*ms_Test + x4*per_Test + x5*tr_Test;
for i = 1:800;
    Z1(i) = abs(Train_Ac(i) - Train_Pr(i))/Train_Ac(i);
end
Z1 = sum(Z1);
Z1 = Z1/800;
for i = 1:400;
    Z2(i) = abs(Test_Ac(i) - Test_Pr(i))/Test_Ac(i);
end
Z2 = sum(Z2);
Z2 = Z2/400;
z = [Z1 Z2]';
end