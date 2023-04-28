tic
load('totflow','totflow')
load('totflowtr','totflowtr')
load('totflowms','totflowms')
load('totflowcm','totflowcm')
load('totflowper','totflowper')
load('fulltimerunfin','fulltimerunfin')


cost=@muli;
[x,fval] = gamultiobj(cost,5);
options = gaoptimset('PopulationSize',60,'PlotFcns',@gaplotpareto);
toc