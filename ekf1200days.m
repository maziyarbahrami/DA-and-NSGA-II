%Integrated multi-source data assimilation and NSGA-II multi objective
%optimization framework for streamflow simulation
%Ensemble Kalman Filter (EKF) data assimilation code
%Catchment: Siakh-darengon
%Area: 1061.8 km^2
clc
clear all
close all
tic
load('baseflow', 'baseflow')
load('fulltimepr','fulltimepr')
load('fulltimerunoff','fulltimerunoff')
load('fulltimetabkhir','fulltimetabkhir')
load('optimizedparameters','optimizedparameters')
load('optimizedstates','optimizedstates')
%change the unit of runoff from m^3/s to the mm
fulltimerunoffmm=(fulltimerunoff*3600*24/(1061.8*10^6))*1000;
fulltimerunoffmm=fulltimerunoff;
%%
%number of particles:100
particleset=zeros(100,14);
%%first day simulation
particleset(:,1)=optimizedparameters(1);
particleset(:,2)=optimizedparameters(2);
particleset(:,3)=optimizedparameters(3);
particleset(:,4)=optimizedparameters(4);
particleset(:,5)=optimizedparameters(5);
particleset(:,6)=optimizedstates(1);
particleset(:,7)=optimizedstates(2);
particleset(:,8)=optimizedstates(3);
particleset(:,9)=optimizedstates(4);
particleset(:,10)=optimizedstates(5);
particleset(:,11)=0;
particleset(:,12)=0;
particleset(:,13)=fulltimepr(1);
particleset(:,14)=fulltimetabkhir(1);
%apply noise on the states
for i=6:10
    for j=1:100
        particleset(j,i)=particleset(j,i)+normrnd(0,0.15*particleset(j,i));
    end
end
meanv=mean(particleset(:,6));
%apply noise on the forcing
for i=13:14
    for j=1:100
        particleset(j,i)=particleset(j,i)+normrnd(0,0.15*particleset(j,i));
    end
end
for i=1:100
    [tot(i),s1(i),s2(i),s3(i),s4(i),s5(i)]=hymodfunction(particleset(i,:));
end
state1vector=s1(:);
state2vector=s2(:);
state3vector=s3(:);
state4vector=s4(:);
state5vector=s5(:);
for i=1:100
    obsrunvector(i,1)=fulltimerunoffmm(1);
end
for i=1:100
    obsrunvector(i,1)=obsrunvector(i,1)+normrnd(0,0.15*obsrunvector(i,1));
end
simrunvector(:,1)=tot(:);
for i=1:100
    simrunvector(i,1)=simrunvector(i,1)+normrnd(0,0.15*simrunvector(i,1));
end
%change streamflow units to mm
obsrunvectormm=(obsrunvector*3600*24/(1061.8*10^6))*1000;
simrunvectormm=(simrunvector*3600*24/(1061.8*10^6))*1000;
%%
a1=cov(state1vector,simrunvector);
a2=cov(state2vector,simrunvector);
a3=cov(state3vector,simrunvector);
a4=cov(state4vector,simrunvector);
a5=cov(state5vector,simrunvector);
b=cov(simrunvector,obsrunvector);
c=var(obsrunvector);
kalmangain1=a1(1,2)/(b(1,2)+c);
kalmangain2=a2(1,2)/(b(1,2)+c);
kalmangain3=a3(1,2)/(b(1,2)+c);
kalmangain4=a4(1,2)/(b(1,2)+c);
kalmangain5=a5(1,2)/(b(1,2)+c);
bb=obsrunvector-simrunvector;

correctedstate1=state1vector+kalmangain1*(obsrunvectormm-simrunvectormm);
correctedstate2=state2vector+kalmangain2*(obsrunvectormm-simrunvectormm);
correctedstate3=state3vector+kalmangain3*(obsrunvectormm-simrunvectormm);
correctedstate4=state4vector+kalmangain4*(obsrunvectormm-simrunvectormm);
correctedstate5=state5vector+kalmangain5*(obsrunvectormm-simrunvectormm);
% figure
% hist(simrunvector)
particleset2=particleset;
particleset2(:,6)=correctedstate1(:);
particleset2(:,7)=correctedstate2(:);
particleset2(:,8)=correctedstate3(:);
particleset2(:,9)=correctedstate4(:);
particleset2(:,10)=correctedstate5(:);
for i=1:100
    [tot2(i),s1(i),s2(i),s3(i),s4(i),s5(i)]=hymodfunction(particleset2(i,:));
end
simrunvector2(:,1)=tot2(:);
% figure
% hist(simrunvector2)
% hold on
% plot(0.129,0)
mean1=mean(tot);
absdif1=abs(obsrunvector-simrunvector);
absdif2=abs(obsrunvector-simrunvector2);
mab1=sum(absdif1)/100;
mab2=sum(absdif2)/100;
mean2=mean(tot2);
aa=simrunvector2-simrunvector;
particleset=particleset2;
%%
priorstatemat=zeros(100,5,4383);
poststatemat=zeros(100,5,4383);
currun=zeros(1,4383);
currunnew=zeros(1,4383);
a1time=zeros(1,4383);
timesimrun=zeros(100,4383);
timeobsrun=zeros(100,4383);

for t=1:1200
    if fulltimepr(t)>=10
        particleset(:,6)=optimizedstates(1);
        particleset(:,7)=optimizedstates(2);
        particleset(:,8)=optimizedstates(3);
        particleset(:,9)=optimizedstates(4);
        particleset(:,10)=optimizedstates(5);
        particleset(:,13)=fulltimepr(t);
        particleset(:,14)=fulltimetabkhir(t);
        %%
        %apply noise on the states
        t
        for i=6:10
            for j=1:100
                particleset(j,i)=particleset(j,i)+normrnd(0,0.15*particleset(j,i));
            end
        end
        
        %apply noise on the forcing
        for i=13:14
            for j=1:100
                particleset(j,i)=particleset(j,i)+normrnd(0,0.15*particleset(j,i));
            end
        end
        %Run Hymod for particles
        for i=1:100
            [tot(i),s1(i),s2(i),s3(i),s4(i),s5(i)]=hymodfunction(particleset(i,:));
        end
        runoff(:,t)=tot(:);
        state1vector=s1(:);
        state2vector=s2(:);
        state3vector=s3(:);
        state4vector=s4(:);
        state5vector=s5(:);
        %%
        timestate1(:,t)=state1vector;
        timestate2(:,t)=state2vector;
        timestate3(:,t)=state3vector;
        timestate4(:,t)=state4vector;
        timestate5(:,t)=state5vector;
        timepar1(:,t)=particleset(:,1);
        prerun(:,t)=tot;
        aaa=mean(tot);
        if fulltimepr(t)>=10
            currunnew(1,t)=0.040*mean(tot);
            
            d=1;
            for s=1:80
                d=d+1;
                currunnew(d,t+s)=-0.000297*s^3+0.05097*s^2-2.826*s+57.29;
            end
        end
        if fulltimepr(t) ~=0
            for i=1:100
                obsrunvector(i,1)=fulltimerunoffmm(t)-currun(t);
                if obsrunvector(i,1)<=0
                    obsrunvector(i,1)=0;
                end
            end
            
        end
        if fulltimepr==0
            obsrunvector(:,1)=0;
        end
        for i=1:100
            obsrunvector(i,1)=obsrunvector(i,1)+normrnd(0,0.15*obsrunvector(i,1));
            if obsrunvector(i,1)==0
                obsrunvector(i,1)=obsrunvector(i,1)+normrnd(0,0.15);
            end
        end
        timeobser(:,t)=obsrunvector(:,1);
        simrunvector(:,1)=tot2(:);
        for i=1:100
            simrunvector(i,1)=simrunvector(i,1)+normrnd(0,0.15*simrunvector(i,1));
        end
        timesimrun(:,t)=simrunvector;
        a1=cov(state1vector,simrunvector,'omitrows');
        a2=cov(state2vector,simrunvector,'omitrows');
        a3=cov(state3vector,simrunvector,'omitrows');
        a4=cov(state4vector,simrunvector,'omitrows');
        a5=cov(state5vector,simrunvector,'omitrows');
        b=cov(simrunvector,obsrunvector,'omitrows');
        c=var(obsrunvector);
        a1time(1,t)=a1(1,2);
        kalmangain1=a1(1,2)/(b(1,2)+c);
        timekal1(t)=kalmangain1;
        kalmangain2=a2(1,2)/(b(1,2)+c);
        kalmangain3=a3(1,2)/(b(1,2)+c);
        kalmangain4=a4(1,2)/(b(1,2)+c);
        kalmangain5=a5(1,2)/(b(1,2)+c);
        bb=obsrunvector-simrunvector;
        correctedstate1=state1vector+kalmangain1*(obsrunvector-simrunvector);
        dif(:,t)=obsrunvector-simrunvector;
        correctedstate2=state2vector+kalmangain2*(obsrunvector-simrunvector);
        correctedstate3=state3vector+kalmangain3*(obsrunvector-simrunvector);
        correctedstate4=state4vector+kalmangain4*(obsrunvector-simrunvector);
        correctedstate5=state5vector+kalmangain5*(obsrunvector-simrunvector);
        %%
        particleset(:,6)=correctedstate1;
        particleset(:,7)=correctedstate2;
        particleset(:,8)=correctedstate3;
        particleset(:,9)=correctedstate4;
        particleset(:,10)=correctedstate5;
        %%
        postcorrstatemat(:,1,t)=particleset(:,6);
        postcorrstatemat(:,2,t)=particleset(:,7);
        postcorrstatemat(:,3,t)=particleset(:,8);
        postcorrstatemat(:,4,t)=particleset(:,9);
        postcorrstatemat(:,5,t)=particleset(:,10);
        %%
        
        for i=1:100
            [tot(i),s1(i),s2(i),s3(i),s4(i),s5(i)]=hymodfunction(particleset(i,:));
        end
        particleset(:,6)=s1(:);
        particleset(:,7)=s2(:);
        particleset(:,8)=s3(:);
        particleset(:,9)=s4(:);
        particleset(:,10)=s5(:);
        %save posterior statevariable matrix
        poststatemat(:,1,t)=particleset(:,6);
        poststatemat(:,2,t)=particleset(:,7);
        poststatemat(:,3,t)=particleset(:,8);
        poststatemat(:,4,t)=particleset(:,9);
        poststatemat(:,5,t)=particleset(:,10);
        %%
        poststate1matvec(:,t)=correctedstate1(:);
        poststate2matvec(:,t)=correctedstate2(:);
        poststate3matvec(:,t)=correctedstate3(:);
        poststate4matvec(:,t)=correctedstate4(:);
        poststate5matvec(:,t)=correctedstate5(:);
        
        postrun(:,t)=tot;
        aaa=mean(tot);
        finalsim(t)=mean(tot);
        if fulltimepr(t)>=450
            currunnew(1,t)=0.040*mean(tot);
            
            d=1;
            for s=1:80
                d=d+1;
                currunnew(d,t+s)=-0.000297*s^3+0.05097*s^2-2.826*s+57.29;
                if s==1
                    currunnew(d,t+s)=currunnew(d,t+s)+0.6*aaa;
                end
            end
        end
    end
end


toc
