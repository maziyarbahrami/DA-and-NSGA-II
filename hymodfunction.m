%Non-linear rainfall-runoff model
%m:inpot vector contain 14 columns
%parameters (5 ones), state variables (5 ones) , precipitation and
%pet;PE=0;EA=0;
%parameters definition (5 ONES)
%m(1):Cmax (Maximum Soil Moisture - storage-  Capacity)(mm)
%m(2):BEXP (degree of spatial variability parameter of soil moisture)
%m(3):ALPHA (the partitioning factor between two series of reservoir tanks)
%m(4):RS (the residence for the time slow-flow tank)(1/day)
%m(5):RQ(the residence for the time quick-flow tank) (1/day)
%state variables definition (5 ONES)
%m(6):SM (Soil Moisture)
%m(7):SF1 (Fast Reservoir 1)
%m(8):SF2 (Fast Reservoir 2)
%m(9):SF3 (Fast Reservoir 3)
%m(10):SS (Slow Reservoir)
%m(11):EA=0;
%m(12):PE=0;
%13th column:precipitation
%14th column:pet
%OUTPUT definitios
%y:total channel flow
%qf3:surface flow
%qs:ground flow
function [y,s1,s2,s3,s4,s5,qs,qf3]=hymod(m)

m(6)=m(6)+m(13)-m(12)-m(11);
m(6)=max(m(6),0);
p_ex=m(6)-m(1);
p_ex=max(p_ex,0);
m(6)=m(6)-p_ex;
p_ef=m(13)-p_ex;
w=min(m(6)/m(1),1);
m(11)=w*m(14);
m(11)=min(m(11),m(6));
def_ea=max(m(14)-m(11),0);
f=1-(1-w)^m(2);
m(12)=f*p_ef;
ea2=min(p_ex,def_ea);
p_ex=p_ex-ea2;
m(11)=m(11)+ea2;
m(7)=m(7)+p_ex+m(3)*m(12)-m(5)*m(7);
m(7)=max(m(7),0);
qf1=m(5)*m(7);
m(8)=m(8)+qf1-m(5)*m(8);
m(8)=max(m(8),0);
qf2=m(5)*m(8);
m(9)=m(9)+qf2-m(5)*m(9);
m(9)=max(m(9),0);
qf3=m(5)*m(9);
m(10)=m(10)+(1-m(3))*m(12)-m(4)*m(10);
m(10)=max(m(10),0);
qs=m(4)*m(10);
y=qs+qf3;
s1=m(6);
s2=m(7);
s3=m(8);
s4=m(9);
s5=m(10);


end