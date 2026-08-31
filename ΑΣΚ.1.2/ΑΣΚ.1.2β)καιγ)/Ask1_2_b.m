clear all;
clc;
N=150;
f0=10^9;
f=0:(3*f0/N):(3*f0);
Zin3=zeros(length(f),1);
ZinC=zeros(length(f),1);
ZL2=zeros(length(f),1);
Zin2=zeros(length(f),1)
ZinB=zeros(length(f),1);
ZL1=zeros(length(f),1);
Zin1=zeros(length(f),1);
Zin=zeros(length(f),1);


Zo=50;
ZL=Zo;
ZoA=101.6 ;
ZoC=98.45;
ZoB=43.6;
bli=bl(f);

for k=1:length(f) %to mhkos-h diastash tou f .Κάνουμε τον υπολογισμο για καθε μια συχνοτητα ξεχωριστά

Zin3(k)=Zinfun(ZL,Zo,bli(k));

ZinC(k)=Zbranch(ZoC,bli(k));

ZL2(k)=Zpar(ZinC(k),Zin3(k));

Zin2(k)=Zinfun(ZL2(k),ZoA,bli(k));

ZinB(k)=Zbranch(ZoB,bli(k));

ZL1(k)=Zpar(ZinB(k),Zin2(k));

Zin1(k)=Zinfun(ZL1(k),ZoA,bli(k));

Zin(k)=Zpar(Zin1(k),ZinC(k)); 

Ref(k)=abs((Zin(k)-Zo)/(Zin(k)+Zo));
Ref_dB(k)=20*log10(Ref(k));

if Ref_dB(k)<-60
   Ref_dB(k)=-60;
end

SWR(k)=(1+Ref(k))/(1-Ref(k));
if SWR(k)>10
    SWR(k)=10;
end
end

yy=2*ones(length(f),1);
kk=-10*ones(length(f),1);

figure(1)
plot(f,Ref_dB,'r',f(2:length(f)),kk(2:length(f)),'b--')
grid on;
xlabel('freq')
ylabel('Reflection-Coefficient (dB)')

 figure(2)
 plot(f(2:length(f)),SWR(2:length(f)),'r',f(2:length(f)),yy(2:length(f)),'b--');
 xlabel('freq')
 ylabel('SWR')
 
function [val]= bl(f)
val= (0.25*f.*pi)/(10^9);
end

function [val]= Zinfun(ZL,Zo,bl_)
val= (Zo*(ZL+1i*Zo*tan(bl_)))/(Zo+1i*ZL.*tan(bl_));
end

function [val]=Zpar(ZinA,ZinB)
val=(ZinA.*ZinB)./(ZinA+ZinB);
end

function [val]= Zbranch(Zoi,bli)
val=Zoi./(1i*tan(bli));
end

