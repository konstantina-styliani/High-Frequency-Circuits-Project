function [val]=Ask1_4fun(p)
clear all;
clc;
N=150;
f0=3*10^8;
f = [(0.01)*f0:0.01*f0:2*f0]; 
normf=[0.01:0.01:2];
p=[0.5000,	0.0100,	0.0100,	0.5000,	0.0100,	0.0909];  % optimized p=[d1,d2,d3,l1,l2,l3]

Zin3=zeros(length(f),1);
ZinC=zeros(length(f),1);
ZL2=zeros(length(f),1);
Zin2=zeros(length(f),1);
ZinB=zeros(length(f),1);
ZL1=zeros(length(f),1);
Zin1=zeros(length(f),1);
Zin=zeros(length(f),1);

Zo=50;
ZL=120+(1i)*60;
bd1=bl(normf,p(1));
bd2=bl(normf,p(2));
bd3=bl(normf,p(3));
bl1=bl(normf,p(4));
bl2=bl(normf,p(5));
bl3=bl(normf,p(6));



for k=1:length(f) %to mhkos-h diastash tou f .Κάνουμε τον υπολογισμο για καθε μια συχνοτητα ξεχωριστά

Zin1(k)=Zinfun(ZL,Zo,bd1(k));

ZinA(k)=Zbranch(Zo,bl1(k));

ZL2(k)=Zpar(ZinA(k),Zin1(k));

Zin2(k)=Zinfun(ZL2(k),Zo,bd2);

ZinB(k)=Zbranch(Zo,bl2(k));

ZL1(k)=Zpar(ZinB(k),Zin2(k));

Zin3(k)=Zinfun(ZL1(k),Zo,bd3(k));

ZinC(k)=Zbranch(Zo,bl3(k));

Zin(k)=Zpar(Zin3(k),ZinC(k)); 

Ref(k)=abs((Zin(k)-Zo)/(Zin(k)+Zo));
end
mean_Ref=mean(Ref);

function [val]= bl(f,l)
val=(l.*2*f.*pi);%/(3*10^8) ;
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
[val]=mean_Ref
plot(f,Ref)
xlabel('freq')
ylabel('Ref')
end
