clc;
clear all;
N=200;
f0=10^9;
f = 0:(2*f0/N):2*f0 ;
Zo=50;
ZL=10+(1i)*15;
Zg=50-(1i)*40;

Zin3=zeros(length(f),1);
ZinA=zeros(length(f),1);
ZL2=zeros(length(f),1);
Zin2=zeros(length(f),1);
Zin1=zeros(length(f),1);

bl1=bl(f,0.139);
bl2=bl(f,0.1);
bl3=bl(f,0.04);
C1=pykn(f,(2.7*(10^(-12))));
%C1=pykn(f,144.68631*0^(-12)/50);
C2=pykn(f,4.93*(10^(-12)));
%C2=pykn(f,222.8169*10^(-12)/50);
C3=pykn(f,3.08*(10^(-12)));
%C3=pykn(f,148.014*10^(-12)/50);

for k=1:length(f)
Zin1(k)=Zinfun(ZL,Zo,bl1(k))+ C1(k); % pyknwths se seira sthn eisodo

ZL2(k)=Zpar(ZL,C2(k));
Zin2(k)=Zinfun(ZL2(k),Zo,bl2(k));    % pyknwths parallhla sto fortio

ZinA(k)=Zinfun(ZL,Zo,bl3(k));
Zin3(k)=Zpar(ZinA(k),C3(k));         % pyknwths parallhla sthn eisodo

%Rin1(k)=Rin(Zin1(k));
%Rin2(k)=Rin(Zin2(k));
%Rin3(k)=Rin(Zin3(k));

PL1(k)=P(Zin1(k),Zg,Zin1(k));
PL2(k)=P(Zin2(k),Zg,Zin2(k));
PL3(k)=P(Zin3(k),Zg,Zin3(k));
end

kk=+4.5*10^(-3)*ones(length(f),1); %Μεταδιδόμενη Ισχύς στο φορτίο 
                                  %τουλάχιστον το 90%  σε σχέση 
                                  %με το αν είχα προσαρμογή.

figure (1)
plot(f,PL1,'r', f,PL2,'b', f,PL3,'m');
hold on
plot(f,kk,'g--')
xlabel('Συχνότητα')
ylabel('PL - Ισχύς στο Φορτίο')



function [val]= Rin(Zin)
val=real(Zin);
end

function [val]= Zinfun(ZL,Zo,bl_)
val= (Zo*(ZL+1i*Zo*tan(bl_)))./(Zo+1i*ZL.*tan(bl_));
end

function [val]=Zpar(ZinA,ZinB)
val=(ZinA.*ZinB)./(ZinA+ZinB);
end

function [val]=pykn(f,c)
val=1./((1i)*2*pi*f.*c) ;
end

function [val]= P(Rin,Zg,Zin)  % H Isxys sto fortio
val=real(Rin)/abs((Zg+Zin)^2);
end

function [val]= bl(f,l)
val=(l.*2*f.*pi)./(10^9) ;
end 