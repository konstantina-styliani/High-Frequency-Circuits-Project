clear all;
clc;
N=200;
f0=10^9;
f = 0:(4*f0/N):(4*f0); 
Zin1=zeros(length(f),1); %δίνουμε διάσταση , 1 στήλης
Zin2=zeros(length(f),1);
Zin3=zeros(length(f),1);
ZB=zeros(length(f),1);
Ref=zeros(length(f),1);  %δημιουργούμε θέσεις 
ZL=100;
Zo=50;
C=2*10^(-12);
%length = το πλήθος του διανύσματος
bl3=bl(f,0.32);
bl2=bl(f,0.24);
bl1=bl(f,0.1);
for k=1:length(f)

Zin3(k)=Zinfun(ZL,Zo,bl3(k));

ZB(k)=Zin3(k)-(1i./(2*pi*f(k)*C));

Zin2(k)=Zinfun(ZB(k),Zo,bl2(k));

Zin1(k)=Zo/(1i*tan(bl1(k)));

%Zin1(k)=Zbranch(Zo,bl1(k));

Zin(k)=(Zin1(k)*Zin2(k))/(Zin1(k)+Zin2(k));
Ref(k)=abs((Zin(k)-Zo)/(Zin(k)+Zo));
Ref_dB(k)=20*log10(Ref(k));

if Ref_dB(k)<-25
    Ref_dB(k)=-25;
end
end

for k= 1:length(f)
    y(k)= 0.316;
    y_dB(k)=-10;
end

figure(1)
plot(f(2:length(f)),Ref(2:length(f)),'r+',f(2:length(f)),y(2:length(f)),'b--')
grid on; %melimetre
xlabel('freq');
ylabel('Reflection-Coefficient');

figure(2)
plot(f(2:length(f)),Ref_dB(2:length(f)),'r|',f(2:length(f)),y_dB(2:length(f)),'b--')
grid on;
xlabel('freq');
ylabel('Reflection-Coefficient (dB)')


function [val]= bl(f,l)
val=(l*2*f.*pi)/(10^9) ;
end 

function [val]= Zinfun(ZL,Zo,bl_)
val= (Zo*(ZL+1i*Zo*tan(bl_)))/(Zo+1i*ZL.*tan(bl_));
end

function [val]= Zbranch(Zo,bl)
val=Zo./(1i*tan(bl));
end
