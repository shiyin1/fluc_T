%clear
%load('HRG.dat')
T=HRG(:,1);
p=HRG(:,2).*T.^4;;
eps=HRG(:,3).*T.^4;
s=HRG(:,4).*T.^3;
smid=(s(2:end)+s(1:end-1))/2;
depsdT=((eps(2:end)-eps(1:end-1)));
Cs=smid./depsdT;
dpdT=((p(2:end)-p(1:end-1)));
plot(T(50:150),dpdT(50:150))
hold on
plot(T(50:150),s(50:150),'r')
hold off