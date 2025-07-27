clear;
V=load('VTOTAL.DAT');
T=load('TMEV.DAT');
p=(V(1)-V)./T.^4;
plot(p);