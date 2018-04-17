clc,clear all;
tic
mu1=5;mu2=10;sigma1=1;sigma2=2;r=-1/2;
xmin=mu1-5*sigma1;xlen=10*sigma1;
ymin=mu2-5*sigma2;ylen=10*sigma2;
mu=[mu1,mu2];sigma=[sigma1,sigma2];
E=[sigma1*sigma1,r*sigma1*sigma2;r*sigma1*sigma2,sigma2*sigma2];
iE=E^(-1);
sdE=sqrt(det(E));

L=50000;
tempa=0;
windL=200;

state=zeros(L,2);
conve=zeros(1,L);
S=[xmin+xlen*rand(),ymin+ylen*rand()];

for i=1:L
    Snext=[xmin+xlen*rand(),ymin+ylen*rand()];
    p=1/(2*pi*sqrt(det(E)))*exp(-0.5*((S-mu)*(E^(-1))*(S-mu)'));
    pnext=1/(2*pi*sdE)*exp(-0.5*((Snext-mu)*iE*(Snext-mu)'));
    alpha=min(1,pnext/p);
    conve(i)=1-alpha;
    if rand()<alpha
        S=Snext;
    end
    state(i,:)=S;
end

x1=1:length(conve);
P=polyfit(x1,conve,1);
x2=1:windL:L;
rtest=zeros(L/windL,1);
for i=1:L/windL
    TempState=state(1:windL*i,:);
    rtest(i)=min(min(corrcoef(TempState)));
end

figure;
plot(state(:,1),state(:,2),'.r')
hold on
plot(state(L,1),state(L,2),'*g',state(1,1),state(1,2),'*b')
xlabel('x');ylabel('y');

figure;
subplot(2,1,1);plot(x2,rtest);
ylim([-1,0]);
title('correlation coefficient');
subplot(2,1,2);plot(x1,polyval(P,x1),'r');
title('Rejection probability fitting');
r=min(min(corrcoef(state)))
toc