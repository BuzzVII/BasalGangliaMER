
bar=waitbar(0,'Repeat simulations');

percenters=[0.005, 0.05, 0.2, 0.5, 0.7, 0.9];
Niterations=20;
plen=length(percenters);

meanS=zeros(1,plen);
maxS=zeros(1,plen);
varianceS=zeros(1,plen);
StoN=zeros(1,plen);
StoNvar=zeros(1,plen);
maxStoN=zeros(1,plen);
j=1;

MaxSignal=zeros(1,Niterations);
RMSnoise=zeros(1,Niterations);
SigToNoise=zeros(1,Niterations);

for percent = percenters
    for i = 1:Niterations
        waitbar(((j-1)*Niterations+i)/(plen*Niterations),bar)
        [t Vout comp_time]=pprateshapesych(0.8,10000,1,10000*percent);
        MaxSignal(i)=max(Vout);
        RMSnoise(i)=sqrt(mean(Vout.^2));
        SigToNoise(i)=MaxSignal(i)/RMSnoise(i);
    end
    meanS(j)=mean(MaxSignal);
    varianceS(j)=var(MaxSignal);
    StoN(j)=mean(SigToNoise);
    StoNvar(j)=var(SigToNoise);
    maxS(j)=max(MaxSignal);
    maxStoN(j)=max(SigToNoise);
    j=j+1;
end

close(bar)


figure(1)
plot(percenters,meanS,'x','markersize',10,'linewidth',10);
h=gca;
set(h,'fontsize',16)
xlabel('% synchronized')
ylabel('Mean Peak Signal (mV)')

figure(2)
plot(percenters,varianceS,'x','markersize',10,'linewidth',10);
h=gca;
set(h,'fontsize',16)
xlabel('% synchronized')
ylabel('Variance of Peak Signal (mV)')