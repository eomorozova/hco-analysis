% load the data with current steps @ 10 and 20 C
load('C:\Users\moroz\Documents\Katya_DynamicClamp\hco_related_codes\GM_FI_Temperature.mat','datafi')


%% calculate mean spike height
V10mean=[];  V20mean=[]; k=1;
for j = 1:12
    meanspkh10(j)=mean(datafi.spkh{1,j});
    meanspkh20(j)=mean(datafi.spkh{2,j});
end
%
%% Wilcoxon signed rank-sum test for spike height
[p,h] = signrank(meanspkh10,meanspkh20)

%% plot mean spike height
clf, plot([1 2],[meanspkh10; meanspkh20],'.-','markersize',20,'linewidth',1,'color','k'), hold on
hold on, plot([0.9 0.9],[mean(meanspkh10); mean(meanspkh10)],'.','markersize',30,'color','b');
hold on, plot([2.1 2.1],[mean(meanspkh20); mean(meanspkh20)],'.','markersize',30,'color','r');
hold on, errorbar(0.9,mean(meanspkh10),std(meanspkh10),'b','linewidth',2)
hold on, errorbar(2.1,mean(meanspkh20),std(meanspkh20),'r','linewidth',2)
xlim([0.8 2.2]); ylim([0 32])
xticks([1 2]); xticklabels({'10^oC','20^oC'}); %ylim([0 5])
%sigstar([1,2],p)
%text(1.6,23,['(p=',num2str(round(p,4)),')'],'fontsize',16,'FontName','Myriad Pro')
ylabel('Spike amplitude, mV'); box off 
text(1.9,25,'N=12','Fontsize',20,'FontName','Arial')
set(gca,'Fontsize',16,'FontName','Arial')


%% plot F-Is
clf,
col=[[0,0,1];[0.1,0.5,0];[1,0,0]];
ii=1; % experiment
k=1;
for n=1:2
    for j = 1:2
        subplot(2,2,n);
        plot(datafi.fi{ii}{j}.Imean{n},datafi.fi{ii}{j}.FR{n},'.-','Markersize',20,'linewidth',2,'color',col(j,:,:)); hold on
        ylabel('Firing rate, Hz'); xlabel('I, nA'); title(['GM',num2str(n)])
        ylim([0 40]); xlim([0 10]);  xticks([0:2:10]); set(gca, 'Fontsize',16)
        %legend('T=10 C','T=15 C','T=20 C','location','best')
        legend('T=10 C','T=20 C','location','best')        
   % plot traces
    i=5; 
    x=1/Fs:1/Fs:length(datafi.fi{ii}{j}.Vseg{n}{i})/Fs;
    subplot(2,4,4+k); plot(x,datafi.fi{ii}{j}.Vseg{n}{i},'linewidth',1,'color',col(j,:,:))
    x=1/Fs:1/Fs:length(datafi.fi{ii}{j}.Vseg{n}{14})/Fs;
    hold on, plot(x,datafi.fi{ii}{j}.Vseg{n}{14},'linewidth',1,'color',col(j,:,:))
    if ~isempty(datafi.fi{ii}{j}.Vseg_neg{n})
        x=1/Fs:1/Fs:length(datafi.fi{ii}{j}.Vseg_neg{n}{3})/Fs;
        hold on, plot(x,datafi.fi{ii}{j}.Vseg_neg{n}{3},'linewidth',1,'color',col(j,:,:))
    end
    ylim([-86 5]); xlim([0 7]);
    ylabel('Vm, mV'); xlabel('Time, sec'); xticks([0:1:5]); set(gca,'Fontsize',16);
    k=k+1;
   end
end


