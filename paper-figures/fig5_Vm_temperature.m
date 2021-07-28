% load the membrane potential data of GM neurons @ 10 and 20 C

load('C:\Users\moroz\Documents\code\hco_analysis\data\GM_Vm_fig5.mat')

% calculate mean membrane potential
V10mean=[];  V20mean=[]; k=1;
for j = 1:numel(data.V)
    for i=1:2 % neuron number
        try
            V10mean(k)=mean(data.V{j}{i}(:,1));
            V20mean(k)=mean(data.V{j}{i}(:,2));
            k=k+1;
        end
    end
end

% Wilcoxon signed rank-sum test
[p,h] = signrank(V10mean,V20mean)

% plot mean Vm @ T=10 and T=20 
clf, plot([1 2],[V10mean; V20mean],'.-','markersize',25,'linewidth',1.5,'color','k');
hold on, plot([0.9 0.9],[mean(V10mean); mean(V10mean)],'.','markersize',30,'color','b');
hold on, plot([2.1 2.1],[mean(V20mean); mean(V20mean)],'.','markersize',30,'color','r');
hold on, errorbar(0.9,mean(V10mean),std(V10mean),'b')
hold on, errorbar(2.1,mean(V20mean),std(V20mean),'r')
xlim([0.8 2.2]); ylim([-74 -47]) 
xticks([1 2]); xticklabels({'T=10^oC','T=20^oC'});
sigstar([1,2],p)
text(1.6,-48,['(p=',num2str(round(p,6)),')'],'fontsize',16,'FontName','Myriad Pro')
ylabel('GM Membrane potential, mV'); box off 
text(1.9,-51,['N=',num2str(length(V10mean))],'Fontsize',20,'FontName','Arial')
set(gca,'Fontsize',16,'FontName','Arial')


