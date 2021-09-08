% load the membrane potential data of GM neurons @ 10 and 20 C

load('C:\Users\moroz\Documents\code\hco_analysis\data\GM_FI_Rin.mat')

%% calculate Rin & FIs

data10=[]; data20=[];
Rin10=[]; Rin20=[];

for j=1:numel(data.V) % experiments
    
    data10{j}=fi(data.V{j}(1,:),data.I{j}(1,:),Fs); % 10oC
    data20{j}=fi(data.V{j}(2,:),data.I{j}(2,:),Fs); % 20oC
    
    Rin10 = [Rin10, data10{j}.Rin1];
    Rin20 = [Rin20, data20{j}.Rin1];
    
end

%% Wilcoxon signed rank-sum test

[p,h] = signrank(Rin10,Rin20)

%% plot input resistance @ 10oC and 20oC

addpath('C:\Users\moroz\Dropbox\ml\plotting')
clf, plot([1 2],[Rin10; Rin20],'.-','markersize',20,'linewidth',1,'color','k'), hold on
xlim([0.9 2.1]); ylim([0 16])
xticks([1 2]); xticklabels({'10^oC','20^oC'});
    ylabel('Input resistance, MOhms'); box off
set(gca,'Fontsize',16,'FontName','Myriad Pro')
display.sigstar([1,2],p)
text(1.3,15.5,['n.s. (p=',num2str(round(p,4)),')'],'fontsize',16,'FontName','Myriad Pro')
text(1.9,15.5,'N=10','fontsize',14,'FontName','Arial')
set(gca,'Fontsize',16,'FontName','Arial')
