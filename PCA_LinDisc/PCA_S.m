%% Plot scores for WT and TG for PCs as set in PCi and PCj
function PCA_S(PCi, PCj, PCI, PCJ, PCAWT, PCATG)
    figure('Name',['PCA scores for PC',PCI,' and PC',PCJ]);
    sc1 = scatter(PCAWT(:,PCi),PCAWT(:,PCj),'b');
    hold on
    sc2 = scatter(PCATG(:,PCi),PCATG(:,PCj),'r');
    legend([sc1(1), sc2(1)], 'WT', 'TG');
    xlabel(['PC ',PCI]);
    ylabel(['PC ',PCJ]);
end