function LinDisc(PCi, PCj, PCI, PCJ, grouping, PCAWT, PCATG, PCAscore)
    ClassName = Define_Class(grouping);
    PCA1 = PCAscore(:,PCi);
    PCA2 = PCAscore(:,PCj);
    scores = [PCA1 PCA2];
    Model = fitcdiscr(scores,ClassName);
    K = Model.Coeffs(1,2).Const;
    L = Model.Coeffs(1,2).Linear;
    f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
    figure('Name',['PCA scores for PC',PCI,' and PC',PCJ,'with linear discrimination']);
    h1 = scatter(PCAWT(:,PCi),PCAWT(:,PCj),'b');
    hold on
    h2 = scatter(PCATG(:,PCi),PCATG(:,PCj),'r');
    legend([h1, h2], 'WT', 'TG');
    h3 = fimplicit(f,[-2 2 -2 2]);
    hold off
end
function ClassName = Define_Class(grouping)
    for i = 1:length(grouping)
        if grouping(i) == 0
            ClassName(i) = "TG";
        elseif grouping(i) == 1
            ClassName(i) = "WT";
        end
    end
    ClassName = ClassName';
end

