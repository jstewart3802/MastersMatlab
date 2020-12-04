function Model = LinDisc(PCi, PCj, PCI, PCJ, grouping, PCAWT, PCATG, PCAscore)
    %call classname fucntion
    ClassName = Define_Class(grouping);
    %get the relevant columns from our data for each PC being used
    PCA1 = PCAscore(:,PCi);
    PCA2 = PCAscore(:,PCj);
    %Combine into one matrix
    scores = [PCA1 PCA2];
    %produce model from data
    Model = fitcdiscr(scores,ClassName);
    %Extract coefficients from model to be used in plotting the model
    K = Model.Coeffs(1,2).Const;
    L = Model.Coeffs(1,2).Linear;
    %create eqauatiojn of line to be plotted
    f = @(x1,x2) K + L(1)*x1 + L(2)*x2;
    %create figure of PCi vs PCj with discrimination shown
    figure('Name',['PCA scores for PC',PCI,' and PC',PCJ,'with linear discrimination']);
    h1 = scatter(PCAWT(:,PCi),PCAWT(:,PCj),'b');
    hold on
    h2 = scatter(PCATG(:,PCi),PCATG(:,PCj),'r');
    h3 = fimplicit(f,[-2 2 -2 2]);
    legend([h1, h2, h3], 'WT', 'TG','Linear Discrimination');
    xlabel(['PC ',PCI]);
    ylabel(['PC ',PCJ]);
    hold off
end
%-------------------Subfunctions-------------------------------------------
function ClassName = Define_Class(grouping)
    for i = 1:length(grouping)
        if grouping(i) == 0
            ClassName(i) = "TG";
        elseif grouping(i) == 1
            ClassName(i) = "WT";
        end
    end
    %transpose to match format of fitcdiscr classname shape
    ClassName = ClassName'; 
end

