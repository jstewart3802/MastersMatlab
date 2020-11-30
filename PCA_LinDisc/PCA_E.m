%% Plots the %variance explained by each PC
function PCA_E(PCAvariances, maxPC)
    MAXPC = int2str(maxPC);
    figure('Name',['Percent variance explained by First ', MAXPC ,' PCs']);
    plot(1:maxPC,100*cumsum(PCAvariances(1:maxPC))/sum(PCAvariances(1:maxPC)),'-ob');
    xlabel('Number of PCs')
    ylabel('Percentage of variance explained in x')
    clear MAXPC;
end