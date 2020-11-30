function [group2wtTG, WT, TG] = sort_spectra(namespectra1, spectra, wavenumbers)
    %% Sorts the normalised_spectra into WT or TG
    group2wtTG = contains(namespectra1, 'wt' );
    groupIndx = repmat(group2wtTG, 1, length(spectra)); %logical index to the same size as the Spectra
    NumWT = numel(group2wtTG(group2wtTG == 1)); % gets how many WT 
    WTb = spectra(groupIndx); 
    WT = reshape(WTb, NumWT, wavenumbers);
    TGb = spectra(groupIndx == 0); 
    NumTG = numel(group2wtTG(group2wtTG == 0));
    TG = reshape(TGb, NumTG, wavenumbers); 
end