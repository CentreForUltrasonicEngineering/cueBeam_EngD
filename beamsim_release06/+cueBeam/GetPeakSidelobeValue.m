function out=GetPeakSidelobeValue(ValueVector)
% use peakdet to estimate peak sidelobe next to the main lobe
try % this can fail if no peaks or main lobe is the sidemost lobe
    [pmax pmin]=cueBeam.peakdet(ValueVector,0.1);
    % find the idx of the highest peak and treat it as main lobe
    [err idx_main]=max(pmax(:,2));
    % get the higher of the two side lobes
    try % try getting right sidelobe
        rsidelobe=pmax(idx_main+1,2);
    catch E
        rsidelobe=NaN;
    end
    try
        lsidelobe=pmax(idx_main-1,2);
    catch E
        lsidelobe=NaN;
    end
    out=max(rsidelobe,lsidelobe);
catch E
    out=NaN;
end

end
