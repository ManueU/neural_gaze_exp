function cor1pt = nangauss_smooth(cor1pt,gaussmean,gaussstd)
    gfilt = fspecial('gaussian',gaussmean,gaussstd);
    nanmask = isnan(cor1pt);
    [r,c] = find(~nanmask);
    [rnan,cnan] = find(nanmask);
    F1 = scatteredInterpolant(c,r,cor1pt(~nanmask),'nearest');
    cor1pt(nanmask) = F1(cnan,rnan);
    cor1pt = conv2(cor1pt,gfilt,'same');
    cor1pt(nanmask) = nan;
end