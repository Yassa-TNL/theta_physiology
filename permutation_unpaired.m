function p = permutation_unpaired(adata, bdata, reps)

numa = length(adata);
numb = length(bdata);
numtot = numa+numb;

if (size(adata,1)>size(adata,2))
    data = vertcat(adata, bdata);
else
    data = horzcat(adata, bdata);
end

real_meandiff = nanmean(adata(:))-nanmean(bdata(:));

permdiff = NaN*ones(1,reps);
for r = 1:reps
    inds = randperm(numtot);
    indsa = inds(1:numa);
    indsb = inds(numa+1:end);
    permdiff(r) = nanmean(data(indsa))-nanmean(data(indsb));
end


if (real_meandiff>0)
    p = sum(real_meandiff>[real_meandiff permdiff])/(reps+1);
else
    p = sum(real_meandiff>=[real_meandiff permdiff])/(reps+1);
end

p = min(p,1-p)*2;
end

