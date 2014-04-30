function inx = geneidxfinder(gene,genes)

len = length(gene);
for i = 1:length(genes)
    if genes(i,1:len) == gene
        inx = i;
        break
    end
end
    