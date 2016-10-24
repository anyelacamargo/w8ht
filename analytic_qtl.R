
# This script prepares data for QTL mapping
rm(list=ls()); # Delete files
cat("\014");
setwd("M:/anyela/repo/w8ht");
source('../senescence_disease/generic.R')



# Creat QTL table from MAGIC output 
# 
createQTLTable = function(mapdata)
{
  setwd('qtl');
  filelist = stackFiles('qtls');
  filelist = checkFileSize(filelist);
  qtldata = stackDataSets(filelist);
  
  imputedata = stackDataSets(stackFiles('imputed'));
  h = merge(qtldata, mapdata, by.x='peak.SNP', by.y='Marker');
  h = h[order(h$phenotype,h$chromosome, h$cM),];
  h$logP = round(h$logP,2);
  h$genomewide.pvalue = round(h$genomewide.pvalue,2);
  h$cM = round(h$cM,2);
  accession_list = as.vector(unique(imputedata$Accession));
  accession = data.frame(matrix(1, nrow=nrow(h), ncol= length(accession_list)));
  colnames(accession) = accession_list;
  qtlfull = data.frame(h, accession);
  
  for(trait in as.vector(h$phenotype))
  {
    for(marker in as.vector(h$peak.SNP))
    {
      
      for(a in accession_list)
      {
     
        m = intersect(which(imputedata$Phenotype == trait),
                      intersect(which(imputedata$SNP == marker), which(imputedata$Accession == a)));
        
        if(length(m) > 0)
        {
          i = which(qtlfull$peak.SNP == marker);
          j = which(qtlfull$phenotype == trait);
          i = intersect(i,j);
          qtlfull[i, a] = round(imputedata[m,'mean'],2)
        }
        
      }
      
    }
  }
  
return(qtlfull);
}

# Get genetic map
get_MAGIC = function(mapfile)
{
  f = list()
  f$map_raw = read.table(mapfile, sep=',', header=TRUE);
  return(f);
}


plotManhattan = function(traitlist)
{
  filelist = stackFiles('scan');
  for(tname in traitlist)
  {
    pdf(paste('qtl_', tname,'.pdf'));
        #,  width = 1400, height = 1080, #, res=200););
    s = filelist[grep(tname, filelist)];
    s = s[order(nchar(s), s)];
    u = loadMarkerFile(s, TRUE, '\t', 'from');
    chrname = gsub("chr", "", unique(u$chrom));
    
    par(mfrow=c(4,2))
    for(cname in colnames(u)[4:ncol(u)])
    {
      manhattan(u[,c('marker','chrom', 'pos',cname)], fdr.level = 0.05, cname, chrname, 0.8)
    }
    dev.off();
  }
}


plotTraitbyDate = function(data)
{
  copydata = t;
  copydata = data.frame(copydata, DAS = sapply(t$phenotype, function(x) strsplit(as.character(x), '_')[[1]][2]));  
  copydata = data.frame(copydata, trait = sapply(t$phenotype, function(x) strsplit(as.character(x), '_')[[1]][1]));
  y = count(copydata, c('DAS', 'trait', 'Chr'))
  y$DAS = as.numeric(as.character(y$DAS));
  y$freq = factor(y$freq);
  colnames(y)[4] = 'QTLperChr';
  y = y[order(y$DAS),];
  y <- resetlevels(y, 'DAS');
  
  
  p <- ggplot(y, aes(Chr, DAS))+ 
    geom_point(aes(size = QTLperChr, colour = QTLperChr)) + facet_grid(. ~ trait)
  return(p);
}

f = get_MAGIC('../senescence_disease/wheat_geno_coordinates.csv'); 
break();
t = createQTLTable(f$map_raw);

write.table(t[,c(2,1,9,10,7,8,12:19)], file='qtltable.csv', sep=',', row.names = F);

tiff(paste('QTLbyDay', '.tiff', sep=''),  width = 2080, height = 1080,res=200);
plotTraitbyDate(t);
dev.off();

plotManhattan(c('Area', 'Height', 'YFL'))
