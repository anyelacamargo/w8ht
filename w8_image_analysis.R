
# This script prepares data for QTL mapping
rm(list=ls()); # Delete files
cat("\014");
setwd("M:/anyela/repo/w8ht");
library(nlme)
#ibrary('mpMap')
#library(rrBLUP)

source('../senescence_disease/generic.R')

#Find parents in pedigree table
processPedrigree = function(cname, pedrigreetable)
{
  l = list()
  cname = as.character(cname);
  
  s = strsplit(cname, '-')[[1]];
  lname = paste(s[1], substr(s[2],1,1), sep='-');
  i = which(lname == pedrigreetable$Gen.0.8.way.F1);
  if(length(i) > 0)
  {
    mother = as.character(pedrigreetable$X4.way.F[i]);
    father = as.character(pedrigreetable$X4.way.M[i]);
    return(paste(mother, father, sep='_'));
  }
  else
  {
    return(cname);
  }
}



processPheno = function(pheno, meta)
{
  
  phenotable = merge(meta[, c('barcode', 'genotype', 'idtype')], pheno, 
                     by.y='Plant', by.x='barcode');
  return(phenotable);
  
}

plotMean = function()
{
  copydata = o;
  p <- ggplot(data = copydata, aes_string(x = 'Group.1', y = 'Area.in.square.mm')) +
    labs(x = "DAS", y = 'ba') + geom_smooth(aes(group = 1), span = 0.3) + 
      stat_summary(fun.y=mean, geom="line", size=0.5) + 
    theme_bw(base_size = 8) + 
    theme(panel.background = element_blank(),
          axis.text.x = element_text(size=8,angle = 90, hjust = 1 )) ;
  print(p)
  
}


#cname= 'Area.in.square.mm'; v1 = d$traitdata$DAS; v2=d$traitdata[[cname]];traitdata = d$traitdata; plabel = c('DAS', cname, 'genotype');
plotAll = function(v1, v2, data, plabel, lgroup)
{
  copydata = data;
  p = qplot(v1, v2, data=copydata, colour= lgroup$genotypes,
            xlab=plabel[1], ylab=plabel[2], ymax=max(v2)) +
    theme_bw() + 
    geom_line(position=position_dodge(width=0)) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.text.x=element_text(angle=90, size = 14)) +
          scale_colour_manual(values= lgroup$cd, guide = guide_legend(title = plabel[3])) +
          scale_x_continuous(breaks = unique(v1)[seq(1,83, by=5)], labels = unique(v1)[seq(1,83, by=5)]);
  
  p
  return(p);
  
}

# Plot average per group
# data = dataset
# traitpheno = array with orig train name, short name, units
# traittime = time variable
# 
plotLongAverage = function(data, traitpheno, traittime, fillname, lgroup)
{
 
  copydata = data;
  p = ggplot(data = copydata,
              aes_string(y = traitpheno, x = traittime, group = 'genogroup')) +
    scale_colour_manual(values = unique(lgroup$cd)) +
    labs(x = traittime, y = traitpheno) +
    theme(legend.key = element_rect(colour = 'transparent', fill = 'transparent')) +
    #stat_summary(fun.y=mean, geom="line", size=0.8) + 
    geom_smooth(aes(colour=genogroup), method='loess', size = 1.2) + 
    theme_bw(base_size = 15) + 
    theme(panel.background = element_blank(),
          axis.text.x = element_text(size=15)); # Face bold
  return(p);
  
}

plotGR = function(data, traitpheno, traittime, fillname, lgroup)
{
  
  copydata = data;
  p = ggplot(data = copydata,
             aes_string(y = traitpheno, x = traittime, group = 'genogroup')) +
    scale_colour_manual(values = unique(lgroup$cd)) +
    labs(x = traittime, y = traitpheno) +
    theme(legend.key = element_rect(colour = 'transparent', fill = 'transparent')) +
    #stat_summary(fun.y=mean, geom="line", size=0.8) + 
    #geom_smooth(aes(colour=genogroup), method='loess', size = 1.2) + 
    geom_line() +
    theme_bw(base_size = 15) + 
    theme(panel.background = element_blank(),
          axis.text.x = element_text(size=15)); # Face bold
  p
  return(p);
  
}

WUECalculate = function(data, riltable, dp)
{
  phenodata = data;
  wue_total = averageValues(phenodata[, c('water_amount')], list(genotype=phenodata$genotype, 
                            genogroup=phenodata$genogroup), 'sum');
  g = phenodata[which(phenodata$DAS == dp), ];
  wue_total = merge(wue_total, g[c('genotype', 'Area.in.square.mm', 'Plant.height.in.mm')], 
                    by.x='genotype', by.y='genotype');
  wue_total = data.frame(wue_total, WUE=sapply(as.numeric(row.names(wue_total)), 
                        function(x) ((wue_total$Area.in.square.mm[x])/1000)/wue_total$x[x]));
  phenotable = merge(riltable, wue_total, by.x='genotype', by.y='genotype');
  phenotable = phenotable[order(phenotable$WUE),];
  phenotable <- resetlevels(phenotable, 'genotype');
  return(phenotable);
}




selecDataPoint = function(data, ril, dplist, trait, traitred)
{
  copydata = data;
  i = which(copydata$DAS == dplist[1]);
  wholedata = copydata[i, which(colnames(copydata) == c('genotype', trait))];
  colnames(wholedata)[2] = paste(traitred, dplist[1], sep='_');
  for(dp in dplist[2:length(dplist)])
  {
    i = which(copydata$DAS == dp);
    subdata = copydata[i, which(colnames(copydata) == c('genotype', trait))];
    colnames(subdata)[2] = paste(traitred, dp, sep='_');
    wholedata = merge(wholedata, subdata, by.x = 'genotype', by.y = 'genotype');
  }
  wholedata = merge(d$riltable, wholedata, by.x='genotype', 'genotype');
  colnames(wholedata)[2] = c('SUBJECT.NAME')
 
    
  return(wholedata)
}


# Find time points with complete set of lines
findCompleteTimePoint  = function(data)
{
  copydata = data;
  mx = max(table(copydata$DAS))
  m = data.frame(which(table(copydata$DAS) == mx));
  return(rownames(m));
}

createGroups = function(data)
{
  k = list();
  copydata = data;
  parent = which(copydata$genogroup == 'RIL');
  k$genotypes = rep('RIL', nrow(copydata));
  k$genotypes[parent] = as.character(copydata$genotype[parent])
  
  
  k$cd = c('RIL' = "pink", 'Alchemy-1A' = "red", 'Avalon' = "gray", 'Brompton-1A' = 'green4', 'Cadenza' = 'black', 
         'Claire-4' = 'yellow', 'Hereward-1A' = 'purple', 'KWS Santiago' = 'brown', 'Rialto-4' = 'deepskyblue1', 
         'Robigus-1' = 'coral1', 'Soissons-1' = 'green', 'Xi-19-1' = 'magenta1', 'Zircon' = 'blue');
  
  return(k);
}

organiseWater = function(data)
{
  watertable = data;
  watertable = process_timestamp(watertable, 'time_stamp',' ', 'timewatered');
  watertable = watertable[, c(-(which(colnames(watertable) == 'time_stamp')))];
  watertable$water_amount[watertable$water_amount == 0] = 0.01;
  watertable = changeDateFormat(watertable, 'date');
  return(watertable);
  
}

plotData = function(data, lgroup, tnames, fillname=NULL)
{
  traitdata = data;
  for(cname in tnames)
  {
    tiff(paste(cname, '.tiff', sep=''),  width = 2080, height = 1080,res=200);
    print(plotLongAverage(traitdata, cname, 'DAS', fillname, lgroup));
    dev.off();
  }
}


plotWUETtotal = function(data)
{
  wuetable = data;
  p <- ggplot(data = wuetable, mapping = aes(x = genotype, y = WUE )) +
    layer(geom = "point") + theme(panel.grid.major = element_blank(),
                                  panel.grid.minor = element_blank(),
                                  panel.background = element_rect(colour = "black", fill = "white"),
                                  axis.text.x=element_blank(),
                                  axis.ticks.x=element_blank(),
                                  axis.text.y=element_text(size = 14));
  return(p);
  
}

loadData = function()
{

  d = list();
  d$imagetable = read_data('w8_HT_avg_pheno.csv', TRUE, ',');
  d$metatable = read_data('w8_metadata.csv', TRUE, ',');
  d$phenotable = processPheno(d$imagetable, d$metatable);
  d$riltable = read_data('ril_name.csv', TRUE, ',');
  d$pedigree = read_data('MAGIC_pedigree.csv', TRUE, sep=',');
  d$traittable = read_data('trait_table.csv', TRUE, sep=',');
  
  
  w = gsub("\n", " ", "SELECT snapshot.id_tag, snapshot.time_stamp, weight_before, weight_after,
  snapshot.water_amount
           FROM public.snapshot WHERE snapshot.id_tag LIKE 'W8_%' AND snapshot.water_amount <> '-1'
           ORDER BY snapshot.time_stamp ASC;");
  
  db_list = list('dbdriver' = 'PostgreSQL', 'dbip'='144.124.105.250', 'dbport' = '5432',
                 'dbname' = 'production', 'dbuser' = 'postgres', 'dbpass' = 'LemnaTec');
  d$watertable = connect2DB(w, db_list);
  
  d$watertable = organiseWater(d$watertable);
  d$datamerged = merge(d$phenotable[,-5], d$watertable, by.x=c('barcode', 'Date'), 
                       by.y=c('id_tag', 'date')); # Deleted average col
  d$traitdata = averageValues(d$datamerged[,5:(ncol(d$datamerged)-1)],
                    list(date=d$datamerged$Date, genotype=d$datamerged$genotype, 
                         genogroup=d$datamerged$idtype));
  
  
  d$traitdata = data.frame(d$traitdata, 
             parents = sapply(d$traitdata$genotype, function(x) processPedrigree(x, d$pedigree)));
  
  DAS = createDAS(ds="2014-10-20", de="2015-04-30");
  d$traitdata = merge(DAS, d$traitdata, by.x='date', by.y='date');
  d$traitdata = d$traitdata[order(c(as.Date(d$traitdata$date, format='%d/%m/%Y')), d$traitdata$genogroup),];
  m = findCompleteTimePoint(d$traitdata);
  j = match(d$traitdata$DAS, m);
  d$traitdata = d$traitdata[!is.na(j),];
  #d$traitdata1 <- resetlevels(d$traitdata, 'date');
  
  return(d);
  
}

exportTraitsforQTL = function(data)
{
  for(i in seq(2,nrow(data$traittable), by=2))
  {
    
    f = selecDataPoint(data$traitdata, data$riltable, 
                       unique(data$traitdata$DAS)[seq(1, length(unique(data$traitdata$DAS)), by=4)],
                       as.character(data$traittable$Ltrait[i]), as.character(data$traittable$Name[i]));
    write.table(f[, -1], file=paste('w8_HTImg', data$traittable$Name[i], '.phenotype', sep=''), 
                sep='\t', row.names = F, quote = F);
  }
}


exportTraitsforQTLRange = function(data, grdata, tname)
{
  
  wholedata = merge(data$riltable, grdata, by.x='genotype', 'genotype');
  head(wholedata)
  colnames(wholedata)[2] = c('SUBJECT.NAME');
  write.table(wholedata[, -1], file=paste(tname, '.phenotype', sep=''), 
              sep='\t', row.names = F, quote = F);
  
}


growthModelLogistic = function(data, genotype_list)
{
  library(lme4);
  g = list();
  pdf('logistic_curves.pdf');
  par(mfrow = c(4,2));
  for(ename in genotype_list)
  {
    sub = data[which(data[['genotype']] == ename),];
    gr <- nls(log(Area.in.square.mm) ~ SSlogis(DAS, phi1, phi2, phi3), data = sub);
    alpha <- coef(gr);  #extracting coefficients
    g[[ename]] = alpha
    #plot(log(Area.in.square.mm) ~ DAS, data = sub, main = ename, 
    #   xlab = "DAS", ylab = paste("Area", '(log)', sep =''))  # Census data
    #curve(alpha[1]/(1 + exp(-(x - alpha[2])/alpha[3])), add = T, col = "red", lwd=2)  # Fitted model
    #text(170, min(log(sub$Area.in.square.mm))+1.5, paste('alpha: ', round(alpha[1],2), 
    #                   '\n  xmid:', round(alpha[2],2), '\n',  ' scale:', round(alpha[3],2), sep=''),cex = 0.8);
  
    with(sub,plot(DAS, log(Area.in.square.mm), main = ename, 
                  xlab = "DAS", ylab = paste("Area", '(log)', sep ='')))
    with(sub,lines(DAS, predict(gr),col="red",lty=2,lwd=3));
    text(170, min(log(sub$Area.in.square.mm))+1.5, paste('alpha: ', round(alpha[1],2), 
                       '\n  xmid:', round(alpha[2],2), '\n',  ' scale:', round(alpha[3],2), sep=''),cex = 0.8);
    
  }
  dev.off();
  
  g = t(data.frame(g));
  g = data.frame(genotype=genotype_list,g, row.names=(1:nrow(g)));
  return(g)
  
}



d = loadData();
gr = growthModelLogistic(d$traitdata, unique(d$traitdata$genotype));


break();
lgroup = createGroups(d$traitdata);
# Plot data
plotData(d$traitdata, lgroup, c('Plant.height.in.mm', 'Area.in.square.mm', 'TYP.area', 'TFL.area',
                                'YFL.area', 'water_amount'), NA);

#Calculate WUE Total
wuetable = WUECalculate(d$traitdata, d$riltable, 189);
lgroup = createGroups(wuetable);
tiff(paste('wue_total', '.tiff', sep=''),  width = 2080, height = 1080,res=200);
plotWUETtotal(wuetable)
dev.off();


wuetableDaily = WUEDaily(d$traitdata);
lgroup = createGroups(wuetableDaily);
tiff(paste('wueDaily', '.tiff', sep=''),  width = 2080, height = 1080,res=200);
plotData(wuetableDaily, lgroup, 'WUE');
dev.off();

#Export traits for QTL analysis
exportTraitsforQTL(d);
