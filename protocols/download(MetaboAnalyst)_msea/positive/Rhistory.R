# PID of current job: 1686941
mSet<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-c("Ethanolamine phosphate","3-Aminobutyric acid","Pyridoxal","O-Succinylhomoserine","3-Hydroxybutyric acid","Adipic acid","Choline","Carnitine","Suberic acid","N-Acetylglucosamine 6-phosphate","Gly-Gly","Lauric acid","Adenosine","N-methyl-L-glutamic Acid","7,8-Dihydrobiopterin","Thioproline","1-Methyl-4-imidazoleacetic acid","GMP","SDMA","2,6-Diaminopimelic acid","3',5'-ADP","Pyridoxamine 5'-phosphate","Pantothenic acid","Adenine","2-Hydroxyisobutyric acid","Pimelic acid","Guanosine","Ophthalmic acid","2-Hydroxyvaleric acid","2-Deoxy-D-glucose 6-phosphate","Hypoxanthine","3-Hydroxy-3-methylglutaric acid","1-Methyladenosine","Betaine aldehyde","Cholic acid","5-Oxoproline","2-Aminobutyric acid","GDP","Phosphocreatine","N-Acetylglucosamine-1-phosphate","Quinolinic acid","CDP-choline","Cytidine","N-Acetyl-beta-alanine","Asp")
mSet<-Setup.MapData(mSet, cmpd.vec);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-Setup.HMDBReferenceMetabolome(mSet, "all_detected.txt");
mSet<-SetMetabolomeFilter(mSet, T);
mSet<-SetCurrentMsetLib(mSet, "smpdb_pathway", 2);
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_0_", "net", "png", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_0_", "png", 72, width=NA)
mSet<-CalculateHyperScore(mSet)
mSet<-PlotORA(mSet, "ora_1_", "net", "png", 72, width=NA)
mSet<-PlotEnrichDotPlot(mSet, "ora", "ora_dot_1_", "png", 72, width=NA)
