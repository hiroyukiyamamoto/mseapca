# PID of current job: 1816139
mSet<-InitDataObjects("conc", "msetora", FALSE)
cmpd.vec<-c("Uridine","Fumaric acid","4-Guanidinobutanoic acid","Thr-Asp","Hexylamine","Sarcosine","2-Methylserine","Ectoine","L-Homoserine","cis-Aconitic acid","Dimethylglycine","Cysteic acid","Urea","Isobutyryl-CoA","L-Acetylcarnitine","Glycolic acid","Isocitric acid","Trehalose 6-phosphate","ADP-ribose","N-Acetylleucine","Ornithine","N-Acetylhistidine","N-Acetylornithine","Dihydroxyacetone phosphate","Malic acid","Glycylleucine","Trimethylamine N-oxide","L-Theanine","Glucosamine","2-Phosphoglyceric acid","Glycerophosphocholine","2,4-Diaminobutyric acid","gamma-Glutamylcysteine","Dyphylline","L-Tryptophan","4-Trimethylammoniobutanoic acid","3-Phosphoglyceric acid","S-Methylmethionine","Phosphorylcholine","Argininosuccinic acid","Gamma-Aminobutyric acid","L-Alanine","Ganciclovir","L-Cysteine","N-Acetyl-L-phenylalanine","Sorbitol 6-phosphate","Xanthopterin","6-Phosphogluconic acid","Trigonelline","Phosphoenolpyruvic acid","L-Arginine","Succinic acid","Homovanillic acid","Spermine","Theobromine","Guanidinosuccinic acid","Spermidine","Glutathione","Pyruvic acid","L-Lactic acid","5-Aminopentanoic acid","Fructose 6-phosphate","Stachydrine","Glucose 6-phosphate","Glucose 1-phosphate","S-Lactoylglutathione","Betonicine")
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
