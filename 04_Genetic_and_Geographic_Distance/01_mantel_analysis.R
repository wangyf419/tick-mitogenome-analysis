# Fig. 4-related species-level Mantel/partial Mantel reanalysis.
# This mode always computes new permutation P-values; no archived P-values replace them.
options(stringsAsFactors=FALSE)
if(.Platform$OS.type=='windows')invisible(try(Sys.setlocale('LC_CTYPE','English_United States.utf8'),silent=TRUE))
invisible(Sys.setlocale('LC_COLLATE','C'))
.src<-vapply(sys.frames(),function(fr)if(is.null(fr$ofile))''else as.character(fr$ofile)[1],character(1))
.arg<-grep('^--file=',commandArgs(FALSE),value=TRUE)
.file<-if(any(nzchar(.src)))tail(.src[nzchar(.src)],1)else if(length(.arg))sub('^--file=','',.arg[1])else''
MODULE<-if(nzchar(.file))dirname(dirname(normalizePath(.file,winslash='/',mustWork=TRUE)))else normalizePath(getwd(),winslash='/')
source(file.path(MODULE,'scripts','lib','release_utils.R'),encoding='UTF-8')
source(file.path(MODULE,'scripts','lib','mantel_functions.R'),encoding='UTF-8')
require_packages(c('vegan','geosphere'))
args<-commandArgs(TRUE)
if(any(!args%in%'--check-reference'))stop('Usage: Rscript scripts/01_mantel_analysis.R [--check-reference]')
site<-read.csv(file.path(MODULE,'data/raw/site.csv'),check.names=FALSE)
required<-c('DNA','Longitude','Latitude','Species','Genus')
if(!all(required%in%names(site)))stop('site.csv requires: ',paste(required,collapse=', '))
site<-site[complete.cases(site[,required]),]
site$DNA<-trimws(site$DNA)
if(anyDuplicated(site$DNA))stop('Sequence IDs in site.csv must be unique.')
if(any(!is.finite(site$Longitude))||any(!is.finite(site$Latitude)))stop('Coordinates must be numeric and finite.')
k80<-read_k80(input_matrix_file(MODULE))
rownames(k80)<-trimws(rownames(k80));colnames(k80)<-trimws(colnames(k80))
common_ids<-intersect(site$DNA,intersect(rownames(k80),colnames(k80)))
site<-site[site$DNA%in%common_ids,]
k80<-k80[common_ids,common_ids,drop=FALSE]
if(length(common_ids)<2L)stop('Fewer than two matched sequences.')
diag_mean<-mean(diag(k80),na.rm=TRUE);min_val<-min(k80,na.rm=TRUE);max_val<-max(k80,na.rm=TRUE)
scale_used<-'distance'
if(diag_mean>.95&&diag_mean<=1.0001&&min_val>=0&&max_val<=1.0001){
  k80<-1-k80;scale_used<-'similarity_0_1 converted to 1-similarity'
}else if(diag_mean>95&&diag_mean<=100.1&&min_val>=0&&max_val<=100.1){
  k80<-1-k80/100;scale_used<-'similarity_0_100 converted to 1-similarity/100'
}
k80<-(k80+t(k80))/2;diag(k80)<-0
matrix_dir<-file.path(MODULE,'outputs','analysis','genus_matrices')
dir.create(matrix_dir,showWarnings=FALSE,recursive=TRUE)
RNGkind('Mersenne-Twister','Inversion','Rejection')
set.seed(20260325)
all_genera<-sort(unique(site$Genus))
results<-do.call(rbind,lapply(all_genera,analyze_one_genus,site_all=site,
  seq_dist_mat=k80,nperm=9999,outdir=matrix_dir,min_n_per_species=2))
write.csv(results,file.path(dirname(matrix_dir),'mantel_results_by_genus.csv'),row.names=FALSE,quote=FALSE)
provenance<-data.frame(P_value_source='seeded_reanalysis',Seed=20260325,
  Permutations_requested=9999,Matched_sequences=length(common_ids),Matrix_interpretation=scale_used,
  RNG_kind=paste(RNGkind(),collapse=';'))
write.csv(provenance,file.path(dirname(matrix_dir),'analysis_provenance.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(dirname(matrix_dir),'sessionInfo.txt'))
if('--check-reference'%in%args){
  reference<-read.csv(file.path(MODULE,'data/reference/mantel_results_by_genus.csv'),check.names=FALSE)
  stopifnot(identical(results$Genus,reference$Genus))
  for(name in c('N_species','N_sequences','Min_samples_per_species','Mantel_r','PartialMantel_r'))
    stopifnot(isTRUE(all.equal(results[[name]],reference[[name]],tolerance=1e-12,check.attributes=FALSE)))
  for(genus in results$Genus[is.finite(results$Mantel_r)])for(suffix in c('species_genetic_distance.csv','species_geographic_distance_km.csv','species_latitude_difference.csv')){
    name<-paste0(genus,'_',suffix)
    a<-read.csv(file.path(matrix_dir,name),row.names=1,check.names=FALSE)
    b<-read.csv(file.path(MODULE,'data/reference/genus_matrices',name),row.names=1,check.names=FALSE)
    stopifnot(identical(dim(a),dim(b)),identical(rownames(a),rownames(b)),identical(colnames(a),colnames(b)))
    stopifnot(isTRUE(all.equal(as.matrix(a),as.matrix(b),tolerance=1e-12,check.attributes=FALSE)))
  }
  comparison<-data.frame(Genus=results$Genus,Archived_Mantel_p=reference$Mantel_p,
    New_Mantel_p=results$Mantel_p,Archived_partial_p=reference$PartialMantel_p,New_partial_p=results$PartialMantel_p)
  write.csv(comparison,file.path(dirname(matrix_dir),'reference_p_comparison.csv'),row.names=FALSE)
  cat('PASS: sample counts, correlations and species distance matrices match the archived analysis.\n')
}
cat('Recomputed Mantel results written to outputs/analysis/.\n')
