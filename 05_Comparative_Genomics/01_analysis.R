# Inputs are the three final master files and nuclear tree in Analysis_inputs/.

options(stringsAsFactors=FALSE,tick.functions.only=TRUE)
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
if(.Platform$OS.type=='windows')invisible(try(Sys.setlocale('LC_CTYPE','English_United States.utf8'),silent=TRUE))
.source_files <- vapply(sys.frames(), function(fr) if(is.null(fr$ofile)) '' else as.character(fr$ofile)[1], character(1))
.script_args <- grep('^--file=', commandArgs(FALSE), value=TRUE)
.entry_file <- if(any(nzchar(.source_files))) tail(.source_files[nzchar(.source_files)],1) else if(length(.script_args)) sub('^--file=','',.script_args[1]) else ''
ROOT <- if(nzchar(.entry_file)) dirname(dirname(normalizePath(.entry_file,winslash='/',mustWork=TRUE))) else normalizePath(getwd(),winslash='/')
source(file.path(ROOT,'scripts','lib','project.R'),encoding='UTF-8')
dir.create(file.path(ROOT,'outputs'),showWarnings=FALSE,recursive=TRUE)
TASK_ROOT<-ROOT

source(file.path(TASK_ROOT,'scripts','lib','model_functions_B.R'),encoding='UTF-8')
CFG$master_dir<-file.path(TASK_ROOT,'data')
CFG$tree_file<-file.path(CFG$master_dir,'nuclear_tree.nwk')
CFG$output_root<-file.path(TASK_ROOT,'outputs','analysis')
dir.create(CFG$output_root,showWarnings=FALSE,recursive=TRUE)
set.seed(CFG$seed)
stopifnot(requireNamespace('ape',quietly=TRUE),requireNamespace('glmmTMB',quietly=TRUE))
if(packageVersion('glmmTMB')<'1.1.14')stop('glmmTMB >=1.1.14 is required for the phylogenetic propto model.')
args<-commandArgs(TRUE)
unknown<-args[!args%in%'--plan-only' & !grepl('^--workers=[0-9]+$',args)]
if(length(unknown))stop('Unknown argument: ',paste(unknown,collapse=', '))
if(sum(grepl('^--workers=',args))>1L)stop('Specify --workers once.')
PLAN_ONLY<-'--plan-only'%in%args
inp<-load_inputs();Cfull<-inp$Cfull
soft_species<-rownames(Cfull)[sub('_.*','',rownames(Cfull))%in%c('alveonasus','argas','ornithodoros')]
hard_species<-setdiff(rownames(Cfull),soft_species)
stopifnot(length(soft_species)==4L,length(hard_species)==17L)
specs<-scenario_plan();jobs<-list();base_jobs<-list()

add_job<-function(z,s,gene,analysis_class,scheme='',taxon='',skip='') {
  id<-sprintf('M%05d',length(jobs)+1L)
  job<-list(id=id,z=z,s=s,gene=gene,analysis_class=analysis_class,scheme=scheme,taxon=taxon,skip=skip)
  jobs[[length(jobs)+1L]]<<-job
  job
}
for(i in seq_len(nrow(specs))) {
  sp<-specs[i,];d<-switch(sp$domain,Geography=inp$geo,REI=inp$rei,Bacteria=inp$path,Virus=inp$path)
  for(gene in genes_for(sp$domain)) {
    reason<-eligibility(d,sp,gene,Cfull);z<-d[!nzchar(reason),,drop=FALSE];z<-z[order(z$Species),,drop=FALSE]
    job<-add_job(z,sp,gene,if(sp$role=='primary')'Primary'else'Measurement sensitivity')
    if(sp$role=='primary')base_jobs[[length(base_jobs)+1L]]<-job
  }
}
stopifnot(length(jobs)==306L,length(base_jobs)==87L)
for(base in base_jobs) {
  z<-base$z;s<-base$s;gene<-base$gene;s$role<-'secondary';s$high_busco<-FALSE;s$adjust_busco<-FALSE
  genus<-sub('_.*','',z$Species)
  for(taxon in z$Species) {
    ss<-s;ss$scenario<-'LOSO';add_job(z[z$Species!=taxon,,drop=FALSE],ss,gene,'Taxonomic sensitivity','LOSO',taxon)
  }
  for(taxon in sort(unique(genus))) {
    ss<-s;ss$scenario<-'LOGO';add_job(z[genus!=taxon,,drop=FALSE],ss,gene,'Taxonomic sensitivity','LOGO',taxon)
  }
  ss<-s;ss$scenario<-'HARD_ONLY';add_job(z[z$Species%in%hard_species,,drop=FALSE],ss,gene,'Taxonomic sensitivity','HARD_ONLY','Ixodidae')
  ss<-s;ss$scenario<-'NO_IXODES';add_job(z[genus!='ixodes',,drop=FALSE],ss,gene,'Taxonomic sensitivity','NO_IXODES','ixodes')
  for(taxon in sort(soft_species)) {
    ss<-s;ss$scenario<-'ONE_SOFT'
    zz<-z[z$Species%in%c(hard_species,taxon),,drop=FALSE]
    add_job(zz,ss,gene,'Taxonomic sensitivity','ONE_SOFT',taxon,
      if(!taxon%in%z$Species)'Retained soft-tick species absent from the model complete cases'else'')
  }
  if(s$domain%in%c('Bacteria','Virus')) {
    ss<-s;ss$scenario<-'NO_PUBLICATIONS';ss$effort<-'None'
    add_job(z,ss,gene,'Publication sensitivity','NO_PUBLICATIONS','Same species as the publication-adjusted primary model')
  }
}
plan<-bind(lapply(jobs,function(j)data.frame(Job_id=j$id,Domain=j$s$domain,Gene=j$gene,
  Role=j$s$role,Analysis_class=j$analysis_class,Scenario=j$s$scenario,Scheme=j$scheme,Taxon=j$taxon,
  N=nrow(j$z),Gene_category=gene_category(j$s$domain,j$gene),Engine=j$s$engine,
  Adjust_BUSCO=j$s$adjust_busco,High_BUSCO_only=j$s$high_busco,Adjust_studies=j$s$adjust_studies,Effort=j$s$effort,
  Included_species=paste(j$z$Species,collapse=';'),Predefined_unavailable=j$skip)))
write_table(plan,CFG$output_root,'00_frozen_analysis_plan.csv')
write_table(gene_class,CFG$output_root,'00_gene_class.csv')
files<-file.path(CFG$master_dir,c('01_REI_MASTER.csv','02_PATHOGEN_MASTER.csv','03_DISTANCE_MASTER.csv','nuclear_tree.nwk'))
write_table(data.frame(File=basename(files),MD5=unname(tools::md5sum(files))),CFG$output_root,'00_input_checksums.csv')
sig_file<-tempfile(pattern='tick_B_inputs_')
writeLines(unname(tools::md5sum(c(files,file.path(TASK_ROOT,'scripts','01_run_B_reanalysis.R'),file.path(TASK_ROOT,'scripts','lib','model_functions_B.R')))),sig_file)
cache_signature<-substr(unname(tools::md5sum(sig_file)),1L,16L);unlink(sig_file)
writeLines(c('Primary hypothesis family: 87 gene/outcome tests.',
 paste0('Secondary hypothesis family: ',sum(plan$Role=='secondary'),' planned tests; BH pooled across all secondary scenarios.'),
 'Measurement/completeness sensitivities require protein BUSCO >90%.',
 'Taxonomic and publication-adjustment sensitivities use the primary complete-case sample to isolate the specified change.',
 'Failed and unavailable fits retain their slots in their hypothesis family.'),file.path(CFG$output_root,'00_analysis_design.txt'))
cat('PLAN:',nrow(plan),'total;',sum(plan$Role=='primary'),'primary;',sum(plan$Role=='secondary'),'secondary\n')
print(with(plan,table(Analysis_class,Scenario)));flush.console()

fit_job<-function(job) {
  t0<-proc.time()[3];z<-job$z;s<-job$s;gene<-job$gene
  a<-if(nzchar(job$skip))list(ok=FALSE,message=job$skip)else fit_one(z,s,gene,Cfull)
  rec<-result_row(a,z,s,gene)
  rec$Job_id<-job$id;rec$Analysis_class<-job$analysis_class;rec$Scheme<-job$scheme;rec$Taxon<-job$taxon
  co<-dg<-NULL
  if(a$ok&&s$role=='primary') {
    f<-a$f
    co<-data.frame(Job_id=job$id,Domain=s$domain,Scenario=s$scenario,Gene=gene,Term=names(f$beta),
      Estimate=unname(f$beta),SE=unname(f$se),CI_lower=unname(f$lo),CI_upper=unname(f$hi),Raw_p=unname(f$p))
    dg<-tryCatch(diagnose_model(a,s),error=function(e)data.frame(Species=z$Species,Observed=z[[s$response]],
      Fitted=NA_real_,Residual=NA_real_,Residual_definition=conditionMessage(e),Residual_flag=NA,
      Conditional_leverage=NA_real_,Leverage_flag=NA,Descriptive_Pearson_dispersion=NA_real_))
    dg<-cbind(Job_id=job$id,Domain=s$domain,Scenario=s$scenario,Gene=gene,Model_status=rec$Status,dg)
  }
  rec$Elapsed_seconds<-unname(proc.time()[3]-t0)
  list(rec=rec,coef=co,diagnostics=dg)
}
if(PLAN_ONLY){
  cat('PLAN ONLY: validated ',nrow(plan),' fits; no models were fitted.\n',sep='')
  quit(save='no')
}
cores<-parallel::detectCores(logical=FALSE)
if(!is.finite(cores))cores<-1L
workers<-min(4L,max(1L,cores-1L))
worker_arg<-grep('^--workers=',args,value=TRUE)
if(length(worker_arg))workers<-as.integer(sub('^--workers=','',worker_arg))
if(!is.finite(workers)||workers<1L)stop('--workers must be a positive integer.')
cat('Workers:',workers,'; glmmTMB',as.character(packageVersion('glmmTMB')),'\n');flush.console()
worker_log<-tempfile('tick_B_workers_',fileext='.log')
cat('Worker log:',worker_log,'\n')
cl<-parallel::makePSOCKcluster(workers,outfile=worker_log)
parallel::clusterEvalQ(cl,{invisible(try(Sys.setlocale('LC_CTYPE','English_United States.utf8'),silent=TRUE));Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1');options(stringsAsFactors=FALSE);NULL})
parallel::clusterExport(cl,c('make_design','profile_pgls','fit_count','fit_one','result_row','diagnose_model',
 'gene_category','fit_job','CFG','VERSION','gene_class','Cfull'),envir=.GlobalEnv)
results<-vector('list',length(jobs));batch_size<-48L
tryCatch({
  for(first in seq.int(1L,length(jobs),by=batch_size)) {
    last<-min(length(jobs),first+batch_size-1L);cache<-file.path(CFG$output_root,paste0('checkpoint_',cache_signature,sprintf('_%05d_%05d.rds',first,last)))
    rr<-if(file.exists(cache))readRDS(cache)else parallel::parLapplyLB(cl,jobs[first:last],fit_job,chunk.size=1L)
    if(!file.exists(cache))saveRDS(rr,cache)
    results[first:last]<-rr
    stat<-table(vapply(rr,function(x)x$rec$Status,character(1)))
    cat(sprintf('%s completed %d/%d (%s)\n',format(Sys.time(),'%H:%M:%S'),last,length(jobs),paste(names(stat),stat,collapse=', ')));flush.console()
  }
},finally=parallel::stopCluster(cl))
models<-bind(lapply(results,`[[`,'rec'));co<-bind(lapply(results,`[[`,'coef'));dg<-bind(lapply(results,`[[`,'diagnostics'))
models$BH_q<-NA_real_;models$FDR_n<-NA_integer_
for(role in c('primary','secondary')) {
  ii<-which(models$Role==role);models$BH_q[ii]<-adjust_planned(models$P_for_FDR[ii]);models$FDR_n[ii]<-length(ii)
}
write_table(models,CFG$output_root,'01_all_models.csv')
write_table(co,CFG$output_root,'02_primary_coefficients.csv')
write_table(dg,CFG$output_root,'03_primary_species_diagnostics.csv')
tax<-models[models$Analysis_class=='Taxonomic sensitivity',]
full<-models[models$Role=='primary',c('Domain','Gene','Beta','SE','Status','BH_q')]
names(full)[3:6]<-c('Full_beta','Full_SE','Full_status','Full_BH_q')
tax<-merge(tax,full,by=c('Domain','Gene'),sort=FALSE)
tax$Delta_beta<-tax$Beta-tax$Full_beta;tax$Delta_over_full_SE<-tax$Delta_beta/tax$Full_SE
tax$Direction_flip<-ifelse(tax$Status=='OK'&tax$Full_status=='OK',sign(tax$Beta)!=sign(tax$Full_beta),NA)
write_table(tax,CFG$output_root,'04_taxonomic_models.csv')
pub<-models[models$Analysis_class=='Publication sensitivity'|models$Role=='primary'&models$Domain%in%c('Bacteria','Virus'),]
write_table(pub,CFG$output_root,'05_publication_comparison.csv')
if(file.exists(file.path(TASK_ROOT,'reference','01_all_models.csv'))) {
  frozen<-read.csv(file.path(TASK_ROOT,'reference','01_all_models.csv'),check.names=FALSE,fileEncoding='UTF-8-BOM')
  cmp<-merge(models[models$Role=='primary',],frozen[frozen$Role=='primary',],by=c('Domain','Scenario','Gene'),suffixes=c('_new','_previous'))
  out<-cmp[,c('Domain','Gene','N_new','N_previous','Status_new','Status_previous','Beta_new','Beta_previous','BH_q_new','BH_q_previous')]
  write_table(out,CFG$output_root,'06_primary_repeatability.csv')
}
capture.output(sessionInfo(),file=file.path(CFG$output_root,'sessionInfo.txt'))
writeLines('COMPLETE',file.path(CFG$output_root,'RUN_STATUS.txt'))
cat('COMPLETE. Primary discoveries:',sum(models$Role=='primary'&models$BH_q<.05,na.rm=TRUE),'\n')
print(with(models,table(Analysis_class,Status)));flush.console()
