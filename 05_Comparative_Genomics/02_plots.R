# Plot regenerated estimates, retaining the original figure numbering and candidates.
options(stringsAsFactors=FALSE,tick.functions.only=TRUE)
if(.Platform$OS.type=='windows')invisible(try(Sys.setlocale('LC_CTYPE','English_United States.utf8'),silent=TRUE))
# Resolve this repository from either Rscript or source(), without a personal path.
.source_files <- vapply(sys.frames(), function(fr) if(is.null(fr$ofile)) '' else as.character(fr$ofile)[1], character(1))
.script_args <- grep('^--file=', commandArgs(FALSE), value=TRUE)
.entry_file <- if(any(nzchar(.source_files))) tail(.source_files[nzchar(.source_files)],1) else if(length(.script_args)) sub('^--file=','',.script_args[1]) else ''
ROOT <- if(nzchar(.entry_file)) dirname(dirname(normalizePath(.entry_file,winslash='/',mustWork=TRUE))) else normalizePath(getwd(),winslash='/')
source(file.path(ROOT,'scripts','lib','project.R'),encoding='UTF-8')
dir.create(file.path(ROOT,'outputs'),showWarnings=FALSE,recursive=TRUE)

args<-commandArgs(TRUE)
unknown<-args[!args%in%c('--diagnostics-only','--main-scatter-only') & !grepl('^--source=(reference|analysis)$',args)]
if(length(unknown))stop('Unknown argument: ',paste(unknown,collapse=', '))
source_arg<-grep('^--source=',args,value=TRUE)
if(length(source_arg)>1L)stop('Specify --source once.')
RESULT_SOURCE<-if(length(source_arg))sub('^--source=','',source_arg)else'reference'
RESULT_DIR<-file.path(ROOT,if(RESULT_SOURCE=='reference')'reference'else file.path('outputs','analysis'))
required_results<-c('01_all_models.csv','02_primary_coefficients.csv','03_primary_species_diagnostics.csv','04_taxonomic_models.csv')
if(!all(file.exists(file.path(RESULT_DIR,required_results))))stop('Missing fitted results in ',RESULT_DIR,'. Run 01_run_B_reanalysis.R before --source=analysis.')
cat('Plotting result source:',RESULT_SOURCE,'\n')
PLOT_FONT<-Sys.getenv('TICK_PLOT_FONT',if(.Platform$OS.type=='windows')'Arial'else'sans')
if(!isTRUE(capabilities('cairo')))stop('This plotting script requires an R build with Cairo graphics support.')

ONLY_DIAGNOSTICS<-'--diagnostics-only'%in%commandArgs(TRUE)
MAIN_SCATTER_ONLY<-'--main-scatter-only'%in%commandArgs(TRUE)
source(file.path(ROOT,'scripts','lib','model_functions_B.R'),encoding='UTF-8')
CFG$master_dir<-file.path(ROOT,'data');CFG$tree_file<-file.path(CFG$master_dir,'nuclear_tree.nwk')
inp<-load_inputs();Cfull<-inp$Cfull
read<-function(n)read.csv(file.path(RESULT_DIR,n),check.names=FALSE,fileEncoding='UTF-8-BOM')
models<-read('01_all_models.csv');primary<-models[models$Role=='primary',];coef<-read('02_primary_coefficients.csv')
tax<-read('04_taxonomic_models.csv');dg<-read('03_primary_species_diagnostics.csv')
dir.create(file.path(ROOT,'outputs','figures'),showWarnings=FALSE)
dir.create(file.path(ROOT,'outputs','audit'),showWarnings=FALSE)
stopifnot(requireNamespace('ggplot2',quietly=TRUE),requireNamespace('ggrepel',quietly=TRUE),requireNamespace('patchwork',quietly=TRUE))
genus_colors<-c(Alveonasus='#D77E52',Amblyomma='#69A6A4',Argas='#C88C8F',Dermacentor='#8F69A9',
 Haemaphysalis='#AD7D87',Hyalomma='#C5AB54',Ixodes='#5F8DA8',Ornithodoros='#6C925A',Rhipicephalus='#B3548C')
display<-function(x){a<-c(His_binding='Histamine-binding',Desat='DESAT',cathepsin_D='Cathepsin D',cathepsin_C='Cathepsin C',
 cathepsin_L='Cathepsin L',Cathepsin_B='Cathepsin B',Cytochrome_P450_monooxygenase='Cytochrome P450',Chitin_synthase='Chitin synthase',
 SCP.TAPS.CAP_family_protein='SCP/TAPS/CAP',SCP_TAPS_CAP_family_protein='SCP/TAPS/CAP',
 'SCP/TAPS/CAP_family_protein'='SCP/TAPS/CAP',Tudor_SN='Tudor SN',MYD88='MyD88');v<-unname(a[x]);v[is.na(v)]<-gsub('_',' ',x[is.na(v)]);v}
fmt<-function(x)ifelse(is.finite(x),ifelse(x<.001,formatC(x,format='e',digits=1),sprintf('%.3f',x)),'NE')
sp_label<-function(x){ss<-strsplit(x,'_',fixed=TRUE);vapply(ss,function(z){g<-switch(z[1],alveonasus='Al.',amblyomma='Am.',argas='Ar.',haemaphysalis='Hae.',hyalomma='Hy.',ornithodoros='Or.',rhipicephalus='R.',dermacentor='D.',ixodes='I.');paste(g,paste(z[-1],collapse=' '))},character(1))}
openpdf<-function(name,w,h,ps=9)cairo_pdf(file.path(ROOT,'outputs','figures',name),width=w,height=h,family=PLOT_FONT,pointsize=ps)
audits<-list();label_records<-list()
get_scatter<-function(domain,gene){
 r<-primary[primary$Domain==domain&primary$Gene==gene,];stopifnot(nrow(r)==1,r$Status=='OK')
 data<-if(domain=='REI')inp$rei else inp$geo;sp<-strsplit(r$Included_species,';',fixed=TRUE)[[1]];data<-data[match(sp,data$Species),]
 X<-cbind(Intercept=1,Gene=log10(data[[gene]]+1),Total=log10(data$Genome_Gene+1))
 b<-coef[coef$Domain==domain&coef$Gene==gene,];beta<-setNames(b$Estimate,b$Term)[colnames(X)]
 y<-log10(data[[r$Response]]+1);C<-Cfull[sp,sp,drop=FALSE];V<-r$Lambda*C;diag(V)<-diag(C);Q<-solve(V)
 Ai<-solve(crossprod(X,Q%*%X));e<-drop(y-X%*%beta);s2<-drop(crossprod(e,Q%*%e))/(nrow(X)-ncol(X));vb<-s2*Ai
 delta<-max(abs(drop(Ai%*%crossprod(X,Q%*%y))-beta));sed<-abs(sqrt(vb[2,2])-r$SE);stopifnot(delta<1e-6,sed<1e-6)
 xx<-seq(min(X[,'Gene']),max(X[,'Gene']),length.out=120);newX<-cbind(1,xx,median(X[,'Total']));mu<-drop(newX%*%beta)
 ci<-qt(.975,nrow(X)-ncol(X))*sqrt(rowSums((newX%*%vb)*newX));gg<-sub('_.*','',data$Species);gg<-paste0(toupper(substr(gg,1,1)),substring(gg,2))
 audits[[length(audits)+1L]]<<-data.frame(Domain=domain,Gene=gene,Beta_difference=delta,SE_difference=sed)
 list(r=r,points=data.frame(x=X[,'Gene'],y=y,Species=data$Species,Label=sp_label(data$Species),Genus=gg),
  curve=data.frame(x=xx,fit=mu,lower=mu-ci,upper=mu+ci))
}
scatter<-function(domain,gene,letter=NULL,small=FALSE){
 z<-get_scatter(domain,gene);rr<-range(c(z$points$y,z$curve$lower,z$curve$upper));dy<-diff(rr)
 if(small){
  stopifnot(nrow(z$points)==z$r$N,all(nzchar(z$points$Label)),!anyDuplicated(z$points$Species))
  label_records[[length(label_records)+1L]]<<-cbind(Domain=domain,Gene=gene,z$points)
 }
 yr<-rr+if(small)c(-.10,.10)*dy else c(-.06,.29)*dy
 xr<-range(z$points$x);dx<-diff(xr);if(dx==0)dx<-1
 xl<-if(small)xr+c(-.12,.12)*dx else c(NA_real_,NA_real_)
 lab<-sprintf('β = %.3f; n = %d\nFDR (BH) = %s',z$r$Beta,z$r$N,fmt(z$r$BH_q))
 p<-ggplot2::ggplot()+ggplot2::geom_ribbon(data=z$curve,ggplot2::aes(x=x,ymin=lower,ymax=upper),fill='#D6D6D6',alpha=.55)+
  ggplot2::geom_line(data=z$curve,ggplot2::aes(x=x,y=fit),linetype='dashed',linewidth=.45,color='#555555')+
  ggplot2::geom_point(data=z$points,ggplot2::aes(x=x,y=y,color=Genus),size=if(small).8 else 1.15)+
  ggrepel::geom_text_repel(data=z$points,ggplot2::aes(x=x,y=y,label=Label,color=Genus),family=PLOT_FONT,fontface='italic',
    size=if(small)1.9 else 2.0,seed=20261105,max.overlaps=Inf,max.time=if(small)Inf else 3,max.iter=if(small)200000 else 50000,
    box.padding=if(small).25 else .18,point.padding=if(small).05 else .07,
    min.segment.length=0,segment.color=if(small)'#B8B8B8'else'#C9C9C9',segment.size=.12,
    force=if(small)20 else 3,force_pull=if(small)0 else .2,verbose=small,
    xlim=xl,ylim=if(small)yr else c(rr[1]-.04*dy,rr[2]+.12*dy))+
  (if(small)NULL else ggplot2::annotate('text',x=xr[1]-.02*dx,y=yr[2]-.02*dy,label=lab,hjust=0,vjust=1,size=2.55,family=PLOT_FONT))+
  ggplot2::scale_color_manual(values=genus_colors)+ggplot2::scale_x_continuous(expand=ggplot2::expansion(mult=c(.08,.08)))+
  ggplot2::coord_cartesian(xlim=if(small)xl else NULL,ylim=yr,expand=!small)+ggplot2::labs(x=paste0('log10(',display(gene),' + 1)'),
    y=if(domain=='REI')'log10(median REI + 1)'else'log10(D95 + 1)',title=if(is.null(letter))NULL else paste0('(',letter,') ',display(gene)),
    subtitle=if(small)lab else NULL)+
  ggplot2::theme_classic(base_size=if(small)7 else 9,base_family=PLOT_FONT)+
  ggplot2::theme(legend.position='none',plot.title=ggplot2::element_text(face='bold',size=9),
    plot.subtitle=ggplot2::element_text(size=6.1,lineheight=.95,margin=ggplot2::margin(b=3)),
    axis.title=ggplot2::element_text(size=if(small)6.7 else 8.5),axis.text=ggplot2::element_text(size=if(small)5.8 else 7.2,color='black'),
    axis.line=ggplot2::element_line(linewidth=.35),axis.ticks=ggplot2::element_line(linewidth=.3),plot.margin=ggplot2::margin(4,5,4,3))
 p
}
scatter_grid<-function(file,domain,genes,ncol,w,h,small=FALSE){
 ps<-lapply(seq_along(genes),function(i)scatter(domain,genes[i],if(small)NULL else LETTERS[i],small))
 p<-patchwork::wrap_plots(ps,ncol=ncol)
 ggplot2::ggsave(file.path(ROOT,'outputs','figures',file),p,width=w,height=h,device=cairo_pdf,bg='white',limitsize=FALSE)
}
# Place every species label in measured physical space. Complete labels take
# precedence over proximity; thin leader lines retain the point association.
place_all_species<-function(px,py,width,height,W,H,seed){
 n<-length(px);set.seed(seed);gap<-.45
 candidates<-lapply(seq_len(n),function(i){
  grid<-expand.grid(x=seq(width[i]/2+.7,W-width[i]/2-.7,length.out=40),
    y=seq(height[i]/2+.7,H-height[i]/2-.7,length.out=34))
  dx<-abs(outer(grid$x,px,'-'));dy<-abs(outer(grid$y,py,'-'))
  pc<-rowSums(dx<width[i]/2+1.5 & dy<height[i]/2+1.5)
  grid$base<-pc*1e7+(grid$x-px[i])^2+(grid$y-py[i])^2;grid
 })
 near<-as.matrix(dist(cbind(px,py)));diag(near)<-Inf
 priority<-rowSums(near<22)+width/100
 choose<-function(i,at){
  c<-candidates[[i]];other<-setdiff(which(is.finite(at[,1])),i);score<-c$base
  if(length(other)){
   ox<-abs(outer(c$x,at[other,1],'-'))<matrix((width[i]+width[other])/2+gap,nrow(c),length(other),byrow=TRUE)
   oy<-abs(outer(c$y,at[other,2],'-'))<matrix((height[i]+height[other])/2+gap,nrow(c),length(other),byrow=TRUE)
   score<-score+rowSums(ox&oy)*1e7
  }
  j<-which.min(score)
  if(is.finite(at[i,1])){
   pc<-sum(abs(at[i,1]-px)<width[i]/2+1.5 & abs(at[i,2]-py)<height[i]/2+1.5)
   lc<-if(length(other))sum(abs(at[i,1]-at[other,1])<(width[i]+width[other])/2+gap & abs(at[i,2]-at[other,2])<(height[i]+height[other])/2+gap)else 0
   current<-(pc+lc)*1e7+(at[i,1]-px[i])^2+(at[i,2]-py[i])^2
   if(current<=score[j])return(at[i,])
  }
  c(c$x[j],c$y[j])
 }
 loss<-function(at){
  ox<-abs(outer(at[,1],at[,1],'-'))<outer(width,width,'+')/2+gap
  oy<-abs(outer(at[,2],at[,2],'-'))<outer(height,height,'+')/2+gap
  pcx<-abs(outer(at[,1],px,'-'))<matrix(width/2+1.5,n,n)
  pcy<-abs(outer(at[,2],py,'-'))<matrix(height/2+1.5,n,n)
  outside<-sum(at[,1]-width/2<.65 | at[,1]+width/2>W-.65 | at[,2]-height/2<.65 | at[,2]+height/2>H-.65)
  overlaps<-sum(ox&oy&upper.tri(ox))+sum(pcx&pcy)+outside
  list(count=overlaps,value=overlaps*1e10+sum((at[,1]-px)^2+(at[,2]-py)^2))
 }
 best<-NULL;best_loss<-Inf
 for(run in seq_len(8L)){
  order0<-order(priority+if(run==1L)0 else runif(n,-2,2),decreasing=TRUE)
  at<-matrix(NA_real_,n,2)
  for(i in order0)at[i,]<-choose(i,at)
  for(pass in seq_len(5L))for(i in if(pass%%2)order0 else rev(order0))at[i,]<-choose(i,at)
  value<-loss(at)
  if(value$value<best_loss){best<-at;best_loss<-value$value}
 }
 for(pass in seq_len(4L)){
  improved<-FALSE
  for(i in seq_len(n-1L))for(j in (i+1L):n){
   trial<-best;trial[c(i,j),]<-trial[c(j,i),];value<-loss(trial)
   if(value$value+1e-8<best_loss){best<-trial;best_loss<-value$value;improved<-TRUE}
  }
  for(i in seq_len(n))best[i,]<-choose(i,best)
  best_loss<-loss(best)$value
  if(!improved)break
 }
 stopifnot(loss(best)$count==0L)
 best
}
main_scatter_grid<-function(file,domain,genes){
 openpdf(file,7.1,2.05,8)
 par(mfrow=c(1,3),mar=c(2.9,3.2,2.4,.7),mgp=c(1.9,.45,0),tcl=-.2,
   cex=1,cex.axis=.72,cex.lab=.82,family=PLOT_FONT)
 for(k in seq_along(genes)){
  z<-get_scatter(domain,genes[k]);q<-z$points
  stopifnot(nrow(q)==z$r$N,all(nzchar(q$Label)),!anyDuplicated(q$Species))
  xr<-range(q$x);dx<-diff(xr);rr<-range(c(q$y,z$curve$lower,z$curve$upper));dy<-diff(rr)
  plot(q$x,q$y,type='n',xlim=xr+c(-.12,.12)*dx,ylim=rr+c(-.10,.10)*dy,xaxs='i',yaxs='i',bty='l',
    xlab=paste0('log10(',display(genes[k]),' + 1)'),ylab=if(domain=='REI')'log10(median REI + 1)'else'log10(D95 + 1)')
  polygon(c(z$curve$x,rev(z$curve$x)),c(z$curve$lower,rev(z$curve$upper)),col=adjustcolor('#D6D6D6',alpha.f=.55),border=NA)
  lines(z$curve$x,z$curve$fit,lty=2,lwd=1.1,col='#555555')
  mtext(sprintf('β = %.3f; n = %d',z$r$Beta,z$r$N),side=3,line=1.0,adj=0,cex=.76)
  mtext(paste0('FDR (BH) = ',fmt(z$r$BH_q)),side=3,line=.12,adj=0,cex=.76)
  u<-par('usr');pin<-par('pin')*72;xs<-pin[1]/diff(u[1:2]);ys<-pin[2]/diff(u[3:4])
  px<-(q$x-u[1])*xs;py<-(q$y-u[3])*ys;ce<-.68
  ww<-strwidth(q$Label,units='inches',cex=ce,font=3)*72+2
  hh<-rep(max(strheight('Ag',units='inches',cex=ce,font=3)*72,8*ce)+1.5,nrow(q))
  at<-place_all_species(px,py,ww,hh,pin[1],pin[2],20261105L+k+if(domain=='REI')10L else 0L)
  dxp<-px-at[,1];dyp<-py-at[,2];t<-pmax(abs(dxp)/(ww/2),abs(dyp)/(hh/2))
  ex<-at[,1]+dxp/t;ey<-at[,2]+dyp/t
  segments(q$x,q$y,u[1]+ex/xs,u[3]+ey/ys,col='#B0B0B0',lwd=.35)
  points(q$x,q$y,pch=16,cex=.40,col=genus_colors[q$Genus])
  text(u[1]+at[,1]/xs,u[3]+at[,2]/ys,q$Label,cex=ce,font=3,col=genus_colors[q$Genus])
  label_records[[length(label_records)+1L]]<<-cbind(Domain=domain,Gene=genes[k],q,
    Label_x_pt=at[,1],Label_y_pt=at[,2],Label_width_pt=ww,Label_height_pt=hh,
    Point_x_pt=px,Point_y_pt=py,Panel_width_pt=pin[1],Panel_height_pt=pin[2])
 }
 dev.off()
}
# C and D retain the original Figure 5 families and their original ordering.
if(!ONLY_DIAGNOSTICS){
main_scatter_grid('Figure5_C.pdf','Geography',c('His_binding','Cystatin','Kunitz'))
main_scatter_grid('Figure5_D.pdf','REI',c('ACC','Desat','Cathepsin_B'))
write.csv(do.call(rbind,label_records),file.path(ROOT,'outputs','audit','Figure5_all_species_labels.csv'),row.names=FALSE)
if(MAIN_SCATTER_ONLY){
 capture.output(sessionInfo(),file=file.path(ROOT,'outputs','audit','plot_sessionInfo.txt'))
 a<-do.call(rbind,audits);path<-file.path(ROOT,'outputs','audit','prediction_audit.csv')
 if(file.exists(path)){old<-read.csv(path);a<-rbind(old[!paste(old$Domain,old$Gene)%in%paste(a$Domain,a$Gene),],a)}
 write.csv(a,path,row.names=FALSE)
 write.csv(do.call(rbind,label_records),file.path(ROOT,'outputs','audit','Figure5_all_species_labels.csv'),row.names=FALSE)
 quit(save='no')
}
scatter_grid('FigureS38.pdf','Geography',c('Cytochrome_P450_monooxygenase','His_binding','Chitin_synthase','ACC'),2,7.1,5.8)
scatter_grid('FigureS40.pdf','REI',c('FAR','FAS','Vg','cathepsin_C','cathepsin_D','cathepsin_L','LAP'),2,7.1,8.3)

# Split bacterial/viral forest panels retain the original E-panel structure.
openpdf('Figure5_E.pdf',1.85,2.05,7.5);layout(matrix(1:2,2,1),heights=c(3.6,1.7))
for(domain in c('Bacteria','Virus')){
 z<-primary[primary$Domain==domain&primary$Status=='OK'&is.finite(primary$BH_q)&primary$BH_q<.05,];yy<-rev(seq_len(nrow(z)))
 par(mai=if(domain=='Bacteria')c(.02,.43,.18,.37)else c(.30,.43,.08,.37),mgp=c(1.3,.35,0),tcl=-.16)
 plot(z$Beta,yy,xlim=c(-8.5,.5),ylim=c(.5,nrow(z)+.5),yaxt='n',ylab='',xlab=if(domain=='Virus')'Gene coefficient (β)'else'',
   xaxt=if(domain=='Bacteria')'n'else's',pch=16,cex=.65,cex.axis=.7,cex.lab=.65,bty='o',col=if(domain=='Bacteria')'#D74D43'else'#559FC4')
 segments(z$CI_lower,yy,z$CI_upper,yy,col=if(domain=='Bacteria')'#D74D43'else'#559FC4');abline(v=0,lty=2,col='grey65')
 axis(2,at=yy,labels=fmt(z$BH_q),las=2,tick=FALSE,cex.axis=.62);axis(4,at=yy,labels=display(z$Gene),las=2,tick=FALSE,cex.axis=.65)
 if(domain=='Bacteria')mtext('FDR (BH)',side=3,line=.25,adj=0,cex=.65)
}
dev.off()

forest<-function(file,domain,h){
 z<-primary[primary$Domain==domain,];z<-z[rev(seq_len(nrow(z))),];yy<-seq_len(nrow(z));openpdf(file,7.1,h,10)
 par(mar=c(4,11,1.5,7),mgp=c(2.4,.6,0),tcl=-.22);cols<-ifelse(z$BH_q<.05,'#C34E43','#737373')
 plot(z$Beta,yy,xlim=range(c(z$CI_lower,z$CI_upper,0)),yaxt='n',ylab='',xlab='Gene coefficient (β) and 95% CI',pch=16,col=cols,bty='l')
 segments(z$CI_lower,yy,z$CI_upper,yy,col=cols,lwd=1.2);abline(v=0,lty=2,col='grey65')
 axis(2,at=yy,labels=display(z$Gene),las=2,tick=FALSE,cex.axis=.8)
 axis(4,at=yy,labels=paste0('P = ',fmt(z$Raw_p),'   FDR = ',fmt(z$BH_q)),las=2,tick=FALSE,cex.axis=.63)
 dev.off()
}
forest('FigureS39.pdf','Geography',6.2);forest('FigureS41.pdf','REI',5.6)

# Publication intensity: same four original gene/outcome examples and same-species comparisons.
openpdf('FigureS43.pdf',8.5,9,10);layout(matrix(c(1,2,3,4),2,2,byrow=TRUE),heights=c(1,1.05))
pubcor<-list()
for(domain in c('Bacteria','Virus')){
 z<-inp$path[is.finite(inp$path[[domain]])&is.finite(inp$path$PubMed_total_count),];x<-log10(z$PubMed_total_count+1);y<-log10(z[[domain]]+1)
 cor<-cor.test(z$PubMed_total_count,z[[domain]],method='spearman',exact=FALSE)
 pubcor[[domain]]<-data.frame(Response=domain,N=nrow(z),Rho=unname(cor$estimate),P=cor$p.value)
 par(mar=c(4.2,4.3,2.8,1),mgp=c(2.6,.6,0),tcl=-.2)
 plot(x,y,pch=16,col=if(domain=='Bacteria')'#B0555E'else'#5D86B1',bty='l',xlab='log10(total publications + 1)',
   ylab=paste0('log10(',if(domain=='Bacteria')'bacterial'else'viral',' breadth + 1)'),main=paste0('(',if(domain=='Bacteria')'A'else'B',')'))
 g<-seq(min(x),max(x),length.out=100);pr<-predict(lm(y~x),newdata=data.frame(x=g),interval='confidence')
 polygon(c(g,rev(g)),c(pr[,2],rev(pr[,3])),border=NA,col=adjustcolor('#BBC7D5',alpha.f=.4));lines(g,pr[,1],col='#63758E')
 points(x,y,pch=16,col=if(domain=='Bacteria')'#B0555E'else'#5D86B1');legend('topleft',legend=paste0('ρ = ',sprintf('%.2f',cor$estimate),'\nP = ',fmt(cor$p.value),'\nn = ',nrow(z)),bty='n',cex=.8)
}
for(domain in c('Bacteria','Virus')){
 genes<-if(domain=='Bacteria')c('MYD88','Pelle')else c('MYD88','Tudor_SN')
 a<-models[models$Analysis_class=='Publication sensitivity'&models$Domain==domain,];a<-a[match(genes,a$Gene),]
 b<-primary[primary$Domain==domain,];b<-b[match(genes,b$Gene),];yy<-rev(seq_along(genes))
 lim<-range(c(0,a$CI_lower[a$Status=='OK'],a$CI_upper[a$Status=='OK'],b$CI_lower[b$Status=='OK'],b$CI_upper[b$Status=='OK']))
 par(mar=c(5.2,5.5,2.8,1));plot(NA,xlim=lim,ylim=c(.2,length(genes)+.8),yaxt='n',ylab='',xlab='NB2 coefficient (β) and 95% CI',bty='l',main=paste0('(',if(domain=='Bacteria')'C'else'D',')'))
 for(k in 1:2){z<-if(k==1)a else b;ok<-z$Status=='OK';y<-yy+if(k==1).12 else -.12;cc<-if(k==1)'#315C89'else'#B56D32';segments(z$CI_lower[ok],y[ok],z$CI_upper[ok],y[ok],col=cc);points(z$Beta[ok],y[ok],pch=if(k==1)16 else 17,col=cc)}
 axis(2,at=yy,labels=display(genes),las=2,tick=FALSE);abline(v=0,lty=2,col='grey65')
 legend('bottomleft',legend=c('Without publication covariate','With publication covariate'),pch=c(16,17),col=c('#315C89','#B56D32'),bty='n',cex=.65)
}
dev.off();write.csv(do.call(rbind,pubcor),file.path(ROOT,'outputs','audit','publication_correlations.csv'),row.names=FALSE)

# Recreate the original A/B/C taxonomic forest layout, now with all five schemes.
sets<-list(A=data.frame(Domain='Geography',Gene=c('His_binding','Kunitz','Cystatin','Chitin_synthase')),
 B=data.frame(Domain='REI',Gene=c('Desat','cathepsin_D','ACC','cathepsin_C','FAR','FAS','LAP','cathepsin_L','Cathepsin_B','Vg','SCP/TAPS/CAP_family_protein')),
 C=rbind(data.frame(Domain='Bacteria',Gene=c('MYD88','Pelle','BSK','STAM','TrCP','p38')),data.frame(Domain='Virus',Gene=c('Tudor_SN','MYD88','TOLLIP'))))
taxpanel<-function(which){
 target<-sets[[which]];key<-paste(target$Domain,target$Gene);q<-primary[paste(primary$Domain,primary$Gene)%in%key,];v<-tax[paste(tax$Domain,tax$Gene)%in%key&tax$Status=='OK',]
 lim<-range(c(0,q$Beta[q$Status=='OK'],v$Beta));pad<-.06*diff(lim);if(pad==0)pad<-.1;lim<-lim+c(-pad,pad)
 par(mar=if(which=='B')c(4,6.4,3,1.1)else c(4,6.8,3,1.1),mgp=c(2.4,.55,0),tcl=-.2)
 yy<-rev(seq_len(nrow(target)));plot(NA,xlim=lim,ylim=c(.35,nrow(target)+.65),yaxt='n',ylab='',xlab='Gene coefficient (β)',bty='l',
   main=paste0('(',which,') ',switch(which,A='D95 geographic range',B='Median REI',C='Reported pathogen breadth')),cex.main=.95,cex.axis=.8,cex.lab=.85)
 abline(v=0,lty=2,col='grey70')
 for(i in seq_len(nrow(target))){
  r<-primary[primary$Domain==target$Domain[i]&primary$Gene==target$Gene[i],];z<-tax[tax$Domain==target$Domain[i]&tax$Gene==target$Gene[i]&tax$Status=='OK',];y<-yy[i]
  if(r$Status=='OK')points(r$Beta,y+.24,pch=if(r$BH_q<.05)18 else 5,col='#214F7C',cex=.85)
  for(scheme in c('HARD_ONLY','NO_IXODES')){a<-z[z$Scheme==scheme,];if(nrow(a))points(a$Beta,y+if(scheme=='HARD_ONLY').1 else -.02,pch=if(scheme=='HARD_ONLY')16 else 15,col=if(scheme=='HARD_ONLY')'#E28629'else'#3C9863',cex=.72)}
  for(scheme in c('ONE_SOFT','LOGO','LOSO')){
   a<-z[z$Scheme==scheme,];if(!nrow(a))next;off<-switch(scheme,ONE_SOFT=-.15,LOGO=-.28,LOSO=-.4);cc<-switch(scheme,ONE_SOFT='#C9503E',LOGO='#B0B0B0',LOSO='#6E6E6E')
   segments(min(a$Beta),y+off,max(a$Beta),y+off,col=cc,lty=if(scheme=='LOSO')2 else 1,lwd=1)
   if(scheme=='ONE_SOFT')points(median(a$Beta),y+off,pch=21,bg='white',col=cc,cex=.8)
  }
 }
 labs<-display(target$Gene);if(which=='C')labs<-paste0(ifelse(target$Domain=='Bacteria','Bacteria: ','Virus: '),labs)
 axis(2,at=yy,labels=labs,las=2,tick=FALSE,cex.axis=if(which=='C').67 else .76)
}
openpdf('FigureS44.pdf',8.5,10.3,10);layout(matrix(c(1,2,3,2,4,4),3,2,byrow=TRUE),widths=c(1,1),heights=c(.34,.51,.12))
taxpanel('A');taxpanel('B');taxpanel('C');par(mar=c(0,0,0,0));plot.new()
legend('center',legend=c('Full model','Hard ticks only','Excluding Ixodes','One soft tick: median/range','LOGO range','LOSO range'),
 pch=c(18,16,15,21,NA,NA),lty=c(NA,NA,NA,1,1,2),pt.bg=c(NA,NA,NA,'white',NA,NA),
 col=c('#214F7C','#E28629','#3C9863','#C9503E','#B0B0B0','#6E6E6E'),ncol=3,bty='n',cex=.75)
dev.off()

# Preserve the original six diagnostics and Q-Q layout. Count models use
}
# conditional randomized quantile residuals rather than Gaussian Pearson Q-Q.
targets<-data.frame(Domain=c('Geography','Geography','REI','REI','Bacteria','Virus'),Gene=c('His_binding','Chitin_synthase','Desat','ACC','MYD88','Tudor_SN'))
openpdf('FigureS45.pdf',8.5,11,12);par(mfrow=c(6,2),cex=1,mar=c(3.4,4.0,2.6,.9),mgp=c(2.0,.55,0),tcl=-.2,cex.axis=.9,cex.lab=.9)
qr_out<-list()
for(i in seq_len(nrow(targets))){
 domain<-targets$Domain[i];gene<-targets$Gene[i];r<-primary[primary$Domain==domain&primary$Gene==gene,];z<-dg[dg$Domain==domain&dg$Gene==gene,]
 stopifnot(r$Status=='OK');flag<-abs(z$Residual)>3|z$Conditional_leverage>=.8
 plot(if(r$Engine=='PGLS')log10(z$Fitted+1)else z$Fitted,z$Residual,pch=1,cex=.7,col=ifelse(flag,'#B94C3D','black'),bty='l',
  xlab=if(r$Engine=='PGLS')'Fitted log10 phenotype'else'Fitted pathogen count',ylab=if(r$Engine=='PGLS')'Studentized residual'else'Pearson residual',
  main='')
 mtext(paste0('(',LETTERS[i],') ',if(domain=='Geography')'D95'else domain,' — ',display(gene)),side=3,line=1.25,cex=.9,font=2)
 abline(h=0,lty=2,col='grey50');mtext(paste0('n = ',r$N,if(r$Engine=='PGLS')paste0('; λ = ',sprintf('%.3f',r$Lambda))else'; NB2'),side=3,line=.15,cex=.65)
 if(r$Engine=='PGLS'){
  data<-if(domain=='REI')inp$rei else inp$geo;data<-data[match(z$Species,data$Species),];X<-cbind(Intercept=1,Gene=log10(data[[gene]]+1),Total=log10(data$Genome_Gene+1))
  b<-coef[coef$Domain==domain&coef$Gene==gene,];beta<-setNames(b$Estimate,b$Term)[colnames(X)];e<-log10(z$Observed+1)-drop(X%*%beta)
  C<-Cfull[z$Species,z$Species,drop=FALSE];V<-r$Lambda*C;diag(V)<-diag(C);w<-drop(backsolve(chol(V),e,transpose=TRUE));w<-w/sqrt(sum(w^2)/(nrow(X)-ncol(X)))
  kind<-'Standardized GLS residual'
 }else{
  set.seed(20261105L+i);low<-pnbinom(z$Observed-1,size=r$NB2_dispersion,mu=z$Fitted);high<-pnbinom(z$Observed,size=r$NB2_dispersion,mu=z$Fitted)
  u<-runif(nrow(z),low,high);w<-qnorm(pmax(1e-12,pmin(1-1e-12,u)));kind<-'Conditional randomized quantile residual'
 }
 qqnorm(w,pch=1,cex=.7,main='',xlab='Theoretical quantiles',ylab='Sample quantiles',bty='l');mtext('Normal Q–Q',side=3,line=1.25,cex=.9,font=2);qqline(w,col='grey40')
 qr_out[[i]]<-data.frame(Domain=domain,Gene=gene,Species=z$Species,QQ_residual=w,Definition=kind)
}
dev.off()
write.csv(do.call(rbind,qr_out),file.path(ROOT,'outputs','audit','QQ_diagnostics.csv'),row.names=FALSE)
if(length(audits))write.csv(do.call(rbind,audits),file.path(ROOT,'outputs','audit','prediction_audit.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(ROOT,'outputs','audit','plot_sessionInfo.txt'))
cat('Publication plots complete; C/D and S38/S40/S43/S45 retain original candidate selections.\n')
