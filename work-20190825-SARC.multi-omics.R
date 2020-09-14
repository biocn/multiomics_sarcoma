#------------------for PAAD multi-omics:
library(iClusterPlus)
library(GenomicRanges)
library(lattice)
library(RTCGAToolbox)
library(bnstruct)
library(ggpubr)
library(NMF)
library(RColorBrewer)
library(data.table)
library(survival)
library(KMsurv)
library(RColorBrewer)
library(glmnet)
library(survivalROC)
library(survminer)
library(pheatmap)
library(survcomp)
library(data.table)
library(clustree)# for hclust to plot tree
draw_survial_curve_custom<-function(myd,column,bk,myd_colors){
	myd[myd[,column]!="",]->myd.rm;
	myd.rm[!is.na(myd.rm[,column]),]->myd.rm;
	if(length(myd.rm$A1_OS)>1){
		Surv(as.numeric(myd.rm$A1_OS),as.numeric(myd.rm$Status))->myd.surv;
	}else if(length(myd.rm$OS)>1){
		Surv(as.numeric(myd.rm$OS),as.numeric(myd.rm$Status))->myd.surv;
	}
	survfit(formula=myd.surv~myd.rm[,column])->myd.fit;
	survdiff(formula=myd.surv~myd.rm[,column],rho=0)->myd.diff;
	table(myd.rm[,column])->myd.table;
	max(myd.rm$A1_OS)+100->max_xlim;
	plot(myd.fit,col=myd_colors,xlab="Time(days)",ylab="Overall Survival(%)",lwd=2,axes=F,main=paste("KM(",colnames(myd.rm)[column],")",sep=""),xlim=c(-100,max_xlim*1.1));
	axis(side=1,at=seq(0,max_xlim,bk),labels=seq(0,max_xlim,bk),pos=0);
	rug(x=seq(0,max_xlim,bk)+bk/2,ticksize=-0.01,side=1,pos=0);
	axis(side=2,at=seq(0,1,0.2),labels=seq(0,100,20),pos=0);
	rug(x=seq(0,0.9,0.2)+0.1,ticksize=-0.01,side=2,pos=0);
	#abline(h=seq(0.2,1,0.2),col=brewer.pal(9,"Greys")[3],lty=3)
	1-pchisq(myd.diff$chisq,df=length(myd.diff$n)-1)->pvalue;
	legend("topright",legend=paste(names(myd.table),paste("(N=",myd.table,")",sep="")),fill=myd_colors,bty="n",cex=1.2);
	if(pvalue<1e-5){
		legend(x=bk/2,y=0.20,legend="p < 1e-5",bty="n",cex=1.2)
	}else{
		legend(x=bk/2,y=0.20,legend=paste("p=",round(pvalue,5),sep=""),bty="n",cex=1.2)
	}
	return(c(pvalue,myd.table));
}
cox_univariant_gene_regr<-function(myd,colums){
	Surv(as.numeric(myd$A1_OS),as.numeric(myd$Status))->myd.surv
	c()->univar_anova_p;
	c()->univar_coxph_HR;
	c()->univar_coxph_low95;
	c()->univar_coxph_high95;
	c()->univar_coxph_logtest;
	c()->univar_coxph_sctest;
	c()->univar_coxph_waldtest;
	c()->fpkm.mean;
	c()->fpkm.median;
	colnames(myd)[colums]->myd.names;
	for(i in myd.names){
		mean(myd[,i],na.rm=T)->tmp.mean;
		median(myd[,i],na.rm=T)->tmp.median;
		c(fpkm.mean,tmp.mean)->fpkm.mean;
		c(fpkm.median,tmp.median)->fpkm.median;
		as.formula(paste("myd.surv~",i))->tmp.formula;
		coxph(formula=tmp.formula,data=myd)->tmp.coxph;
		summary(tmp.coxph)->tmp.coxph.summary;
		c(univar_anova_p,tmp.coxph.summary$coefficients[,5])->univar_anova_p;
		c(univar_coxph_HR,tmp.coxph.summary$coefficients[,2])->univar_coxph_HR;
		c(univar_coxph_low95,tmp.coxph.summary$conf.int[,3])->univar_coxph_low95;
		c(univar_coxph_high95,tmp.coxph.summary$conf.int[,4])->univar_coxph_high95;
		c(univar_coxph_logtest,tmp.coxph.summary$logtest[3])->univar_coxph_logtest;
		c(univar_coxph_sctest,tmp.coxph.summary$sctest[3])->univar_coxph_sctest;
		c(univar_coxph_waldtest,tmp.coxph.summary$waldtest[3])->univar_coxph_waldtest;	
	}
	data.frame("gName"=myd.names,"Pvalue"=univar_anova_p,"HR"=univar_coxph_HR,"Low(95%CI)"=univar_coxph_low95,"High(95%CI)"=univar_coxph_high95,"Logrank"=univar_coxph_logtest,"Sctest"=univar_coxph_sctest,"Waldtest"=univar_coxph_waldtest,"fpkm_median"=fpkm.median,"fpkm_mean"=fpkm.mean)->myd.coxph.df;
	#myd.coxph.df[myd.coxph.df$fpkm_median>=1,]->myd.coxph.df;
	myd.coxph.df[!is.na(myd.coxph.df$Logrank),]->myd.coxph.df;
	myd.coxph.df[order(myd.coxph.df$Logrank),]->myd.coxph.df;
	return(myd.coxph.df);
}
#---read exp data 
setwd("d:/TCGA-SARC/")
c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(8,"Dark2"),brewer.pal(12,"Set3"),brewer.pal(8,"Accent"),brewer.pal(11,"Spectral")[c(1,4,10,11)])->myd_colors;
15->g_start_column;
#--------------------read CNV data
read.table("CNV/nocnv/Copy Number Variation.nocnv.merge.txt",header=T,sep="\t",stringsAsFactors=F)->myd.cnv
myd.cnv[myd.cnv$Chromosome!="X",]->myd.cnv#***********
as.numeric(myd.cnv$Chromosome)->myd.cnv$Chromosome;#************
CNregions(seg=myd.cnv,epsilon=0,adaptive=FALSE,frac.overlap=0.5,rmSmallseg=TRUE,nProbes=5)->myd.cnv_merged
data.frame("Samples"=rownames(myd.cnv_merged),myd.cnv_merged)->myd.cnv_merged
myd.cnv_merged[grep("-10",myd.cnv_merged$Samples,invert=T),]->myd.cnv_merged_rm
fwrite(myd.cnv_merged_rm,"multi-omics/SARC.CNV_merged.txt",quote=F,sep="\t",row.names=F)
#-
readLines("multi-omics/common_samples")->myd_samples;
read.table("multi-omics/cnv_interval2genes.txt",header=T,sep="\t")->cnv_interval2genes
myd.cnv_merged[myd_samples,]->myd.cnv_merged;
########read Methylation file
fread("multi-omics/SARC_450k_0.7cutoff.rmCrossMap.factors.res",header=T,sep="\t",stringsAsFactors=F)->myd.methy
as.data.frame(myd.methy)->myd.methy
which(gsub("\\.","-",colnames(myd.methy))%in%myd_samples)->methy_index;
myd.methy[,c(1:9,methy_index)]->myd.methy
#knn impute method to handle NA values
t(myd.methy[,10:ncol(myd.methy)])->myd.methy_t;
knn.impute(myd.methy_t,k=10,cat.var=1:ncol(myd.methy_t),to.impute=1:nrow(myd.methy_t),using=1:nrow(myd.methy_t))->myd.methy_t
data.frame(myd.methy[,c(1:9)],t(myd.methy_t))->myd.methy;
write.table(myd.methy,"multi-omics/SARC_450k_07.methy_filter.txt",quote=F,sep="\t",row.names=F)
########read EXP file
preprocess_expd_v2<-function(expd,gStartColumn,aliveEvent,log_trans){
	#-----------remove OS ==NA
	expd[!is.na(expd$A1_OS),]->expd;
	expd[!is.na(expd$A2_Event),]->expd;
	expd[expd$A1_OS!=0,]->expd;
	expd[expd$A1_OS>=30,]->expd;
	#----------remove gene value==1 or 0 in all samples-------
	c()->filter_colums;
	for(i in gStartColumn:ncol(expd)){
		#length(expd[expd[,i]<1,i])->failed_samples;#for TCGA
		length(expd[expd[,i]==0,i])->failed_samples;#for GEO
		if(failed_samples/nrow(expd)<0.5){
			c(filter_colums,i)->filter_colums
		}
	}	
	#----change colnames: '-' to '.'
	gsub("-",".",colnames(expd)[filter_colums])->names(expd)[filter_colums];
	#--------------------------------
	if(log_trans=="yes"){
		data.frame(expd[,1:(gStartColumn-1)],log2(expd[,filter_colums]+1))->expd.filter;
	}else{
		data.frame(expd[,1:(gStartColumn-1)],expd[,filter_colums])->expd.filter;
	}
	print(length(filter_colums));
	flush.console();
	#---------status: 0->alive,1->death---------
	c()->status;
	for(i in 1:nrow(expd.filter)){
		if(expd.filter$A2_Event[i]==aliveEvent){
			c(status,0)->status;
		}else{
			c(status,1)->status; 
		}
	}
	status->expd.filter$Status
	return(expd.filter);
}
preprocess_expd<-function(expd,gStartColumn,aliveEvent){
	#-----------remove OS ==NA
	expd[!is.na(expd$A1_OS),]->expd;
	expd[!is.na(expd$A2_Event),]->expd;
	expd[expd$A1_OS!=0,]->expd;
	expd[expd$A1_OS>=30,]->expd;
	#----------remove IRGP value==1 or 0 in all samples-------
	c()->filter_colums;
	for(i in gStartColumn:ncol(expd)){
		length(expd[expd[,i]<1,i])->failed_samples;
		if(failed_samples/nrow(expd)<0.5){
			c(filter_colums,i)->filter_colums
		}
	}
	#----change colnames: '-' to '.'
	gsub("-",".",colnames(expd)[filter_colums])->names(expd)[filter_colums];
	#--------------------------------
	expd[,c(1:(gStartColumn-1),filter_colums)]->expd.filter;
	print(length(filter_colums));
	flush.console();
	#---------status: 0->alive,1->death---------
	c()->status;
	for(i in 1:nrow(expd.filter)){
		if(expd.filter$A2_Event[i]==aliveEvent){
			c(status,0)->status;
		}else{
			c(status,1)->status; 
		}
	}
	status->expd.filter$Status
	return(expd.filter);
}
read.table("drug-mining/HTSeq-FPKM_SARC.merge.SYMBOL.factors.txt",header=T,sep="\t",stringsAsFactors=F)->myd.exp
myd.exp$A0_Samples->rownames(myd.exp)
myd.exp[myd_samples,]->myd.exp;
#----------low fpkm filter: fpkm<1 samples more than 50% total samples will be removed-------
#----for new events(metastasis),not overall surival;
myd.exp->myd.exp_newEvent;
myd.exp_newEvent$NewEvent->myd.exp_newEvent$A2_Event;
myd.exp_newEvent$NewEventTime->myd.exp_newEvent$A1_OS;
preprocess_expd(myd.exp_newEvent,g_start_column,"0")->myd.exp_newEvent_processed;#should be "0"
preprocess_expd_v2(myd.exp_newEvent,g_start_column,"0","yes")->myd.exp_newEvent_processed_log2;#should be "0"
#-----------------
preprocess_expd(myd.exp,g_start_column,"Alive")->myd_exp_processed;
preprocess_expd_v2(myd.exp,g_start_column,"Alive","yes")->myd_exp_processed_log2;
myd_exp_processed$A0_Samples->myd_samples;
#------------------------------------univariant cox regression: new events not Overall survival time
cox_univariant_gene_regr(myd.exp_newEvent_processed,c(g_start_column:c(ncol(myd.exp_newEvent_processed)-1)))->x#->myd.exp_processed.coxph;
myd.exp_processed.coxph[!is.na(myd.exp_processed.coxph$Pvalue),]->myd.exp_processed.coxph
myd.exp_processed.coxph[myd.exp_processed.coxph$Pvalue<0.05,"gName"]->coxph_filtered_genes;
as.character(coxph_filtered_genes)->coxph_filtered_genes;
#----------do correlation calculate for CNV_EXP
#------------------------------------------------------------------for regions merged to one gene
fread("multi-omics/SARC.CNV_merged_regions2symbol.txt",header=T,sep="\t",stringsAsFactors=F)->myd.cnv_region2symbol;
as.data.frame(myd.cnv_region2symbol)->myd.cnv_region2symbol;
myd.cnv_region2symbol$Sample->rownames(myd.cnv_region2symbol);
myd.cnv_region2symbol[myd_samples,]->myd.cnv_region2symbol;
cal_correlation_CNV_EXP_region2symbol<-function(cnvd,expd,cnvd_start_column){
	#-------calculate correlation
	cnvd$Sample->rownames(cnvd);
	expd$A0_Samples->rownames(expd);
	cnvd[myd_samples,]->cnvd;
	expd[myd_samples,]->expd;
	colnames(cnvd)->cnvd_names;
	c()->cnvd_expd_cor;
	c()->cnvd_expd_cor_pvalue;
	c()->cnvd_genes;
	c()->cnvd_regions;
	c()->expd_genes;
	for(i in cnvd_start_column:ncol(cnvd)){
		cnvd_names[i]->g_name;
		as.character(g_name)->g_name;
		if(length(is.na(expd[1,g_name]))==0){
			next;
		}
		tryCatch({
			cor.test(cnvd[,i],expd[,g_name])->tmp_cor;
			c(cnvd_regions,cnvd_names[i])->cnvd_regions;
			c(cnvd_genes,as.character(g_name))->cnvd_genes;
			c(expd_genes,as.character(g_name))->expd_genes;
			c(cnvd_expd_cor,tmp_cor$estimate)->cnvd_expd_cor;
			c(cnvd_expd_cor_pvalue,tmp_cor$p.value)->cnvd_expd_cor_pvalue;
			},error=function(e){cat("Error:",g_name,"\n");}
		)
	}
	data.frame("CnvRegion"=cnvd_regions,"CnvGene"=cnvd_genes,"ExpGene"=expd_genes,"Cor"=cnvd_expd_cor,"Cor.pvalue"=cnvd_expd_cor_pvalue)->res_df;
	return(res_df);
}
t(myd.cnv_region2symbol[,-1])->myd.cnv_region2symbol_t;
apply(myd.cnv_region2symbol_t,1,function(x){2^x})->test;
data.frame("Samples"=myd.cnv_region2symbol$Sample,test)->myd.cnv_region2symbol_logTrans;
colnames(myd.cnv_region2symbol)->colnames(myd.cnv_region2symbol_logTrans);

cal_correlation_CNV_EXP_region2symbol(myd.cnv_region2symbol_logTrans,myd_exp_processed,2)->CNV_region2symbol_EXP.cor;
unlist(lapply(CNV_region2symbol_EXP.cor[,4],function(x){log((1+x)/(1-x))*0.5}))->CNV_region2symbol_EXP.cor$Z_value;
CNV_region2symbol_EXP.cor[CNV_region2symbol_EXP.cor$Cor.pvalue<1e-5,]->CNV_region2symbol_EXP.cor_filter;
CNV_region2symbol_EXP.cor_filter[CNV_region2symbol_EXP.cor_filter$Z_value>0,]->CNV_region2symbol_EXP.cor_filter;
CNV_region2symbol_EXP.cor_filter[order(CNV_region2symbol_EXP.cor_filter$Z_value,decreasing=T),]->CNV_region2symbol_EXP.cor_filter;
table(CNV_region2symbol_EXP.cor_filter$CnvGene)->CnvGene_table;
names(CnvGene_table[CnvGene_table!=0])->CNVCor_region2symbol_genes;
write.table(CNV_region2symbol_EXP.cor,"multi-omics/CNV_region2symbol_EXP.cor.txt",row.names=F,sep="\t",quote=F)
#-------do correlation calculate for MET_EXP
cal_correlation_MET_EXP<-function(metd,expd,metd_start_column,expd_start_column){
	#--------map gene to probes
	c()->g_probes;
	names(table(metd$Symbol))->col1_table;
	for(s in col1_table){
		which(metd$Symbol%in%s)->tmp_probes;
		c(g_probes,paste(tmp_probes,collapse="_"))->g_probes;
	}
	col1_table->names(g_probes);
	#-------prepare matrix;
	t(metd[,metd_start_column:ncol(metd)])->metd_t;
	gsub("\\.","-",rownames(metd_t),perl=T)->rownames(metd_t);
	metd_t[myd_samples,]->metd_t;
	metd$Probe->colnames(metd_t);
	#print(metd_t[1:10,1:10])
	#flush.console();
	#-------calculate correlation
	colnames(expd)->expd_names;
	c()->metd_expd_cor;
	c()->metd_expd_cor_pvalue;
	c()->metd_probes;
	c()->metd_genes;
	c()->expd_genes;
	for(i in expd_start_column:ncol(expd)){
		expd_names[i]->g_name;
		g_probes[g_name]->g_name_probes;
		if(is.na(g_name_probes)){
			next;
		}
		as.numeric(unlist(strsplit(g_name_probes,"_")))->g_name_probes;
		for(p in g_name_probes){	
			cor.test(expd[,i],metd_t[,p])->tmp_cor;
			c(metd_probes,as.character(colnames(metd_t)[p]))->metd_probes;
			c(metd_genes,g_name)->metd_genes;
			c(expd_genes,g_name)->expd_genes;
			c(metd_expd_cor,tmp_cor$estimate)->metd_expd_cor;
			c(metd_expd_cor_pvalue,tmp_cor$p.value)->metd_expd_cor_pvalue;
		}
	}
	data.frame("MetProbe"=metd_probes,"MetGene"=metd_genes,"ExpGene"=expd_genes,"Cor"=metd_expd_cor,"Cor.pvalue"=metd_expd_cor_pvalue)->res_df;
	return(res_df);
}
cal_correlation_MET_EXP(myd.methy,myd_exp_processed,10,g_start_column)->MET_EXP.cor
unlist(lapply(MET_EXP.cor[,4],function(x){log((1+x)/(1-x))*0.5}))->MET_EXP.cor$Z_value;
MET_EXP.cor[MET_EXP.cor$Cor.pvalue<1e-5,]->MET_EXP.cor_filter
MET_EXP.cor_filter[MET_EXP.cor_filter$Z_value<0,]->MET_EXP.cor_filter
MET_EXP.cor_filter[order(MET_EXP.cor_filter$Z_value),]->MET_EXP.cor_filter;
table(MET_EXP.cor_filter$MetGene)->MetGene_table;
names(MetGene_table[MetGene_table!=0])->METCor_genes;
write.table(MET_EXP.cor,"multi-omics/MET_EXP.cor.txt",row.names=F,sep="\t",quote=F)
#------------------------------------------------------for CNVCor genes:
CNVCor_genes_count<-function(myd,intervals){
	names(table(as.character(myd$CnvGene)))->genes;
	lapply(genes,function(x){which(intervals$Gene%in%x)[1]})->cnv_intervals_index;
	intervals[unlist(cnv_intervals_index),]->res;
	lapply(res$CNV_region,function(x){unlist(strsplit(as.character(x),split="\\."))->res;res[1]})->cnv_genes_chr;
	unlist(cnv_genes_chr)->res$Chr;
	merge(myd,res,by.x="CnvGene",by.y="Gene")->res;
	return(res);
}
CNVCor_genes_count(CNV_region2symbol_EXP.cor_filter1,cnv_interval2genes)->CNVCor_genes_chr_count;
as.numeric(gsub("chr","",CNVCor_genes_chr_count$Chr))->CNVCor_genes_chr_count$ChrN
table(CNVCor_genes_chr_count$ChrN)->chr_table;
chr_table/sum(chr_table)->chr_table_ratio;
paste("Chr",names(chr_table),sep="")->names(chr_table_ratio);
layout(matrix(c(1,2,2),3,1,byrow=T))
barplot(chr_table_ratio,las=2,col=brewer.pal(12,"Set3")[12],cex.names=1.5);
boxplot(Cor~ChrN,data=CNVCor_genes_chr_count,col=brewer.pal(12,"Set3")[12],pch=20,cex.names=2.5,names=paste("Chr",names(chr_table),sep=""),lty=1)
#ggviolin(CNVCor_genes_chr_count,x="Chr",y="Cor",add=c("mean_se"),fill=brewer.pal(11,"Set3")[2],order=c(paste("chr",seq(1,22,1),sep=""),"chrX"))+theme(axis.text.x=element_text(angle=90))+geom_hline(yintercept=c(0.7,0.8),linetype=2)

#--------------------------------------------------------------------------------------------------------------
METCor_genes_count<-function(myd,myd_cor){
	merge(myd_cor,myd[,1:9],by.x="MetProbe",by.y="Probe")->myd_merged;
	table(myd_merged$MetGene)->met_genes;
	met_genes[met_genes!=0]->met_genes;
	names(met_genes)->met_genes;
	c()->filter_index;
	for(g in met_genes){
		for(i in 1:nrow(myd_merged)){
			if(g==myd_merged[i,2]){
				c(filter_index,i)->filter_index;
				break;
			}
		}
	}
	return(myd_merged[filter_index,]);
}
METCor_genes_count(myd.methy,MET_EXP.cor_filter1)->METCor_genes_summary_count;
ggboxplot(METCor_genes_summary_count,x="Chr",y="Cor",fill=brewer.pal(12,"Set3")[12],order=c(paste("chr",seq(1,22,1),sep=""),"chrX"))+theme(axis.text.x=element_text(angle=90))+geom_hline(yintercept=c(-0.5,-0.6),linetype=2)
ggboxplot(METCor_genes_summary_count,x="Chr",y="TSSDist",fill=brewer.pal(12,"Set3")[12],order=c(paste("chr",seq(1,22,1),sep=""),"chrX"))+theme(axis.text.x=element_text(angle=90))+geom_hline(yintercept=c(-200,500,1000,1500),linetype=2)
#par(mfrow=c(2,1))
#boxplot(Cor~Chr,data=METCor_genes_summary_count,col=brewer.pal(11,"Set3")[1],pch=20)
#boxplot(TSSDist~Chr,data=METCor_genes_summary_count,col=brewer.pal(11,"Set3")[2],pch=20)

draw_METCor_genes_summary<-function(myd){
	layout(matrix(c(1,1,1,1,2,2,2,3),2,4,byrow=T));
	table(myd$Chr)->chr_table;
	c(paste("chr",seq(1,22,1),sep=""),"chrX")->chr_order;
	chr_table/sum(chr_table)->chr_table_ratio;
	names(chr_table)->names(chr_table_ratio);
	chr_table_ratio[chr_order]->chr_table_ratio;
	barplot(chr_table_ratio,col=brewer.pal(12,"Set3")[12],cex.names=1.5,las=2);
	table(myd$GeneType)->gtype_table;
	gtype_table[!is.na(names(gtype_table))]->gtype_table;
	gtype_table/sum(gtype_table)->gtype_table_ratio;
	names(gtype_table)->names(gtype_table_ratio);
	gtype_table_ratio[order(gtype_table_ratio,decreasing=T)[1:5]]->gtype_table_ratio;
	barplot(gtype_table_ratio,horiz=T,names.arg=F,axes=F,col=brewer.pal(12,"Set3")[12])->gtype_mp;
	axis(1);
	text(x=0.01,y=gtype_mp,labels=names(gtype_table_ratio),cex=1.5,pos=4);
	#axis(2,at=gtype_mp,mgp=c(0,-15,0),cex=2,tick=F,labels=names(gtype_table_ratio),las=2)
	table(myd$CGI)->cgi_table;
	which(names(cgi_table)%in%"")->null_index;
	c(null_index,which(names(cgi_table)%in%"NA"))->null_index;
	cgi_table[-null_index]->cgi_table;
	cgi_table/sum(cgi_table)->cgi_table_ratio;
	names(cgi_table)->names(cgi_table_ratio);
	barplot(cgi_table_ratio,col=brewer.pal(12,"Set3")[12],names.arg=F,axes=F,horiz=T)->cgitype_mp;
	axis(1);
	text(y=cgitype_mp,x=0.01,labels=names(cgi_table_ratio),cex=1.5,pos=4);

}
draw_METCor_genes_summary(METCor_genes_summary_count)
########################
read.table("D:\\TCGA-OV\\multi-omics\\gencode.v22.annotation.genes",header=T,sep="\t",stringsAsFactors=F)->total_genes;
CNV_region2symbol_EXP.cor_filter1$CnvGene->genes;
lapply(genes,function(x){which(total_genes$GeneName%in%x)[1]})->cnv_intervals_index;
table(total_genes[unlist(cnv_intervals_index),1])->CNVCor_chr_table
total_genes[total_genes$GeneType=="protein_coding",]->total_genes_pro;
table(total_genes_pro$Chr)->total_genes_pro_chr;
do_multiple_ks_test<-function(dat1,dat2){
	dat2[names(dat1)]->dat2_;
	sum(dat1)->dat1_sum;
	sum(dat2_)->dat2_sum;
	c()->chr_pvalue;
	for(i in 1:length(dat1)){
		matrix(c(dat1[i],dat2[i],dat1_sum-dat1[i],dat2_sum-dat2_[i]),nrow=2,byrow=T)->dat1_dat2_matrix;
		fisher.test(dat1_dat2_matrix)->ks_test;#chisq test is also good!
		c(chr_pvalue,ks_test$p.value)->chr_pvalue;
	}
	data.frame("Chr"=names(dat1),"gCount"=as.numeric(dat1),"FisherP"=chr_pvalue)->res_df;
	res_df[order(res_df$FisherP),]->res_df;
	p.adjust(res_df$FisherP)->res_df$Padj;
	return(res_df);
}
paste("chr",names(chr_table),sep="")->names(chr_table)
do_multiple_ks_test(chr_table,total_genes_pro_chr)->CNVCor_chr_table.fisher_test;
do_multiple_ks_test(table(METCor_genes_summary_count$Chr),total_genes_pro_chr)->METCor_chr_table. ;

#------------------------------------
intersect(coxph_filtered_genes,CNVCor_region2symbol_genes)->CNVCor_genes_coxph;
intersect(coxph_filtered_genes,METCor_genes)->METCor_genes_coxph;

CNV_region2symbol_EXP.cor_filter[which(CNV_region2symbol_EXP.cor_filter$CnvGene%in%CNVCor_genes_coxph),]->CNV_region2symbol_EXP.cor_filter1;
MET_EXP.cor_filter[which(MET_EXP.cor_filter$MetGene%in%METCor_genes_coxph),]->MET_EXP.cor_filter1;
########--------NMF cluster:
library(NMF);
library(IntNMF)
#---using CNVCor_genes_coxph------------------------------------------------
myd_exp_processed[,CNVCor_genes_coxph]->CNVCor_genes_EXP;
myd_exp_processed$A0_Samples->rownames(CNVCor_genes_EXP);
nmf(t(CNVCor_genes_EXP),2:10,nrun=50,seed=12345)->CNVCor_genes_nmf;
plot(CNVCor_genes_nmf);
consensusmap(CNVCor_genes_nmf,labCol=NA,labRow=NA,tracks=NA)
#---------------
myd_exp_processed[,METCor_genes_coxph]->METCor_genes_EXP;
myd_exp_processed$A0_Samples->rownames(METCor_genes_EXP)
nmf(t(METCor_genes_EXP),2:10,nrun=50,seed=12345)->METCor_genes_nmf;
plot(METCor_genes_nmf);
consensusmap(METCor_genes_nmf,labCol=NA,labRow=NA,tracks=NA)
################################################################################################
#------------------------------select the best rank: CNV->4,MET->4;using CNVCor_region2symbol_genes
retrive_cluster_names<-function(myd,myd_consensusmap,hvalue){
	rownames(myd)->sample_names;
	lapply(cut(myd_consensusmap$Colv,hvalue)$lower, function(l)rapply(l,function(i)i))->myd_cut_list;
	c()->cluster_sample_names;
	c()->tmp_cluster;
	1->c_index;
	for(i in myd_cut_list){
		c(cluster_sample_names,as.character(sample_names[unlist(i)]))->cluster_sample_names;
		c(tmp_cluster,rep(paste("C",c_index,sep=""),length(unlist(i))))->tmp_cluster;
		c_index+1->c_index;		
	}
	data.frame("Sample"=cluster_sample_names,"Cluster"=tmp_cluster)->cluster_df;
	return(cluster_df);
}
nmf(t(CNVCor_genes_EXP),4,nrun=50,seed=12345)->CNVCor_genes_nmf_3;#5 clusters
consensusmap(CNVCor_genes_nmf_3,labCol=NA,labRow=NA,tracks=NA)->CNVCor_genes_nmf_3_consensusmap
retrive_cluster_names(CNVCor_genes_EXP,CNVCor_genes_nmf_3_consensusmap,0.9)->CNVCor_genes_nmf_3_cluster;
table(CNVCor_genes_nmf_3_cluster$Cluster)
c("Sample","CNVCor_C")->colnames(CNVCor_genes_nmf_3_cluster)
paste("CNVCor",CNVCor_genes_nmf_3_cluster$CNVCor_C,sep="")->CNVCor_genes_nmf_3_cluster$CNVCor_C;

nmf(t(METCor_genes_EXP),4,nrun=50,seed=12345)->METCor_genes_nmf_4;#4 clusters
consensusmap(METCor_genes_nmf_4,labCol=NA,labRow=NA,tracks=NA)->METCor_genes_nmf_4_consensusmap;
#
nmf(t(METCor_genes_EXP),3,nrun=50,seed=12345)->METCor_genes_nmf_3;#3 clusters
consensusmap(METCor_genes_nmf_3,labCol=NA,labRow=NA,tracks=NA)->METCor_genes_nmf_3_consensusmap;

retrive_cluster_names(METCor_genes_EXP,METCor_genes_nmf_3_consensusmap,0.7)->METCor_genes_nmf_3_cluster;
table(METCor_genes_nmf_3_cluster$Cluster)
#"C1"->METCor_genes_nmf_3_cluster[METCor_genes_nmf_3_cluster$Cluster=="C4",2]
c("Sample","METCor_C")->colnames(METCor_genes_nmf_3_cluster)
paste("METCor",METCor_genes_nmf_3_cluster$METCor_C,sep="")->METCor_genes_nmf_3_cluster$METCor_C;

#---for each cluster samples: draw survival curve;
merge(myd.exp_newEvent_processed,CNVCor_genes_nmf_3_cluster,by.x="A0_Samples",by.y="Sample")->myd.exp_filter_CNVCor_cluster;
merge(myd.exp_newEvent_processed,METCor_genes_nmf_3_cluster,by.x="A0_Samples",by.y="Sample")->myd.exp_filter_METCor_cluster;
par(mfrow=c(2,1))
draw_survial_curve_custom(myd.exp_filter_CNVCor_cluster,ncol(myd.exp_filter_CNVCor_cluster),1000,myd_colors)#p<0.01
draw_survial_curve_custom(myd.exp_filter_METCor_cluster,ncol(myd.exp_filter_METCor_cluster),1000,myd_colors)#p=0.011
######################for skewness test: Performs D'Agostino test for skewness in normally distributed data.
library(moments)
agostino.test(CNV_region2symbol_EXP.cor$Z_value, alternative = c("two.sided", "less", "greater"))#p-value < 2.2e-16,skew=0.67815
agostino.test(sample(MET_EXP.cor$Z_value,size=40000), alternative = c("two.sided", "less", "greater"))#p-value < 2.2e-16,skew=-0.27318
#################################################################
library(VennDiagram)
library(ggpubr);
library(limma)
library(ggplotify)
library(gridExtra);
rep("CNV",nrow(CNV_region2symbol_EXP.cor))->CNV_region2symbol_EXP.cor$Type;
rep("MET",nrow(MET_EXP.cor))->MET_EXP.cor$Type;
rbind(CNV_region2symbol_EXP.cor[,3:7],MET_EXP.cor[,3:7])->CNV_MET.cor_combined;
ggdensity(CNV_MET.cor_combined,x="Z_value",fill="Type",add="mean",)->density_p1;
#--------
data.frame("Type"=c(rep("CNVCor",length(CNVCor_genes_coxph)),rep("METCor",length(METCor_genes_coxph))),"Symbol"=c(CNVCor_genes_coxph,METCor_genes_coxph))->CNVCor_METCor.df
table(CNVCor_METCor.df$Symbol,CNVCor_METCor.df$Type)->CNVCor_METCor.df_table;
vennCounts(CNVCor_METCor.df_table)->CNVCor_METCor.df_table.vennCount;
as.grob(~vennDiagram(CNVCor_METCor.df_table.vennCount,names=c("CNVCor","METCor"),cex=1.5,lwd=1.5,circle.col=brewer.pal(9,"Set1")))->venn_p1;
grid.arrange(density_p1,venn_p1,nrow=2);
############################################################################################
####-----compare CNVCor_genes cluster with METCor_genes cluster:
draw_samples_colinear<-function(myd,factors,myd_colors,legend_f){
	if(length(factors)<2){
		return(NULL);
	}
	factors[1]->f1;
	list()->seg_colors;
	#-------
	names(table(myd[,f1]))->x_names;
	c()->x_names_index;
	c()->x_names_colors;
	1->index;
	for(xn in x_names){
		which(myd[,f1]==xn)->xn_index;
		c(x_names_index,xn_index)->x_names_index;
		c(x_names_colors,rep(myd_colors[index],length(xn_index)))->x_names_colors;
		index+1->index;
	}
	myd[x_names_index,]->myd_sort;
	seq(1,nrow(myd_sort),1)->myd_sort$x_pos;
	c(seg_colors,list(x_names_colors))->seg_colors;
	#-------draw 
	5->x_start;5->rect_width;15->xy_width;
	#----------
	seq(1,nrow(myd_sort),1)->counts;
	plot(c(1,60),c(0,nrow(myd_sort)*1.1),type='n',xlab="",ylab="")
	rect(x_start,counts,x_start+rect_width,counts+1,col=x_names_colors,border=NA)
	#---for other factors:
	1->xy_index;
	for(fx in factors[-1]){
		names(table(myd[,fx]))->y_names;
		c()->y_names_pos;
		c()->y_names_colors;
		1->index;
		for(yn in y_names){
			which(myd_sort[,fx]==yn)->yn_index;
			myd_sort$x_pos[yn_index]->yn_pos;
			yn_pos[order(yn_pos)]->yn_pos;
			c(y_names_pos,yn_pos)->y_names_pos;
			c(y_names_colors,rep(myd_colors[index],length(yn_index)))->y_names_colors;
			index+1->index;
		}
		c(seg_colors,list(y_names_colors))->seg_colors;
		rect(x_start+xy_width*xy_index,counts,x_start+xy_width*xy_index+rect_width,counts+1,col=y_names_colors,border=NA);
		#----------x_start;
		myd_sort[y_names_pos,]->myd_sort;
		myd_sort$x_pos->x_names_pos;
		seg_colors[[xy_index]][x_names_pos]->xy_seg_colors;
		#---
		seq(1,nrow(myd_sort),1)->y_names_pos;
		segments(x_start+rect_width+xy_width*(xy_index-1),x_names_pos,x_start+xy_width*xy_index,y_names_pos,col=xy_seg_colors)
		#--
		y_names_pos->myd_sort$x_pos;
		xy_index+1->xy_index;
	}
	legend("top",legend=factors,fill=myd_colors[1:length(factors)],horiz=T,border=NA);
	#------
	table(myd[,legend_f])->legend_f_table;
	names(legend_f_table)->legend_f_names;
	paste(legend_f_names,legend_f_table,sep="(n=")->legend_labels;
	paste(legend_labels,")",sep="")->legend_labels;
	legend("right",legend=legend_labels,fill=myd_colors[1:length(legend_labels)],border=NA,title=legend_f)
	#return(list(x_names_pos,y_names_pos));
}
merge(CNVCor_genes_nmf_3_cluster,METCor_genes_nmf_3_cluster,by.x="Sample",by.y="Sample")->CNVCor_METCor_nmf_cluster
merge(myd_exp_processed[,c("A0_Samples","HistologicalType")],CNVCor_METCor_nmf_cluster,by.x="A0_Samples",by.y="Sample")->CNVCor_METCor_nmf_cluster;
draw_samples_colinear(CNVCor_METCor_nmf_cluster,c("CNVCor_C","HistologicalType","METCor_C"),myd_colors,"HistologicalType")->y_res;
#table(CNVCor_METCor_nmf_cluster$CNVCor_C,CNVCor_METCor_nmf_cluster$METCor_C)->CNVCor_METCor_nmf_cluster_table;
#layout(matrix(c(1,1,2,3,3,3),nrow=2,byrow=T))
#barplot(apply(CNVCor_METCor_nmf_cluster_table,2,function(x){x/sum(x)}),col=brewer.pal(9,"Set1")[1:4],beside=F)
#plot.new()
#legend("topleft",legend=rownames(CNVCor_METCor_nmf_cluster_table),fill=brewer.pal(9,"Set1")[1:4])
library("gplots")
balloonplot(CNVCor_METCor_nmf_cluster_table,main="CNVCor genes clustering subset \n overlap with METCor genes clustering subset",xlab ="", ylab="",label = FALSE, show.margins = FALSE)
chisq.test(CNVCor_METCor_nmf_cluster_table)$p.value->table_pvalue;
if(table_pvalue<1e-5){
	legend("topleft",legend=paste("chisq-p:","<1e-5"))
}else{
	legend("topleft",legend=paste("chisq-p:",round(table_pvalue,5)))
}

##################################################################################################################################
##################################################################################################################################
library(gplots)
library(iClusterPlus);
myd.cnv_region2symbol_logTrans[,as.character(CNVCor_genes_coxph)]->myd.cnv_merged_CNVCor;
myd.cnv_region2symbol_logTrans$Sample->rownames(myd.cnv_merged_CNVCor);
myd.cnv_merged_CNVCor[myd_samples,]->myd.cnv_merged_CNVCor;

t(myd.methy[,10:ncol(myd.methy)])->myd.methy_t;
myd.methy$Probe->colnames(myd.methy_t);
gsub("\\.","-",rownames(myd.methy_t),perl=T)->rownames(myd.methy_t);
myd.methy_t[,as.character(MET_EXP.cor_filter1$MetProbe)]->myd.methy_t_METCor;
myd.methy_t_METCor[myd_samples,]->myd.methy_t_METCor;

unique(c(CNVCor_genes_coxph,METCor_genes_coxph))->CNVCor_METCor_genes;
myd_exp_processed[,CNVCor_METCor_genes]->CNVCor_METCor_genes_EXP;
myd_exp_processed$A0_Samples->rownames(CNVCor_METCor_genes_EXP)
CNVCor_METCor_genes_EXP[myd_samples,]->CNVCor_METCor_genes_EXP;
save.image("multi-omics/work-iclusterplus2.RData")#
#------------------------------------single lambda iclusterplusï¼5 clusters or 6 clusters
get_best_lambda<-function(myd_tune){
	myd_tune$fit->myd_tune_fit;
	myd_tune$lambda->myd_tune_lambda;
	which.min(unlist(lapply(myd_tune_fit,function(x){x$BIC})))->min_BIC_index;
	myd_tune_lambda[min_BIC_index,]->myd_tune_best_lambda;
	return(myd_tune_best_lambda);
}
lapply(2:6,function(x){paste("icluster_tune_fit",x,sep="")->x_fit;get_best_lambda(eval(as.symbol(x_fit)))})#2: 0.651351351 0.002702703 0.856756757
##K=2:0.305405405 0.008108108 0.716216216
#K=3:0.61351351 0.01891892 0.43513514
c()->fit_list1;
for(i in 1:20){
	set.seed(i*100+i*3+i);
	iClusterPlus(dt1=as.matrix(myd.cnv_merged_CNVCor),dt2=as.matrix(myd.methy_t_METCor),dt3=as.matrix(CNVCor_METCor_genes_EXP),type=c("gaussian","gaussian","gaussian"),lambda=c(0.305405405,0.008108108,0.716216216),K=2,maxiter=20)->single_fit;#three clusters
	c(fit_list1,list(single_fit))->fit_list1;
}
c()->fit_list2;
for(i in 1:20){
	set.seed(i*100+i*3+i);
	iClusterPlus(dt1=as.matrix(myd.cnv_merged_CNVCor),dt2=as.matrix(myd.methy_t_METCor),dt3=as.matrix(CNVCor_METCor_genes_EXP),type=c("gaussian","gaussian","gaussian"),lambda=c(0.61351351,0.01891892,0.43513514),K=3,maxiter=20)->single_fit;#four clusters
	c(fit_list2,list(single_fit))->fit_list2;
}
generate_icluster_groups<-function(icluster_fit){
	icluster_fit$cluster->ic_cluster;
	table(ic_cluster)->ic_table;
	paste("iC",names(ic_table),sep="")->ic_table_groups;
	names(ic_table)[order(ic_table)]->ic_table_names;
	c()->new_clusters;
	for(j in 1:length(ic_cluster)){
		for(i in 1:length(ic_table_names)){
			as.numeric(ic_table_names[i])->ic;
			if(ic_cluster[j]==ic){
				c(new_clusters,ic_table_groups[i])->new_clusters;
			}
		}
	}
	return(new_clusters);
}
draw_icluster_repeat20_survival_curve<-function(fit_list,bk,cols){
	c()->rank_repeat20_summary;
	par(mfrow=c(4,5))
	for(i in 1:20){
		fit_list[[i]]->single_fit;
		generate_icluster_groups(single_fit)->single_fit_ic;
		c(rank_repeat20_summary,table(single_fit_ic))->rank_repeat20_summary;
		data.frame("Sample"=rownames(myd.cnv_merged_CNVCor),single_fit_ic)->iCluster_df;
		merge(myd.exp_newEvent_processed,iCluster_df,by.x="A0_Samples",by.y="Sample")->iCluster_merged;
		draw_survial_curve_custom(iCluster_merged,ncol(iCluster_merged),bk,cols)
	}
	matrix(rank_repeat20_summary,ncol=length(table(single_fit_ic)),byrow=T)->rank_repeat20_summary;
	names(table(single_fit_ic))->colnames(rank_repeat20_summary);
	data.frame("Repeat"=paste("Rep",1:20,sep=""),rank_repeat20_summary)->rank_repeat20_summary;
	return(rank_repeat20_summary);
}
draw_icluster_repeat20_survival_curve(fit_list1,1000,myd_colors)->rank2_repeat20_summary;#three clusters,,,
draw_icluster_repeat20_survival_curve(fit_list2,1000,myd_colors)->rank3_repeat20_summary;#four clusters;not good
#k=3: 4 clusters, i=1;
fit_list2[[1]]->single_fit;
data.frame("Sample"=rownames(myd.cnv_merged_CNVCor),"iCluster"=generate_icluster_groups(single_fit))->iCluster_df;
change_group<-function(iclusters,oldGroups,newGroups){
	c()->new_clusters;
	for(i in 1:nrow(iclusters)){
		for(ic in 1:length(oldGroups)){
			if(iclusters[i,2]==oldGroups[ic]){
				c(new_clusters,newGroups[ic])->new_clusters;
			}
		}
	}
	new_clusters->iclusters$iCluster;
	return(iclusters);
}
merge(myd.exp_newEvent_processed,iCluster_df,by.x="A0_Samples",by.y="Sample")->iCluster_merged;
par(mfrow=c(2,2))
draw_survial_curve_custom(iCluster_merged,ncol(iCluster_merged),1000,myd_colors)#0.0104
subset_cluster<-function(myd,sub_clusters){
	c()->sub_clusters_index;
	for(i in sub_clusters){
		which(myd$iCluster==i)->tmp_index;
		c(sub_clusters_index,tmp_index)->sub_clusters_index;
	}
	data.frame(myd[sub_clusters_index,])->myd_sub;
	as.character(myd_sub$iCluster)->myd_sub$iCluster;
	return(myd_sub);
}
subset_cluster(iCluster_merged,c("iC2","iC1"))->test_;#log rank p-value:045549
draw_survial_curve_custom(test_,ncol(iCluster_merged),1000,myd_colors[c(1,2)])
subset_cluster(iCluster_merged,c("iC2","iC3"))->test_;#log rank p-value:0.00138
draw_survial_curve_custom(test_,ncol(iCluster_merged),1000,myd_colors[c(2,3)])
subset_cluster(iCluster_merged,c("iC2","iC4"))->test_;#log rank p-value:0.00409
draw_survial_curve_custom(test_,ncol(iCluster_merged),1000,myd_colors[c(2,4)])
#---merge CNVCluster, METCluster and iCluster:
merge(CNVCor_METCor_nmf_cluster,iCluster_df,by.x="A0_Samples",by.y="Sample")->CNVCor_METCor_iC_cluster;
draw_samples_colinear(CNVCor_METCor_iC_cluster,c("CNVCor_C","HistologicalType","iCluster","METCor_C"),myd_colors,"HistologicalType")->y_res;
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------
preoder_samples<-function(myd,sample_clusters.df,cluster_type){
	c()->myd_order;
	myd[as.character(sample_clusters.df$Sample),]->myd;
	as.numeric(substr(sample_clusters.df$iCluster,3,3))->sample_clusters;
	for(i in cluster_type){
		which(sample_clusters==i)->tmp_index;
		hclust(dist(myd[tmp_index,]))->tmp_index_hclust;
		tmp_index[tmp_index_hclust$order]->tmp_index;
		c(myd_order,tmp_index)->myd_order;	
	}
	return(myd[myd_order,]);
}
CNVCor_METCor_iC_cluster[,2:5]->col_CNVCor_METCor_iC_cluster
CNVCor_METCor_iC_cluster$A0_Samples->rownames(col_CNVCor_METCor_iC_cluster);
generate_color<-function(x){
	brewer.pal(9,"Set1")->myd_colors;
	names(table(x))->x_table;
	myd_colors[1:length(x_table)]->x_colors;
	x_table->names(x_colors);
	return(x_colors);
}
generate_gaps<-function(myd_clusters){
	table(myd_clusters)->cluster_count;
	c()->myd_gaps;
	0->tmp_gap;
	for(i in 1:(length(cluster_count)-1)){
		tmp_gap+cluster_count[i]->tmp_gap;
		c(myd_gaps,tmp_gap)->myd_gaps;
	}
	return(myd_gaps);
}
generate_gaps(CNVCor_METCor_iC_cluster$iCluster)->gaps_col;
list("HistologicalType"=generate_color(CNVCor_METCor_iC_cluster[,2]),"CNVCor_C"=generate_color(CNVCor_METCor_iC_cluster[,3]),"METCor_C"=generate_color(CNVCor_METCor_iC_cluster[,4]),"iCluster"=generate_color(CNVCor_METCor_iC_cluster[,5]))->col_colors;

preoder_samples(myd.cnv_merged_CNVCor,iCluster_df,c(1,2,3,4))->test1;
#myd.cnv_merged_CNVCor
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(250),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(150))->fill_colors;
pheatmap(t(log2(test1)),cluster_rows=T,cluster_cols=F,show_rownames=F,show_colnames=F,color=fill_colors,gaps_col=gaps_col,annotation_col=col_CNVCor_METCor_iC_cluster,annotation_colors=col_colors)->test1_p1#CNV
#myd.methy_t_METCor
preoder_samples(myd.methy_t_METCor,iCluster_df,c(1,2,3,4))->test2;
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(240),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(150))->fill_colors;
pheatmap(t(test2),cluster_rows=T,cluster_cols=F,show_rownames=F,show_colnames=F,color=fill_colors,gaps_col=gaps_col,annotation_col=col_CNVCor_METCor_iC_cluster,annotation_colors=col_colors)->test2_p1#Methylation
#-
grid.arrange(test1_p1$gtable,test2_p1$gtable,nrow=2)
#######################################################################################################
#-------------------------------------for CNV genes gain or loss, MET hyper or hypo
#myd.cnv_region2symbol[,CNVCor_genes_coxph]->myd.cnv_region2symbol_CNVCor;
apply(myd.cnv_region2symbol[,-1],1,function(x){c(length(x[x>=0.3]),length(x[x<=(-0.3)]))})->myd.cnv_merged_gainLoss;
as.data.frame(t(myd.cnv_merged_gainLoss))->myd.cnv_merged_gainLoss;
c("Gain","Loss")->colnames(myd.cnv_merged_gainLoss);
myd.cnv_region2symbol$Sample->myd.cnv_merged_gainLoss$Samples;
#for Methylation:
t(myd.methy[,10:ncol(myd.methy)])->myd.methy_t;
myd.methy$Probe->colnames(myd.methy_t);
gsub("\\.","-",rownames(myd.methy_t),perl=T)->rownames(myd.methy_t);
apply(myd.methy_t,1,function(x){c(length(x[x>0.8]),length(x[x<0.2]))})->myd.methy_t_hyperHypo;
as.data.frame(t(myd.methy_t_hyperHypo))->myd.methy_t_hyperHypo;
c("MetHyper","MetHypo")->colnames(myd.methy_t_hyperHypo)
rownames(myd.methy_t)->myd.methy_t_hyperHypo$Samples;
merge(myd.cnv_merged_gainLoss,myd.methy_t_hyperHypo,by.x="Samples",by.y="Samples")->myd.CNV_MET_abnormal_frequency;
write.table(myd.CNV_MET_abnormal_frequency,"multi-omics/CNV_MET_abnormal_frequency.txt",quote=F,sep="\t",row.names=F)
#_-------draw points :
draw_genes_exp_point_v1<-function(expd,g1,g2){
	plot(log2(expd[,g1]+1),log2(expd[,g2]+1),xlab=paste(c("log2(",g1,")"),collapse=""),ylab=paste(c("log2(",g2,")"),collapse=""),pch=20,cex=0.5,bty="o");
	#-------------
	expd[,c(g1,g2)]->expd_g1_g2;
	expd_g1_g2[expd_g1_g2[,1]!=0,]->expd_g1_g2;
	expd_g1_g2[expd_g1_g2[,2]!=0,]->expd_g1_g2;
	log2(expd_g1_g2+1)->expd_g1_g2;
	#------------------
	cor.test(expd_g1_g2[,g1],expd_g1_g2[,g2],method = c("spearman"))->g1_g2_cortest;
	g1_g2_cortest$estimate->cortest_cor;
	g1_g2_cortest$p.value->cortest_p;
	#--------
	paste(g2,g1,sep="~")->g1_g2_formula;
	as.formula(g1_g2_formula)->g1_g2_formula;
	lm(g1_g2_formula,data=expd_g1_g2)->g1_g2_lmfit;
	abline(g1_g2_lmfit,lty=3,col="red",lwd=2)
	if(cortest_p < 0.001){
		"p < 0.001"->cortest_p
	}else{
		paste("p = ",round(cortest_p,3),sep="")->cortest_p;
	}
	paste("R = ",round(cortest_cor,3),sep="")->cortest_cor;
	paste(cortest_cor,cortest_p,sep=",")->leg;
	mtext(leg,cex=1.1,font=3);
}
calculate_correlation_fdr<-function(expd,genes,f){
	c()->tmp_cor;
	c()->tmp_p;
	for(g in genes){
		#-------------
		expd[,c(f,g)]->expd_g1_g2;
		expd_g1_g2[expd_g1_g2[,1]!=0,]->expd_g1_g2;#remove x=0 
		expd_g1_g2[expd_g1_g2[,2]!=0,]->expd_g1_g2;#remove y=0
		cor.test(expd_g1_g2[,f],expd_g1_g2[,g],method = c("spearman"))->g1_g2_cortest;
		g1_g2_cortest$estimate->cortest_cor;
		g1_g2_cortest$p.value->cortest_p;
		c(tmp_cor,cortest_cor)->tmp_cor;
		c(tmp_p,cortest_p)->tmp_p;
	}
	data.frame("genes"=paste(genes,f,sep="~"),"Corr"=tmp_cor,"CorrP"=tmp_p)->res;
	res[!is.na(res$CorrP),]->res;
	res[order(res$CorrP),]->res;
	p.adjust(res$CorrP)->res$FDR;
	return(res);
}
draw_genes_exp_point_fdr<-function(expd,g1,g2,cortest){
	plot(log2(expd[,g1]+1),log2(expd[,g2]+1),xlab=paste(c("log2(",g1,")"),collapse=""),ylab=paste(c("log2(",g2,")"),collapse=""),pch=20,cex=0.5,bty="o");
	paste(g2,g1,sep="~")->g1_g2_names;
	cortest[which(cortest$genes==g1_g2_names),2]->cortest_cor;
	cortest[which(cortest$genes==g1_g2_names),4]->cortest_padj;
	#-------------
	expd[,c(g1,g2)]->expd_g1_g2;
	expd_g1_g2[expd_g1_g2[,1]!=0,]->expd_g1_g2;#remove x=0 
	expd_g1_g2[expd_g1_g2[,2]!=0,]->expd_g1_g2;#remove y=0
	log2(expd_g1_g2+1)->expd_g1_g2;
	paste(g2,g1,sep="~")->g1_g2_formula;
	as.formula(g1_g2_formula)->g1_g2_formula;
	lm(g1_g2_formula,data=expd_g1_g2)->g1_g2_lmfit;
	abline(g1_g2_lmfit,lty=3,col="red",lwd=2)
	if(cortest_padj < 0.001){
		"FDR < 0.001"->cortest_padj
	}else{
		paste("FDR = ",round(cortest_padj,3),sep="")->cortest_padj;
	}
	paste("Rho = ",round(cortest_cor,3),sep="")->cortest_cor;
	paste(cortest_cor,cortest_padj,sep=",")->leg;
	mtext(leg,cex=1,font=3);
}
c("Gain","Loss","MetHyper","MetHypo")->selected_geneSet;
data.frame("genes"=NULL,"Corr"=NULL,"CorrP"=NULL,"FDR"=NULL)->expd_cortest;
for(g1 in selected_geneSet){
	calculate_correlation_fdr(myd.CNV_MET_abnormal_frequency,selected_geneSet,g1)->tmp_cortest;#-------for Gain vs others
	rbind(expd_cortest,tmp_cortest)->expd_cortest;
}
#corrplot(corr_table.cor,p.mat=corr_table.res1$p,order="hclust",addrect=3,sig.level=0.05,insig="blank",tl.col="black",tl.srt=45,col=rev(col2(200)))->corr_p;
par(mfrow=c(2,3))
draw_genes_exp_point_fdr(myd.CNV_MET_abnormal_frequency,"Loss","Gain",expd_cortest);
draw_genes_exp_point_fdr(myd.CNV_MET_abnormal_frequency,"MetHyper","Gain",expd_cortest);
draw_genes_exp_point_fdr(myd.CNV_MET_abnormal_frequency,"MetHypo","Gain",expd_cortest);
draw_genes_exp_point_fdr(myd.CNV_MET_abnormal_frequency,"MetHyper","Loss",expd_cortest);
draw_genes_exp_point_fdr(myd.CNV_MET_abnormal_frequency,"MetHypo","Loss",expd_cortest);
draw_genes_exp_point_fdr(myd.CNV_MET_abnormal_frequency,"MetHypo","MetHyper",expd_cortest);
#######################################################################################################################
#---------------immune landscape analysis with clinical features:
#######################################################################################################################
#-----------------------------immune landscape for C1-C3: immune infiltrating,Leukocyte fraction,SNV-Neoantigen,Indel-Neoantigen
#---------draw scores boxplot:
draw_boxplot_genes_by_factors<-function(expd,genes,f,myd_colors){
	expd[,f]->f_values;
	as.character(genes)->genes;
	intersect(genes,colnames(expd))->genes;
	table(expd[,f])->f_table;
	names(f_table)->f_names;
	paste(f_names,f_table,sep="(")->f_names_legend;
	paste(f_names_legend,")",sep="")->f_names_legend;
	lapply(f_names,function(x){which(f_values==x)->res;res})->f_names_index;
	lapply(genes,function(x){expd[,x]->x_values;lapply(f_names_index,function(y){x_values[y]->res;res})})->genes_f_values;
	list()->genes_f_values_list;
	c()->rank_p_values;
	for(i in 1:length(genes)){
		c()->f_levels;
		c()->f_level_values;
		for(j in 1:length(f_names)){
			c(f_levels,rep(f_names[j],length(genes_f_values[[i]][[j]])))->f_levels;
			c(f_level_values,genes_f_values[[i]][[j]])->f_level_values;
			c(genes_f_values_list,list(genes_f_values[[i]][[j]]))->genes_f_values_list
		}
		data.frame("F"=f_levels,"V"=f_level_values)->f_level_df;
		factor(f_level_df[,1],levels=f_names)->f_level_df[,1];
		paste(c("V","F"),collapse="~")->tmp_formula;
		c(rank_p_values,kruskal.test(as.formula(tmp_formula),data=f_level_df)$p.value)->rank_p_values;
	}
	c()->sig_symbol;
	for(x in rank_p_values){
		if(x <= 1e-5){
			c(sig_symbol,"***")->sig_symbol
		}else if(x > 1e-5 && x <= 0.01){
			c(sig_symbol,"**")->sig_symbol
		}else if(x > 0.01 && x <= 0.05){
			c(sig_symbol,"*")->sig_symbol
		}else{
			c(sig_symbol,"")->sig_symbol
		}
	}
	data.frame("Factors"=genes,"KruskalP"=rank_p_values)->kruskal_test_res;
	seq(length(f_names)+1,(length(genes)+1)*(length(f_names)+1)-1,length(f_names)+1)->seg_index;
	seg_index-1->seg_end_index;
	seg_index-length(f_names)->seg_start_index;
	c()->x_at;
	c()->x_label_at;
	for(i in 1:length(seg_start_index)){
		seq(seg_start_index[i],seg_end_index[i])->res;
		c(x_label_at,mean(res))->x_label_at;
		c(x_at,res)->x_at;
	}
	boxplot(genes_f_values_list,plot=F)->genes_f_values.boxplot;
	min(genes_f_values.boxplot$stats)->boxplot_min;
	max(genes_f_values.boxplot$stats)->boxplot_max;
	boxplot(genes_f_values_list,at=x_at,axes=F,ylab="Values",boxwex=0.5,pch=20,cex=0.5,col=myd_colors[1:length(f_names)],ylim=c(round(boxplot_min),boxplot_max)*1.5)#->genes_f_values.boxplot;
	(boxplot_max-round(boxplot_min))/5->bk;
	axis(side=1,at=x_label_at,labels=genes,tick=T,las=2)
	if(is.wholenumber(boxplot_max) && abs(boxplot_max) < 1){
		axis(side=2,at=seq(round(boxplot_min),boxplot_max*1.1,bk),labels=seq(round(boxplot_min),boxplot_max*1.1,bk))
	}else{
		axis(side=2,at=seq(round(boxplot_min),boxplot_max*1.1,bk),labels=round(seq(round(boxplot_min),boxplot_max*1.1,bk),3))
	}
	text(x=x_label_at,y=boxplot_max*1.3,labels=sig_symbol,cex=1.3)
	print(bk);
	flush.console();
	#legend(x=mean(seg_index),y=round(boxplot_max)*1.4,horiz=T,fill=myd_colors[1:length(f_names)],legend=f_names,title=f,inset=0,xjust=0.5,yjust=0.5)
	legend("top",horiz=T,fill=myd_colors[1:length(f_names)],legend=f_names_legend,title=f,inset=0,xjust=0.5,yjust=0.5)
	legend("topright",horiz=T,legend=c("*** p<1e-5","** p<0.01","* p<0.05"),title="Kruskal-Wallis Test")
	box();
	return(kruskal_test_res);
}
draw_boxplot_genes_by_factors_v2<-function(expd,genes,f,myd_colors,log_trans){
	expd[,f]->f_values;
	as.character(genes)->genes;
	intersect(genes,colnames(expd))->genes;
	table(expd[,f])->f_table;
	names(f_table)->f_names;
	paste(f_names,f_table,sep="(")->f_names_legend;
	paste(f_names_legend,")",sep="")->f_names_legend;
	lapply(f_names,function(x){which(f_values==x)->res;res})->f_names_index;
	if(log_trans=="no"){
		lapply(genes,function(x){expd[,x]->x_values;lapply(f_names_index,function(y){x_values[y]->res;res})})->genes_f_values;
	}else{
		lapply(genes,function(x){log2(expd[,x]+1)->x_values;lapply(f_names_index,function(y){x_values[y]->res;res})})->genes_f_values;
	}
	list()->genes_f_values_list;
	c()->rank_p_values;
	for(i in 1:length(genes)){
		c()->f_levels;
		c()->f_level_values;
		for(j in 1:length(f_names)){
			c(f_levels,rep(f_names[j],length(genes_f_values[[i]][[j]])))->f_levels;
			c(f_level_values,genes_f_values[[i]][[j]])->f_level_values;
			c(genes_f_values_list,list(genes_f_values[[i]][[j]]))->genes_f_values_list
		}
		data.frame("F"=f_levels,"V"=f_level_values)->f_level_df;
		factor(f_level_df[,1],levels=f_names)->f_level_df[,1];
		paste(c("V","F"),collapse="~")->tmp_formula;
		c(rank_p_values,kruskal.test(as.formula(tmp_formula),data=f_level_df)$p.value)->rank_p_values;
	}
	c()->sig_symbol;
	for(x in rank_p_values){
		if(x <= 1e-5){
			c(sig_symbol,"***")->sig_symbol
		}else if(x > 1e-5 && x <= 0.01){
			c(sig_symbol,"**")->sig_symbol
		}else if(x > 0.01 && x <= 0.05){
			c(sig_symbol,"*")->sig_symbol
		}else{
			c(sig_symbol,"")->sig_symbol
		}
	}
	data.frame("Factors"=genes,"KruskalP"=rank_p_values)->kruskal_test_res;
	seq(length(f_names)+1,(length(genes)+1)*(length(f_names)+1)-1,length(f_names)+1)->seg_index;
	seg_index-1->seg_end_index;
	seg_index-length(f_names)->seg_start_index;
	c()->x_at;
	c()->x_label_at;
	for(i in 1:length(seg_start_index)){
		seq(seg_start_index[i],seg_end_index[i])->res;
		c(x_label_at,mean(res))->x_label_at;
		c(x_at,res)->x_at;
	}
	boxplot(genes_f_values_list,plot=F)->genes_f_values.boxplot;
	min(genes_f_values.boxplot$stats)->boxplot_min;
	max(genes_f_values.boxplot$stats)->boxplot_max;
	boxplot(genes_f_values_list,at=x_at,axes=F,ylab="Values",boxcol=myd_colors[1:length(f_names)],boxwex=0.5,pch=20,xaxt="n",lty=1,cex=0.5,col=myd_colors[1:length(f_names)],ylim=c(round(boxplot_min),round(boxplot_max)*1.5))#->genes_f_values.boxplot;
	(round(boxplot_max)-round(boxplot_min))/5->bk;
	axis(side=1,at=x_label_at,labels=genes,tick=T,las=2)
	axis(side=2,at=seq(round(boxplot_min),round(boxplot_max)*1.1,bk),labels=seq(round(boxplot_min),round(boxplot_max)*1.1,bk))
	#legend(x=mean(seg_index),y=round(boxplot_max)*1.2,horiz=T,fill=myd_colors[1:length(f_names)],legend=f_names,title=f,inset=0,xjust=0.5,yjust=0.5)
	text(x=x_label_at,y=boxplot_max*1.3,labels=sig_symbol,cex=1.3)
	legend("top",horiz=T,fill=myd_colors[1:length(f_names)],legend=f_names_legend,title=f,inset=0,xjust=0.5,yjust=0.5)
	legend("topright",horiz=T,legend=c("*** p<1e-5","** p<0.01","* p<0.05"),title="Kruskal-Wallis Test")
	#abline(v=seg_index,lty=3)
	box();
	return(kruskal_test_res);
}
#------------
draw_boxplot_factors_by_gene_v1<-function(expd,clusters,f,f_column,myd_colors){
	which(expd[,f_column]==f)->f_index;
	expd[f_index,]->expd;
	list()->f_values;
	for(cl in clusters){
		c(f_values,list(expd[,cl]))->f_values;
	}
	clusters->names(f_values)
	data.frame("F"=rep(clusters,each=length(f_index)),"V"=unlist(f_values))->f_level_df;
	factor(f_level_df[,1],levels=clusters)->f_level_df[,1];
	paste(c("V","F"),collapse="~")->tmp_formula;
	kruskal.test(as.formula(tmp_formula),data=f_level_df)$p.value->rank_p_values;
	#------------------
	rank_p_values->x;
	if(x <= 1e-3){
		"p < 0.001"->sig_symbol
	}else{
		round(x,4)->sig_symbol
	}
	#-----------
	paste(f," Count: ",length(f_index),sep="")->f_count;
	paste("Rank Test ",sig_symbol,sep="")->f_pvalue;
	boxplot(f_values,col=myd_colors[1:length(clusters)],boxcol=myd_colors[1:length(clusters)],boxwex=0.4,ylab="log2(FPKM)",pch=20,cex=0.7,xaxt="n",lty=1)
	#axis(1,at=seq(1,length(clusters),1),labels=clusters,tick=F,cex=0.3,las=2)
	legend("top",legend=clusters,fill=myd_colors[1:length(clusters)],bty='n',horiz=T)
	mtext(text=paste(f_count,f_pvalue,sep=","),font=3);
}
read.table("d:/script/data-set/TIMER-data/TCGA-immuneEstimation.txt",header=T,sep="\t",stringsAsFactors=F)->immuneEstimation;
merge(immuneEstimation,CNVCor_METCor_iC_cluster,by.x="barcode",by.y="A0_Samples")->iC_cluster_immunescore;
write.table(iC_cluster_immunescore,"multi-omics/iC_cluster_immunescore.txt",quote=F,sep="\t",row.names=F)

#-------heatmap:
draw_immune_score_heatmap<-function(myd){
	brewer.pal(9,"Set1")->myd_colors;
	names(table(myd$iCluster))->iC_names;
	c()->myd_order;
	for(i in iC_names){
		which(myd$iCluster==i)->tmp_index;
		hclust(dist(myd[tmp_index,2:7]))->tmp_index_hclust;
		tmp_index[tmp_index_hclust$order]->tmp_index;
		c(myd_order,tmp_index)->myd_order;	
	}
	myd[myd_order,]->myd_sort;
	return(myd_sort)
}
draw_immune_score_heatmap(iC_cluster_immunescore)->iC_cluster_immunescore_sort;

list("HistologicalType"=generate_color(iC_cluster_immunescore_sort[,8]),"CNVCor_C"=generate_color(iC_cluster_immunescore_sort[,9]),"METCor_C"=generate_color(iC_cluster_immunescore_sort[,10]),"iCluster"=generate_color(iC_cluster_immunescore_sort[,11]))->col_colors;
generate_gaps(as.numeric(substr(iC_cluster_immunescore$iCluster,3,3)))->gaps_col;
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(100),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(200))->fill_colors;
pheatmap(t(iC_cluster_immunescore_sort[,2:7]),cluster_rows=F,color=fill_colors,cluster_cols=F,gaps_col=gaps_col,show_colnames=F,cellheight=30,border_color=NA,annotation_col=iC_cluster_immunescore_sort[,8:11],annotation_colors=col_colors)->immuneScore_p1;
grid.arrange(immuneScore_p1$gtable,ncol=2)
#par(mar=c(5,4,4,4),mfrow=c(1,2))
colnames(iC_cluster_immunescore)[2:7]->genes;
draw_boxplot_genes_by_factors_v2(iC_cluster_immunescore,genes,"iCluster",myd_colors,"no")
#---------------------
read.table("d:/script/data-set/immune-landscape.txt",header=T,sep="\t",stringsAsFactors=F)->immune_landscape;
retrive_TCGA_immuneValues<-function(samples,immuneSig){
	unlist(strsplit(as.character(samples[1]),split="-"))->s_split;
	if(length(s_split)>3){
		unlist(lapply(samples,function(x){as.character(x)->x;substr(x,1,nchar(x)-3)->s_res;s_res}))->tmp_samples;
	}
	tmp_samples->names(samples);
	unlist(lapply(tmp_samples,function(x){which(immune_landscape[,1]==x)->x_index;x_index}))->s_index;
	immune_landscape[s_index,immuneSig]->res;
	immune_landscape[s_index,1]->immune_landscape_s;
	data.frame(samples[immune_landscape_s],res)->res;
	c("Sample",immuneSig)->colnames(res);
	return(res);
}
#retrive_TCGA_immuneValues(CNVCor_METCor_iC_cluster$Sample,c("Leukocyte.Fraction","BCR.Shannon","TCR.Shannon"))->sample_immune_landscape;
retrive_TCGA_immuneValues(CNVCor_METCor_iC_cluster$A0_Samples,colnames(immune_landscape)[-c(1,2)])->sample_immune_landscape;
merge(sample_immune_landscape,CNVCor_METCor_iC_cluster,by.x="Sample",by.y="A0_Samples")->iCluster_landscape;
write.table(iCluster_landscape,"multi-omics/iCluster_landscape.txt",quote=F,sep="\t",row.names=F)
colnames(iCluster_landscape)[9:13]->genes;
par(mar=c(10,4,4,4))
draw_boxplot_genes_by_factors_v2(iCluster_landscape,genes,"iCluster",myd_colors,"no")
#grep("Immune.Subtype",colnames(iCluster_landscape))
#draw_survial_curve_custom(iCluster_landscape,2,500,myd_colors)#
###################################################################################
#---compare clinical features in iC3 and iC4:-------------------------------------
library(dgof)
factors_barplot_clustering<-function(myd,sample_clusters,column){
	merge(myd,sample_clusters,by.x="A0_Samples",by.y="A0_Samples")->myd.merged;
	myd.merged[!is.na(myd.merged[column]),]->myd.merged;
	table(myd.merged$iCluster,myd.merged[,column])->riskType_vs_cluster.table;
	dim(riskType_vs_cluster.table)->table_size;
	print(riskType_vs_cluster.table)
	#ks.test(riskType_vs_cluster.table[1,],riskType_vs_cluster.table[2,])$p.value->riskType_vs_cluster.pvalue;
	chisq.test(riskType_vs_cluster.table[1,],riskType_vs_cluster.table[2,])$p.value->riskType_vs_cluster.pvalue;
	max(riskType_vs_cluster.table)*0.9->max_ylim;
	barplot(riskType_vs_cluster.table,beside=T,border=F,col=brewer.pal(9,"Set1")[1:table_size[1]],ylab="Count")#for high and low risk boxplot 
	legend("topright",legend=rownames(riskType_vs_cluster.table),fill=brewer.pal(9,"Set1")[1:table_size[2]],horiz=F)
	text(x=1,y=max_ylim,labels=paste("chisq p:",round(riskType_vs_cluster.pvalue,5),sep=""),pos=4)
	mtext(colnames(myd)[column])
}

factors_barplot_ratio_clustering<-function(myd,sample_clusters,column){
	merge(myd,sample_clusters,by.x="A0_Samples",by.y="A0_Samples")->myd.merged;
	myd.merged[!is.na(myd.merged[column]),]->myd.merged;
	table(myd.merged$iCluster,myd.merged[,column])->riskType_vs_cluster.table;
	dim(riskType_vs_cluster.table)->table_size;
	c()->table_ratio;
	for(i in 1:nrow(riskType_vs_cluster.table)){
		for(j in 1:ncol(riskType_vs_cluster.table)){
			riskType_vs_cluster.table[i,j]/sum(riskType_vs_cluster.table[i,])->res;
			c(table_ratio,res)->table_ratio;
		}
	}
	matrix(table_ratio,nrow=table_size[2],byrow=F)->table_ratio;
	rownames(riskType_vs_cluster.table)->colnames(table_ratio)
	colnames(riskType_vs_cluster.table)->rownames(table_ratio)
	print(table_ratio)
	ks.test(riskType_vs_cluster.table[1,],riskType_vs_cluster.table[2,])$p.value->riskType_vs_cluster.pvalue;
	#chisq.test(riskType_vs_cluster.table)$p.value->riskType_vs_cluster.pvalue;
	barplot(table_ratio,beside=F,border=F,col=brewer.pal(9,"Set1")[1:table_size[2]],ylab="Ratio")#for high and low risk boxplot 
	legend(x="bottom",legend=colnames(riskType_vs_cluster.table),fill=brewer.pal(9,"Set1")[1:table_size[2]],horiz=T,bty="n",inset=c(0,-0.15),xpd=T)
	paste("Chisq-p:",round(riskType_vs_cluster.pvalue,5),sep="")->chisq_p;
	mtext(paste(colnames(myd)[column],chisq_p,sep="\n"))
}
#factors_barplot_ratio_clustering(tmp.filter,CNVCor_METCor_iC_cluster,4)
#---change_values for clinical features;
#merge(myd_exp_processed,CNVCor_METCor_iC_cluster,by.x="A0_Samples",by.y="Sample")->myd_exp_processed.merged;
myd_exp_processed->tmp.filter;
change_values<-function(myd,column,values){
	which(is.na(myd[,column]))->na_index;
	if(length(na_index)>0){
		myd[na_index,]->myd_na;
		"Un"->myd_na[is.na(myd_na[,column]),column];
	}else{
		NULL->myd_na;
	}
	myd[!is.na(myd[,column]),]->myd
	colnames(myd)[column]->column_name;
	table(myd[,column])->N_stage.table;
	as.factor(names(N_stage.table))->N_stage.names;
	data.frame(column_name=N_stage.names,"Value"=values)->N_stage.df;
	c()->tmp.value;
	for(i in 1:nrow(myd)){
		for(j in 1:nrow(N_stage.df)){
			if(myd[i,column]==N_stage.df[j,1]){
				c(tmp.value,as.character(N_stage.df[j,2]))->tmp.value;
			}
		}
	}
	tmp.value->myd[,column];
	if(!is.null(myd_na)){
		rbind(myd,myd_na)->myd;
	}
	return(myd);
}
change_values(tmp.filter,12,c("Primary","DistantMetastasis","LocoregionalDisease","LocoregionalRecurrence","NewPrimaryTumor"))->tmp.filter;#new events
#------------------------------------------------------------------
split_factor<-function(myd,column,values){
	c()->cut_values;
	c()->range_name;
	for(k in 1:(length(values)-1)){
		c(values[k],values[k+1])->row.value;
		rbind(cut_values,row.value)->cut_values;
		c(range_name,paste(values[k],values[k+1],sep="~"))->range_name;
	}
	c("Start","End")->colnames(cut_values);
	as.data.frame(cut_values)->cut_values;
	range_name->cut_values$Name;
	c()->test.values;
	for(j in 1:nrow(myd)){
		for(i in 1:nrow(cut_values)){
			if(myd[j,column]>=cut_values[i,1] && myd[j,column]<cut_values[i,2]){
				c(test.values,cut_values[i,3])->test.values;
			}
		}
	}
	test.values->myd[,column]
	return(myd);
}
0->tmp.filter[is.na(tmp.filter[,4]),4]
split_factor(tmp.filter,4,c(0,50,60,70,80,100))->tmp.filter;

#------------------------------------------------------------------
par(mfrow=c(2,3))
for(i in c(3,4,5,6,12)){
	factors_barplot_ratio_clustering(tmp.filter,CNVCor_METCor_iC_cluster,i)
}
#-----------------------------------------------------------------------------------------------------------
Sample_info_summary<-function(myd,columns){
	list()->summary.res;
	for(i in columns){
		as.data.frame(table(myd[i]))->tmp.df;
		c(colnames(myd)[i],"Values")->colnames(tmp.df);
		print(tmp.df);
	}
}
merge(tmp.filter,CNVCor_METCor_iC_cluster,by.x="A0_Samples",by.y="A0_Samples")->tmp_filter_merged;
Sample_info_summary(tmp_filter_merged,c(4,5,6,9,12));
Sample_info_summary(subset_cluster(tmp_filter_merged,sub_clusters=c("iC1")),c(3,4,5,6,12));
#---merge C2+C3+C4:
change_values(CNVCor_METCor_iC_cluster,5,c("iC134","iC2","iC134","iC134"))->CNVCor_METCor_iC_cluster.merged;
par(mfrow=c(2,3))
for(i in c(3,4,5,6,12)){
	factors_barplot_ratio_clustering(tmp.filter,CNVCor_METCor_iC_cluster.merged,i)
}

list(GSE21050_exp_processed,tmp.filter,GSE71118_exp_processed)->datasets;
lapply(datasets,function(x){
	summary(x[which(x$Status==0),"A1_OS"])->x1;
	summary(x[which(x$Status!=0),"A1_OS"])->x2;
	list(x1,x2);
});
####################################################################################
###################################################################################
#--------------------------DESeq
library(DESeq2)
library(parallel)
read.table("drug-mining/HTSeq-Count_SARC.merge.SYMBOL.factors.txt",header=T,sep="\t",stringsAsFactors=F)->myd_ht_count
as.matrix(myd_ht_count[,-c(1:(g_start_column-1))])->myd_readcount_raw_matrix;
myd_ht_count$A0_Samples->rownames(myd_readcount_raw_matrix)
t(myd_readcount_raw_matrix)->myd_readcount_raw_matrix;
data.frame("SampleID"=CNVCor_METCor_iC_cluster$A0_Samples,"Condition"=CNVCor_METCor_iC_cluster$iCluster)->iC_group_condition;
#---
intersect(iC_group_condition$SampleID,myd_ht_count$A0_Samples)->ht_count_shared_samples;
iC_group_condition[which(iC_group_condition$SampleID%in%ht_count_shared_samples),]->iC_group_condition;
#-----
generate_group<-function(myd_clusters,target_clusters){
	c()->target_index;
	for(i in target_clusters){
		
		which(myd_clusters$Condition==i)->tmp_index;
		c(target_index,tmp_index)->target_index;
	}
	myd_clusters[target_index,]->myd_clusters_subset;
	data.frame("SampleID"=myd_clusters_subset[,1],"Condition"=as.character(myd_clusters_subset[,2]))->myd_clusters_subset;
	return(myd_clusters_subset);
}
generate_group(iC_group_condition,c("iC2","iC3"))->iC_filter_group;
myd_readcount_raw_matrix[,as.character(iC_filter_group$SampleID)]->myd_readcount_matrix;##log2(iC3/iC2)
detectCores()->no_cores;
makeCluster(no_cores-1)->c1;	
clusterExport(c1,c("myd_readcount_matrix"));#
parSapply(c1,1:nrow(myd_readcount_matrix),function(i){length(myd_readcount_matrix[i,myd_readcount_matrix[i,]<5])->failed_samples;if(failed_samples/ncol(myd_readcount_matrix)<0.5){i}})->filter_colums;
stopCluster(c1);
unlist(filter_colums)->filter_colums;
myd_readcount_matrix[filter_colums,]->myd_readcount_matrix;
#---------------------------------do DESeq
do_DESeq<-function(expd,colData){
	DESeqDataSetFromMatrix(countData=expd,colData=colData,design=~Condition)->dds
	keep<-rowMeans(counts(dds))>=5
	dds[keep,]->dds.keep
	DESeq(dds.keep)->dds.keep.deseq
	results(dds.keep.deseq)->dds_res
	dds_res[order(dds_res$pvalue),]->dds_res_order;
	as.data.frame(dds_res_order)->dds_res_order_df;
	data.frame("TransID"=rownames(dds_res_order_df),dds_res_order_df)->dds_res_order_df;##log2((iC2 exp)/(iC1 exp))
	return(dds_res_order_df);
}
sample_distance_estimate<-function(expd,colData){
	DESeqDataSetFromMatrix(countData=expd,colData=colData,design=~Condition)->dds
	keep<-rowMeans(counts(dds))>=5
	dds[keep,]->dds.keep
	assay(vst(dds.keep),blind = FALSE)->dds_keep_vst;
	dist(t(dds_keep_vst))->sampleDists
	as.matrix(sampleDists)->sampleDistMatrix
	colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)->colors
	data.frame("Condition"=colData$Condition)->annot_col;
	colData$SampleID->rownames(annot_col)
	brewer.pal(9,"Set1")[c(1,2)]->annot_colors;
	names(table(annot_col))->names(annot_colors)
	list("Condition"=annot_colors)->annot_colors;
	pheatmap(sampleDistMatrix,clustering_distance_rows = sampleDists,clustering_distance_cols = sampleDists,col = colors,show_colnames=F,show_rownames=F,annotation_col=annot_col,annotation_row=annot_col,annotation_names_row=F,annotation_names_col=F,annotation_colors=annot_colors)->p;
	return(p);
}
do_DESeq(myd_readcount_matrix,iC_filter_group)->C3_C2_deseq.res;
sample_distance_estimate(myd_readcount_matrix,iC_filter_group)->myd_readcount_matrix_distance;
print(myd_readcount_matrix_distance)
write.table(C3_C2_deseq.res,"multi-omics/C3_C2_dds_res_order_df.txt",quote=F,sep="\t",row.names=F)
##########----------------------for C1~C2------------------
generate_group(iC_group_condition,c("iC1","iC2"))->iC_filter_group;
myd_readcount_raw_matrix[,as.character(iC_filter_group$SampleID)]->myd_readcount_matrix;##log2(iC2/iC1)
detectCores()->no_cores;
makeCluster(no_cores-2)->c1;	
clusterExport(c1,c("myd_readcount_matrix"));#
parSapply(c1,1:nrow(myd_readcount_matrix),function(i){length(myd_readcount_matrix[i,myd_readcount_matrix[i,]<5])->failed_samples;if(failed_samples/ncol(myd_readcount_matrix)<0.5){i}})->filter_colums;
stopCluster(c1);
unlist(filter_colums)->filter_colums;
myd_readcount_matrix[filter_colums,]->myd_readcount_matrix;
#---------------------------------do DESeq
do_DESeq(myd_readcount_matrix,iC_filter_group)->C2_C1_deseq.res;
sample_distance_estimate(myd_readcount_matrix,iC_filter_group)->myd_readcount_matrix_distance;
print(myd_readcount_matrix_distance)
write.table(C2_C1_deseq.res,"multi-omics/C2_C1_dds_res_order_df.txt",quote=F,sep="\t",row.names=F)
##########----------------------for C2~C4------------------
generate_group(iC_group_condition,c("iC2","iC4"))->iC_filter_group;
myd_readcount_raw_matrix[,as.character(iC_filter_group$SampleID)]->myd_readcount_matrix;##log2(iC4/iC2)
detectCores()->no_cores;
makeCluster(no_cores-2)->c1;	
clusterExport(c1,c("myd_readcount_matrix"));#
parSapply(c1,1:nrow(myd_readcount_matrix),function(i){length(myd_readcount_matrix[i,myd_readcount_matrix[i,]<5])->failed_samples;if(failed_samples/ncol(myd_readcount_matrix)<0.5){i}})->filter_colums;
stopCluster(c1);
unlist(filter_colums)->filter_colums;
myd_readcount_matrix[filter_colums,]->myd_readcount_matrix;
#---------------------------------do DESeq
do_DESeq(myd_readcount_matrix,iC_filter_group)->C4_C2_deseq.res;
sample_distance_estimate(myd_readcount_matrix,iC_filter_group)->myd_readcount_matrix_distance;
print(myd_readcount_matrix_distance)
write.table(C4_C2_deseq.res,"multi-omics/C4_C2_dds_res_order_df.txt",quote=F,sep="\t",row.names=F)
########################################################################
#--------------------draw volcano figures for DEGs:
draw_volcano_deseq_figure<-function(deseq,cut_p,cut_lfc,p_title){
	brewer.pal(9,"Set1")->myd_colors;
	1->deseq[is.na(deseq$padj),"padj"];
	c()->filter_index;
	c()->deseq_color;
	for(i in 1:nrow(deseq)){
		if(deseq$log2FoldChange[i]<(-cut_lfc) && deseq$padj[i]<cut_p){
			c(filter_index,i)->filter_index;
			c(deseq_color,myd_colors[2])->deseq_color;
		}else if(deseq$log2FoldChange[i]>cut_lfc && deseq$padj[i]<cut_p){
			c(filter_index,i)->filter_index;
			c(deseq_color,myd_colors[1])->deseq_color;
		}else{
			c(deseq_color,myd_colors[9])->deseq_color
		}
	}
	#print(table(deseq_color));flush.console();
	plot(x=deseq$log2FoldChange,y=-log(deseq$padj),pch=20,cex=0.5,col=deseq_color,xlim=c(-5,5),xlab="Log2 fold change",ylab="-Log10(padj)")
	legend("topright",legend=c("Down","Up"),pch=20,cex=1,col=myd_colors[c(2,1)],title=p_title,border=NA)
	return(deseq[filter_index,]);
}
draw_volcano_limma_figure<-function(resd,cut_p,cut_lfc,p_title){
	brewer.pal(9,"Set1")->myd_colors;
	c()->filter_index;
	c()->deseq_color;
	for(i in 1:nrow(resd)){
		if(resd$logFC[i]<(-cut_lfc) && resd$adj.P.Val[i]<cut_p){
			c(filter_index,i)->filter_index;
			c(deseq_color,myd_colors[2])->deseq_color;
		}else if(resd$logFC[i]>cut_lfc && resd$adj.P.Val[i]<cut_p){
			c(filter_index,i)->filter_index;
			c(deseq_color,myd_colors[1])->deseq_color;
		}else{
			c(deseq_color,myd_colors[9])->deseq_color
		}
	}
	#print(table(deseq_color));flush.console();
	plot(x=resd$logFC,y=-log(resd$adj.P.Val),pch=20,cex=0.5,col=deseq_color,xlab="Log2 fold change",ylab="-Log10(padj)")
	legend("topright",legend=c("Down","Up"),pch=20,cex=1,col=myd_colors[c(2,1)],title=p_title,border=NA)
	return(resd[filter_index,]);
}
par(mfrow=c(2,2))
draw_volcano_deseq_figure(C2_C1_deseq.res,0.05,1,"C2/C1")->C2_C1_deseq.res_filter;
mtext("C2/C1")
draw_volcano_deseq_figure(C3_C2_deseq.res,0.05,1,"C3/C2")->C3_C2_deseq.res_filter
mtext("C3/C2")
draw_volcano_deseq_figure(C4_C2_deseq.res,0.05,1,"C4/C3")->C4_C2_deseq.res_filter;
mtext("C4/C3")
#-------------------------compare C2_C1 and C3_C2 , C4_C2 DGEs: keep shared DGEs
library(limma)#bioconductor package
generate_DEGs<-function(deseq_list,x){
	c()->tmp_res;
	for(y in deseq_list){
		"F"->res;
		for(j in 1:nrow(y)){	
			if(y[j,1]== as.character(x) && abs(y[j,3])>cut_lfc && y[j,7]<cut_p){
				"T"->res;
			}
		}
		c(tmp_res,res)->tmp_res;
	}
	return(tmp_res);
}
prepare_vennCount<-function(deseq_list,cut_p,cut_lfc){
	#unlist(lapply(deseq_list,function(x){x$TransID}))->total_transIDs;
	#unique(total_transIDs)->total_transIDs;
	#detectCores()->no_cores;
	#makeCluster(no_cores-1)->c1;	
	#clusterExport(c1,c("generate_DEGs","deseq_list","cut_p","cut_lfc"));#
	#parLapply(c1,total_transIDs,function(x)generate_DEGs(deseq_list,x))->res_list;
	#stopCluster(c1);
	#unlist(res_list)->res_list;
	c()->res_list;
	c()->res_group;
	1->index;
	for(y in deseq_list){
		y[y$padj<cut_p,]->y;
		y[abs(y$log2FoldChange)>cut_lfc,"TransID"]->y_transIDs;
		c(res_list,as.character(y_transIDs))->res_list;
		c(res_group,rep(names(deseq_list)[index],length(y_transIDs)))->res_group;
		index+1->index;
	}
	data.frame("ConditionGroup"=res_group,"TransID"=res_list)->res;
	return(res);
}
list("C2_C1_deseq"=C2_C1_deseq.res_filter,"C3_C2_deseq"=C3_C2_deseq.res_filter,"C4_C2_deseq"=C4_C2_deseq.res_filter)->deseq_list;
#save.image("work-20190424-TME-2.RData")#use linux service to do this work;
prepare_vennCount(deseq_list,0.05,1)->deseq_list.df;
table(deseq_list.df$TransID,deseq_list.df$ConditionGroup)->deseq_list.table
vennCounts(deseq_list.table)->deseq_list.vennCount;
vennDiagram(deseq_list.vennCount,names=c("C2/C1","C3/C2","C4/C2"),cex=1.5,lwd=1.5,circle.col=brewer.pal(9,"Set1"));
#------------------------compare: 
#---------------------------------------------------------------------------------
c()->deseq_transid_gnames;
for(i in 1:nrow(deseq_list.table)){
	if(sum(deseq_list.table[i,])==3){
		rownames(deseq_list.table)[i]->g;
		c(deseq_transid_gnames,g)->deseq_transid_gnames;
	}
}
intersect(deseq_transid_gnames,colnames(myd_exp_processed))->CNVCor_METCor_df_genes.shared;
###################################################################step5: GO and KEGG-------------------------
#-------------------------GO and KEGG
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot);#upsetplot,enrichMap...
library("ggplotify")
library("grid")
library("gridExtra")
draw_GO_KEGG_plot<-function(myd,showCategory){
	if(length(myd)>0){
		c(brewer.pal(9,"Set1"),brewer.pal(11,"Set3"),brewer.pal(11,"Spectral"))->myd.colors;
		ifelse(nrow(myd)<showCategory,nrow(myd),showCategory)->showCategory;
		myd[1:showCategory,]->myd;
		myd$Count->myd.count;
		myd.count->myd_count_color;
		if(max(myd.count)>30){
			round(myd.count/2)->myd_count_color;
		}
		-log10(myd$pvalue)->myd.pvalue;
		plot(x=myd.count,y=myd.pvalue,pch=19,cex=-log10(myd$pvalue),col=myd.colors[myd_count_color],xlab="Count",ylab="-log10(P value)",xlim=c(min(myd.count)-1,max(myd.count)+1),ylim=c(min(myd.pvalue)-0.2,max(myd.pvalue)+0.2));
		text(x=myd.count,y=myd.pvalue,labels=myd$Description,offset=0.5,pos=3,cex=0.8,adj=c(0,1))
	}else{
		plot.new();
	}
}
do_GO_KEGG<-function(genes,gtitle){
	as.character(genes)->GO_KEGG.selected_genes;
	groupGO(gene=GO_KEGG.selected_genes,OrgDb=org.Hs.eg.db,ont="BP",keyType='SYMBOL')->GO_KEGG.selected_genes.group;
	as.data.frame(GO_KEGG.selected_genes.group)->GO_KEGG.selected_genes.group;
	write.table(GO_KEGG.selected_genes.group,paste(gtitle,".groupGO.txt",sep=""),quote=F,row.names=F,sep="\t");
	enrichGO(gene=GO_KEGG.selected_genes,OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont="BP",pAdjustMethod="BH",pvalueCutoff=0.1)->GO_KEGG.selected_genes.GO;
	GO_KEGG.selected_genes.GO->GO.enrichment;
	simplify(GO_KEGG.selected_genes.GO)->GO_KEGG.selected_genes.GO;
	as.data.frame(GO_KEGG.selected_genes.GO)->GO_KEGG.selected_genes.GO
	write.table(GO_KEGG.selected_genes.GO,paste(gtitle,".GOEnrichment.txt",sep=""),quote=F,row.names=F,sep="\t")
	bitr(GO_KEGG.selected_genes,fromType="SYMBOL",toType=c("ENSEMBL","ENTREZID"),OrgDb=org.Hs.eg.db)->myd.gName.df;
	enrichKEGG(gene=myd.gName.df[,3],organism="hsa",keyType="ncbi-geneid",pvalueCutoff=0.5)->myd.gName.df.KEGG
	write.table(as.data.frame(myd.gName.df.KEGG),paste(gtitle,".KEGGEnrichment.txt",sep=""),quote=F,row.names=F,sep="\t")
	#if(nrow(GO_KEGG.selected_genes.GO)>0 && nrow(myd.gName.df.KEGG)>0){
	#	par(mfrow=c(1,2))
	#	draw_GO_KEGG_plot(GO_KEGG.selected_genes.GO,10)
	#	draw_GO_KEGG_plot(myd.gName.df.KEGG,10)
	#}else if(nrow(GO_KEGG.selected_genes.GO)>0 && nrow(myd.gName.df.KEGG)==0){
	#	par(mfrow=c(1,1))
	#	draw_GO_KEGG_plot(GO_KEGG.selected_genes.GO,10)
	#}else if(nrow(GO_KEGG.selected_genes.GO)==0 && nrow(myd.gName.df.KEGG)>0){
	#	par(mfrow=c(1,1))
	#	draw_GO_KEGG_plot(myd.gName.df.KEGG,10)
	#}
	return(list("enrichGO"=GO.enrichment,"enrichKEGG"=myd.gName.df.KEGG));
}
do_GO_KEGG(CNVCor_METCor_df_genes.shared,"multi-omics/CNVCor_METCor_df_genes.shared")->CNVCor_METCor_df_genes.shared_enrichment;
#upsetplot(sig1_enrichment[[1]])
#upsetplot(sig2_enrichment[[1]])
#upsetplot(sig3_enrichment[[1]])
#upsetplot(sig4_enrichment[[1]])
as.grob(dotplot(CNVCor_METCor_df_genes.shared_enrichment$enrichGO,showCategory=10,color="qvalue"))->GO_p1
as.grob(dotplot(CNVCor_METCor_df_genes.shared_enrichment$enrichKEGG,showCategory=8,color="qvalue"))->KEGG_p1
grid.arrange(GO_p1,KEGG_p1,nrow=2);

#----------------compare iC1 and iC3 CNV: CNVCor_METCor_iC_cluster+CNVCor_METCor_df_genes+myd.cnv_merged
which(colnames(myd.cnv_region2symbol)%in%CNVCor_METCor_df_genes.shared)->cnv_intervals_index;
myd.cnv_region2symbol[,cnv_intervals_index]->myd.cnv_merged_changed;
for(i in 1:length(cnv_intervals_index)){
	for(j in 1:nrow(myd.cnv_merged_changed)){
			if(myd.cnv_merged_changed[j,i]>0.3){
				1->myd.cnv_merged_changed[j,i];
			}else if(myd.cnv_merged_changed[j,i]<(-0.3)){
				(-1)->myd.cnv_merged_changed[j,i];
			}else{
				0->myd.cnv_merged_changed[j,i];
			}
	}
}
data.frame("Sample"=myd.cnv_region2symbol$Sample,myd.cnv_merged_changed)->cnv_genes_changed;
merge(cnv_genes_changed,CNVCor_METCor_iC_cluster,by.x="Sample",by.y="A0_Samples")->cnv_genes_changed_merged;
lapply(c("iC1","iC2","iC3","iC4"),function(x){
	which(CNVCor_METCor_iC_cluster$iCluster==x)->x_index;
	as.character(CNVCor_METCor_iC_cluster$A0_Samples[x_index])->x_samples;
	0->x_sum;
	for(i in x_samples){
		length(which(cnv_genes_changed[i,-1]!=0))->i_length;
		x_sum+i_length->x_sum;
	}
	x_sum/length(x_index);
})
#-------------------------------------------------
#----------------compare MET: CNVCor_METCor_iC_cluster+CNVCor_METCor_df_genes+myd.methy_t
myd.methy[which(myd.methy$Symbol%in%CNVCor_METCor_df_genes.shared),]$Probe->methy_probes;
myd.methy_t[,as.character(methy_probes)]->myd.methy_t_changed;
for(i in 1:ncol(myd.methy_t_changed)){
	for(j in 1:nrow(myd.methy_t_changed)){
			if(myd.methy_t_changed[j,i]>0.8){
				1->myd.methy_t_changed[j,i];##methylation hyper
			}else if(myd.methy_t_changed[j,i]<0.2){
				(-1)->myd.methy_t_changed[j,i];##methylation hypo
			}else{
				0->myd.methy_t_changed[j,i];##normal
			}
	}
}
data.frame("Sample"=rownames(myd.methy_t),myd.methy_t_changed)->methy_t_changed;
merge(methy_t_changed,CNVCor_METCor_iC_cluster,by.x="Sample",by.y="A0_Samples")->methy_t_changed_merged;
lapply(c("iC1","iC2","iC3","iC4"),function(x){
	which(CNVCor_METCor_iC_cluster$iCluster==x)->x_index;
	as.character(CNVCor_METCor_iC_cluster$A0_Samples[x_index])->x_samples;
	0->x_sum;
	for(i in x_samples){
		length(which(methy_t_changed[i,-1]!=0))->i_length;
		x_sum+i_length->x_sum;
	}
	x_sum/length(x_index);
})


compare_icluster_cnv_count(methy_t_changed_merged,2:1582,c("iC1","iC4"),1582)->CNVCor_METCor_df_genes_metCount;
myd.methy[which(myd.methy$Probe%in%CNVCor_METCor_df_genes_metCount$ID),]->myd.methy_METCor_filter;
#######################################################################################
generate_color<-function(x,myd_colors){
	names(table(x))->x_table;
	myd_colors[1:length(x_table)]->x_colors;
	x_table->names(x_colors);
	return(x_colors);
}
preorder_cnvcor_metcor_df_genes_heatmap<-function(myd,sample_clusters,clusters,geneOrder){
	myd[,as.character(geneOrder)]->myd_order;
	data.frame("Sample"=myd$Sample,myd_order)->myd_order;
	merge(myd_order,sample_clusters,by.x="Sample",by.y="A0_Samples")->myd_merged;
	c()->order_index;
	for(i in clusters){
		which(myd_merged$iCluster==i)->tmp_index;
		hclust(dist(myd_merged[tmp_index,as.character(geneOrder)]))->tmp_index_hclust;
		tmp_index[tmp_index_hclust$order]->tmp_index;
		c(order_index,tmp_index)->order_index;
	}
	myd_merged[order_index,]->myd_merged;
	myd_merged$Sample->rownames(myd_merged);
	return(myd_merged)
}
brewer.pal(11,"RdYlBu")[c(11,7,1)]
#--------for EXP:
iCluster_merged[,CNVCor_METCor_df_genes.shared]->CNVCor_METCor_df_genes.matrix;
data.frame("Sample"=iCluster_merged$A0_Samples,CNVCor_METCor_df_genes.matrix)->CNVCor_METCor_df_genes.matrix

list("HistologicalType"=generate_color(iC_cluster_immunescore_sort[,8],myd_colors),"CNVCor_C"=generate_color(iC_cluster_immunescore_sort[,9],myd_colors),"METCor_C"=generate_color(iC_cluster_immunescore_sort[,10],myd_colors),"iCluster"=generate_color(iC_cluster_immunescore_sort[,11],myd_colors),"Description"=generate_color(CNVCor_METCor_df_genes.shared.GO$Description,myd_colors))->col_colors;
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(100),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(150))->fill_colors;
preorder_cnvcor_metcor_df_genes_heatmap(CNVCor_METCor_df_genes.matrix,CNVCor_METCor_iC_cluster,c("iC1","iC2","iC3","iC4"),CNVCor_METCor_df_genes.shared)->myd_exp_processed.df_order;
generate_gaps(as.numeric(substr(myd_exp_processed.df_order$iCluster,3,3)))->gaps_col;
pheatmap(t(log2(myd_exp_processed.df_order[,CNVCor_METCor_df_genes.shared]+1)),color=fill_colors,cluster_rows=T,cluster_cols=F,gaps_col=gaps_col,show_colnames=F,show_rownames=F,cellwidth=2,border_color=NA,annotation_col=myd_exp_processed.df_order[,c("HistologicalType","iCluster")],annotation_colors=col_colors)->exp_heatmap_gtable
CNVCor_METCor_df_genes.shared[exp_heatmap_gtable$tree_row$order]->CNVCor_METCor_df_genes.shared_sort;
as.character(myd_exp_processed.df_order$Sample)->shared_samples_order;
#-----for top GO enrichment genes:
##------------------------------------------------------------------map gene to GO and KEGG map:
get_enrich_genes<-function(enrichd){
	unlist(lapply(enrichd$geneID,function(x){
		unlist(strsplit(as.character(x),split="/"))->x_genes;
		x_genes;
	}))->res_genes;
	unique(res_genes);
	return(res_genes);
}
get_des_genes<-function(enrichd,des_list){
	c()->tmp_genes;
	c()->tmp_des;
	as.character(enrichd$ID)->tmp_terms;
	enrichd$Description->names(tmp_terms);
	for(des in des_list){
		which(enrichd$Description%in%des)->des_index;
		get_enrich_genes(enrichd[des_index,])->des_genes;
		c(tmp_genes,des_genes)->tmp_genes;
		c(tmp_des,rep(des[1],length(des_genes)))->tmp_des;
	}
	unique(tmp_genes)->genes_uniq;
	unlist(lapply(genes_uniq,function(x){
		which(tmp_genes==x)->x_index;
		x_index[1];
	}))->tmp_genes_filter;
	data.frame("gName"=tmp_genes[tmp_genes_filter],"Description"=tmp_des[tmp_genes_filter])->res_df;
	tmp_terms[as.character(res_df$Description)]->res_df$ID;
	return(res_df);
}
map_gene_to_GO_KEGG<-function(enrichd,go_des_list,kegg_des_list){
	enrichd$enrichGO->go_res;
	enrichd$enrichKEGG->kegg_res;
	get_enrich_genes(go_res)->go_res_genes;
	get_enrich_genes(kegg_res)->kegg_res_genes;
	#-----
	get_des_genes(go_res,go_des_list)->go_res.df;
	c("gName","GO_des","GO_id")->colnames(go_res.df)
	get_des_genes(kegg_res,kegg_des_list)->kegg_res.df;
	c("gName","KEGG_des","KEGG_id")->colnames(kegg_res.df)
	print(kegg_res.df);flush.console();
	merge(go_res.df,kegg_res.df,by.x="gName",by.y="gName",all.x=T,all.y=T)->go_kegg_res.df;
	return(go_kegg_res.df);
}
list(c("positive regulation of cell adhesion"),c("mesenchyme development","mesenchymal cell differentiation"),c("regulation of inflammatory response"),c("negative regulation of cell growth"))->go_des_list;
list(c("Human T-cell leukemia virus 1 infection"),c("Prion diseases"))->kegg_des_list
#map_gene_to_GO_KEGG(CNVCor_METCor_df_genes.shared_enrichment,go_des_list,kegg_des_list)->tt
get_des_genes(CNVCor_METCor_df_genes.shared_enrichment$enrichGO,go_des_list)->CNVCor_METCor_df_genes.shared.GO
c(colorRampPalette(rev(brewer.pal(11,"RdYlBu")[7:11]))(100),colorRampPalette(rev(brewer.pal(11,"RdYlBu")[1:6]))(150))->fill_colors;
pheatmap(t(log2(myd_exp_processed.df_order[,as.character(CNVCor_METCor_df_genes.shared.GO$gName)]+1)),color=fill_colors,cluster_rows=T,cluster_cols=F,gaps_col=gaps_col,show_colnames=F,show_rownames=T,border_color=NA,annotation_col=myd_exp_processed.df_order[,c("HistologicalType","iCluster")],annotation_row=CNVCor_METCor_df_genes.shared.GO[,c("Description","Description")],annotation_colors=col_colors)
#-----------------------------------------------draw CNVCor_METCor_df_genes heatmap
intersect(CNVCor_METCor_df_genes.shared_sort,colnames(cnv_genes_changed))->CNVCor_df_genes.shared_sort
preorder_cnvcor_metcor_df_genes_heatmap(cnv_genes_changed,CNVCor_METCor_iC_cluster,c("iC1","iC2","iC3","iC4"),CNVCor_df_genes.shared_sort)->cnv_genes_changed_order;
cnv_genes_changed_order[shared_samples_order,]->cnv_genes_changed_order;
generate_gaps(as.numeric(substr(cnv_genes_changed_order$iCluster,3,3)))->gaps_col;
pheatmap(t(cnv_genes_changed_order[,CNVCor_df_genes.shared_sort]),color=c("#313695","#BABABA","#A50026"),cluster_rows=F,cluster_cols=F,gaps_col=gaps_col,show_colnames=F,cellwidth=2,show_rownames=F,border_color=NA,annotation_col=cnv_genes_changed_order[,c("HistologicalType","iCluster")],annotation_colors=col_colors,silent=T)->cnv_heatmap_gtable;
cnv_heatmap_gtable#,border_color=brewer.pal(9,"Greys")[2]
#--------for MET: 
unlist(lapply(CNVCor_METCor_df_genes.shared_sort,function(x){which(myd.methy_METCor_filter$Symbol==x)->res;res}))->probe_order;
as.character(myd.methy_METCor_filter[probe_order,"Probe"])->METCor_df_genes.shared_sort
preorder_cnvcor_metcor_df_genes_heatmap(methy_t_changed,CNVCor_METCor_iC_cluster,c("iC1","iC2","iC3","iC4"),METCor_df_genes.shared_sort)->met_genes_changed_order;
met_genes_changed_order[shared_samples_order,]->met_genes_changed_order;
pheatmap(t(met_genes_changed_order[,METCor_df_genes.shared_sort]),color=c("#313695","#BABABA","#A50026"),cluster_rows=F,cluster_cols=F,gaps_col=gaps_col,show_colnames=F,cellwidth=2,show_rownames=F,border_color=NA,annotation_col=met_genes_changed_order[,c("HistologicalType","iCluster")],annotation_colors=col_colors)->met_heatmap_gtable;

grid.arrange(cnv_heatmap_gtable$gtable,met_heatmap_gtable$gtable,exp_heatmap_gtable$gtable,nrow=3)
#--------------coxph:
which(colnames(myd_exp_processed_log2)%in%CNVCor_METCor_df_genes.shared)->myd_exp_processed.index;
cox_univariant_gene_regr(myd_exp_processed_log2,myd_exp_processed.index)->CNVCor_METCor_df_genes.shared.coxph_df;
CNVCor_METCor_df_genes.shared.coxph_df[CNVCor_METCor_df_genes.shared.coxph_df$Pvalue<0.05,"gName"]->myd_exp.filtered_gnames;

##########################################------------------------for GEO data:
#install.packages("BiocManager")
#BiocManager::install("illuminaHumanv2.db", version = "3.8")
library(hgu133plus2.db)
library(hgu133a.db);
library(annotate)
library("illuminaHumanv3.db")#for illuminaHumanWG-v3,Illumina HumanHT12v3 
library("illuminaHumanv2.db")#for Illumina human-6 v2.0
library(illuminaHumanv4.db)#for Illumina HumanHT-12 V4.0 expression beadchip
library(affy)
#source("https://bioconductor.org/biocLite.R")
#biocLite("illuminaHumanv3.db")
process_cel_data<-function(celpath,method){
	ReadAffy(celfile.path=celpath)->celdat;
	if(method=="mas5"){
		mas5(celdat)->celdat_normalized;
	}else if(method=="rma"){
		rma(celdat)->celdat_normalized;
	}
	exprs(celdat_normalized)->celdat_normalized_exp
	colnames(celdat_normalized_exp)->expd_names;
	unlist(lapply(expd_names,function(x){unlist(strsplit(x,split="\\."))->res;res[1]}))->expd_names;
	expd_names->colnames(celdat_normalized_exp);
	data.frame("DATA"=rownames(celdat_normalized_exp),celdat_normalized_exp)->celdat_normalized_exp;
	#scale(celdat_normalized_exp,center=rep(1000,ncol(celdat_normalized_exp)))->celdat_mas5_exp_scaled;
	return(celdat_normalized_exp);
}
map_illuminaHuman2SYMBOL<-function(expd,db){
	as.character(expd[,1])->expd_probes;
	if(db=="illuminaHumanv3SYMBOL"){
		unlist(as.list(illuminaHumanv3SYMBOL))->illuminaHumanv4SYMBOL_list;
	}else if(db=="illuminaHumanv4SYMBOL"){
		unlist(as.list(illuminaHumanv4SYMBOL))->illuminaHumanv4SYMBOL_list;
	}else{
		unlist(as.list(illuminaHumanv2SYMBOL))->illuminaHumanv4SYMBOL_list;
	}
	illuminaHumanv4SYMBOL_list[as.character(expd_probes)]->expd_probes_symbol;
	#select(illuminaHumanv4.db,keys = expd_probes,columns=c("SYMBOL"),keytype="PROBEID")->expd_probes_symbol;
	which(is.na(expd_probes_symbol))->na_index;
	expd_probes[na_index]->expd_probes_symbol[na_index]
	data.frame("Probe"=expd_probes,"gName"=expd_probes_symbol,expd[,-1])->expd_t;
	return(expd_t);
}
map_affy_probe2SYMBOL<-function(expd,db){
	as.character(expd$DATA)->expd_probes;	
	getSYMBOL(expd_probes,db)->expd_probes_substr_SYMBOL;
	which(is.na(expd_probes_substr_SYMBOL))->na_index;
	names(expd_probes_substr_SYMBOL)[na_index]->expd_probes_substr_SYMBOL[na_index]
	data.frame("Probe"=expd[,1],"gName"=expd_probes_substr_SYMBOL,expd[,-1])->expd;
	return(expd)
}
read_geo_data<-function(geof){
	fread(geof,header=T,sep="\t",stringsAsFactors=F,fill=T)->geof_factors;
	as.data.frame(geof_factors)->geof_factors;#log2 transformed
	grep("_at",colnames(geof_factors))->probe_index;
	if(length(probe_index)>0){
		geof_factors[,-probe_index]->geof_factors;
	}
	return(geof_factors);
}
#-----------------------------------------for GSE21050 data: 310 samples, Affymetrix Human Genome U133 Plus 2.0 Array
#-----GSE21050: metastasis time information
read_geo_data("drug-mining/GSE21050_exp_symbol.merged.factors.txt")->GSE21050_exp_factors;
preprocess_expd_v2(GSE21050_exp_factors,9,"no","no")->GSE21050_exp_processed;
#-----GSE71118: metastasis time information
read.table("GSE71118_family/MergeExpro_contrib1-GPL570.txt",header=T,sep="\t",stringsAsFactors=F)->GSE71118_exp;
map_affy_probe2SYMBOL(GSE71118_exp,"hgu133plus2.db")->GSE71118_exp_symbol;
write.table(GSE71118_exp_symbol,"multi-omics/GSE71118_exp_symbol.txt",quote=F,sep="\t",row.names=F)
#--
read_geo_data("multi-omics/GSE71118_exp_symbol.merged.factors.txt")->GSE71118_exp_factors;
preprocess_expd_v2(GSE71118_exp_factors,9,"No","no")->GSE71118_exp_processed;
#----draw metastasis free time vs histological type:
layout(matrix(c(1,2,2),nrow=1))
draw_boxplot_genes_by_factors_v2(GSE21050_exp_processed,c("A1_OS"),"HistologicalType",myd_colors,"yes")
draw_survial_curve_custom(GSE21050_exp_processed,4,1000,myd_colors)#4.485264e-05
#-
draw_boxplot_genes_by_factors_v2(GSE71118_exp_processed,c("A1_OS"),"HistologicalType",myd_colors,"yes")
draw_survial_curve_custom(GSE71118_exp_processed,5,1000,myd_colors)#2.344183e-05
Sample_info_summary(GSE71118_exp_processed,c(3,5));
#-
#----draw metastasis free time vs histological type:
layout(matrix(c(1,2,2),nrow=1))
draw_boxplot_genes_by_factors_v2(myd_exp_processed,c("A1_OS","NewEventTime"),"HistologicalType",myd_colors,"yes")
draw_survial_curve_custom(iCluster_merged,6,1000,myd_colors)#0.00368

######################################################################
#------------------------------------GSE21050 coxph:
GSE21050_exp_processed->expd_processed;
which(colnames(expd_processed)%in%CNVCor_METCor_df_genes.shared)->expd_processed.index;
cox_univariant_gene_regr(expd_processed,expd_processed.index)->GSE21050.coxph_df;
GSE21050.coxph_df[GSE21050.coxph_df$Pvalue<0.05,"gName"]->GSE21050.filtered_gnames;
#------------------------------------GSE71118 coxph:
GSE71118_exp_processed->expd_processed;
which(colnames(expd_processed)%in%CNVCor_METCor_df_genes.shared)->expd_processed.index;
cox_univariant_gene_regr(expd_processed,expd_processed.index)->GSE71118.coxph_df;
GSE71118.coxph_df[GSE71118.coxph_df$Pvalue<0.05,"gName"]->GSE71118.filtered_gnames;
#--------------------
intersect(intersect(myd_exp.filtered_gnames,GSE21050.filtered_gnames),GSE71118.filtered_gnames)
"APBB1IP"->g;#"APBB1IP"   "ACVRL1" "ENO1"
CNVCor_METCor_df_genes.shared.coxph_df[which(CNVCor_METCor_df_genes.shared.coxph_df$gName==g),]
GSE21050.coxph_df[which(GSE21050.coxph_df$gName==g),]
GSE71118.coxph_df[which(GSE71118.coxph_df$gName==g),]
###################################
#--------------------------------------------------------for survival analysis 
draw_gene_survial_curve_custom<-function(myd,g,bk){
	myd[myd[,g]!="",]->myd.rm;
	myd.rm[!is.na(myd.rm[,g]),]->myd.rm;
	median(myd[,g],rm.na=T)->g_median;
	unlist(lapply(myd[,g],function(x){if(x>g_median){"H"}else{"L"}}))->g_HL;
	g_HL->myd.rm$ExpLevel;
	if(length(myd.rm$A1_OS)>1){
		Surv(as.numeric(myd.rm$A1_OS),as.numeric(myd.rm$Status))->myd.surv;
	}else if(length(myd.rm$OS)>1){
		Surv(as.numeric(myd.rm$OS),as.numeric(myd.rm$Status))->myd.surv;
	}
	survfit(formula=myd.surv~ExpLevel,data=myd.rm)->myd.fit;
	survdiff(formula=myd.surv~ExpLevel,rho=0,data=myd.rm)->myd.diff;
	table(myd.rm$ExpLevel)->myd.table;
	max(myd.rm$A1_OS)+100->max_xlim;
	plot(myd.fit,col=brewer.pal(length(myd.table),"Set1"),xlab="Time(days)",ylab="Overall Survival(%)",lwd=2,axes=F,main=paste("Method Kaplan Meier(",g,")",sep=""),xlim=c(0,max_xlim));
	axis(side=1,at=seq(0,max_xlim,bk),labels=seq(0,max_xlim,bk));
	axis(side=2,at=seq(0,1,0.5),labels=seq(0,100,50));
	1-pchisq(myd.diff$chisq,df=length(myd.diff$n)-1)->pvalue;
	legend("topright",legend=paste(names(myd.table),paste("(n=",myd.table,")",sep="")),fill=brewer.pal(length(myd.table),"Set1"),bty="n");
	if(pvalue<1e-5){
		legend("bottomleft",legend="p < 1e-5",bty="n")
	}else{
		legend("bottomleft",legend=paste("p=",round(pvalue,5),sep=""),bty="n")
	}
	return(c(pvalue,myd.table));
}
draw_gene_survial_curve_custom_cut<-function(myd,g,bk,cuts,myd_colors){
	myd[myd[,g]!="",]->myd.rm;
	myd.rm[!is.na(myd.rm[,g]),]->myd.rm;
	as.numeric(myd.rm[,g])->myd.rm[,g];
	sort(myd.rm[,g],na.last=NA)->g_x;
	c(min(g_x,na.rm=T))->g_cuts;
	c("L0")->g_cuts_level;
	for(i in 1:(cuts-1)){
		round(length(g_x)*(i/cuts))->i_index;
		c(g_cuts,g_x[i_index])->g_cuts;
		c(g_cuts_level,paste("L",i,sep=""))->g_cuts_level;
	}
	c(g_cuts,max(g_x,na.rm=T))->g_cuts;
	c(g_cuts_level,paste("L",cuts,sep=""))->g_cuts_level;
	c()->g_HL;
	for(i_g in myd.rm[,g]){
		for(i in 1:(length(g_cuts)-1)){
			if(i_g>=g_cuts[i] && i_g<=g_cuts[i+1]){
				c(g_HL,g_cuts_level[i+1])->g_HL;
				break;
			}
		}
	}
	g_HL->myd.rm$ExpLevel;
	if(length(myd.rm$A1_OS)>1){
		Surv(as.numeric(myd.rm$A1_OS),as.numeric(myd.rm$Status))->myd.surv;
	}else if(length(myd.rm$OS)>1){
		Surv(as.numeric(myd.rm$OS),as.numeric(myd.rm$Status))->myd.surv;
	}
	survfit(formula=myd.surv~ExpLevel,data=myd.rm)->myd.fit;
	survdiff(formula=myd.surv~ExpLevel,rho=0,data=myd.rm)->myd.diff;
	table(myd.rm$ExpLevel)->myd.table;
	max(myd.rm$A1_OS)+100->max_xlim;
	plot(myd.fit,col=myd_colors,xlab="Time(days)",ylab="Overall Survival(%)",lwd=2,axes=F,main=paste("Method Kaplan Meier(",g,")",sep=""),xlim=c(0,max_xlim));
	axis(side=1,at=seq(0,max_xlim,bk),labels=seq(0,max_xlim,bk),pos=0);
	rug(x=seq(0,max_xlim,bk)+bk/2,ticksize=-0.01,side=1,pos=0);
	axis(side=2,at=seq(0,1,0.2),labels=seq(0,100,20),pos=0);
	rug(x=seq(0,0.9,0.2)+0.1,ticksize=-0.01,side=2,pos=0);
	1-pchisq(myd.diff$chisq,df=length(myd.diff$n)-1)->pvalue;
	legend("topright",legend=paste(names(myd.table),paste("(n=",myd.table,")",sep="")),fill=myd_colors,bty="n");
	if(pvalue<1e-5){
		text(x=bk/2,y=0.1,labels="p < 1e-5",bty="n",cex=1.1,pos=4,adj=0.5,font=3)
	}else{
		text(x=bk/2,y=0.1,labels=paste("p=",round(pvalue,5),sep=""),bty="n",cex=1.1,pos=4,adj=0.5,font=3)
	}
	return(c(pvalue,myd.table));
}
draw_genes_two_exp_point_v1<-function(expd1,g1,expd2,g2,log_trans){
	as.character(expd1$A0_Samples)->g1_samples;
	as.character(expd2$Sample)->g2_samples;
	intersect(g1_samples,g2_samples)->g1_g2_samples;
	#---------same order
	expd1[,g1]->g1_values;
	g1_samples->names(g1_values);
	g1_values[g1_g2_samples]->g1_values;
	expd2[,g2]->g2_values;
	g2_samples->names(g2_values);
	g2_values[g1_g2_samples]->g2_values;
	#-------------
	if(log_trans=="yes"){
		plot(log2(g1_values+1),log2(g2_values+1),xlab=paste(c("log2(",g1,")"),collapse=""),ylab=paste(c("log2(",g2,")"),collapse=""),pch=20,cex=0.5,bty="o");
	}else{
		plot(g1_values,g2_valuesxlab=paste(c("",g1),collapse=""),ylab=paste(c("",g2),collapse=""),pch=20,cex=0.5,bty="o");
	}
	#-------------
	if(g1==g2){
		paste(g1,"_1",sep="")->g1;
		paste(g2,"_2",sep="")->g2;
	}
	data.frame(g1_values,g2_values)->expd_g1_g2;
	c(g1,g2)->colnames(expd_g1_g2)
	print(head(expd_g1_g2));flush.console();
	expd_g1_g2[expd_g1_g2[,1]!=0,]->expd_g1_g2;
	expd_g1_g2[expd_g1_g2[,2]!=0,]->expd_g1_g2;
	#------------------
	cor.test(expd_g1_g2[,g1],expd_g1_g2[,g2],method = c("spearman"))->g1_g2_cortest;
	g1_g2_cortest$estimate->cortest_cor;
	g1_g2_cortest$p.value->cortest_p;
	#--------
	paste(g2,g1,sep="~")->g1_g2_formula;
	as.formula(g1_g2_formula)->g1_g2_formula;
	if(log_trans=="yes"){
		lm(g1_g2_formula,data=log2(expd_g1_g2+1))->g1_g2_lmfit;
	}else{
		lm(g1_g2_formula,data=expd_g1_g2)->g1_g2_lmfit;
	}
	abline(g1_g2_lmfit,lty=3,col="red",lwd=2)
	if(cortest_p < 0.001){
		"p < 0.001"->cortest_p
	}else{
		paste("p = ",round(cortest_p,3),sep="")->cortest_p;
	}
	paste("Rho = ",round(cortest_cor,3),sep="")->cortest_cor;
	paste(cortest_cor,cortest_p,sep=",")->leg;
	mtext(leg,cex=1.1,font=3);
}
#
data.frame("Sample"=rownames(myd.methy_t),myd.methy_t)->myd.methy_t_sample
par(mfrow=c(1,4))
draw_genes_two_exp_point_v1(myd_exp_processed,g,myd.cnv_region2symbol_logTrans,g,"yes")#
MET_EXP.cor[which(MET_EXP.cor$MetGene==g),]
draw_genes_two_exp_point_v1(myd_exp_processed,g,myd.methy_t_sample,"cg12159575","yes")
#
draw_boxplot_genes_by_factors_v2(iCluster_merged,g,"iCluster",myd_colors,"yes")
draw_gene_survial_curve_custom_cut(myd_exp_processed,g,1000,2,brewer.pal(11,"Spectral")[c(4,10,11)])->tmp_res;#HR=1.21,p<0.05
draw_gene_survial_curve_custom_cut(GSE21050_exp_processed,g,1000,2,brewer.pal(11,"Spectral")[c(4,10,11)])->tmp_res;#
draw_gene_survial_curve_custom_cut(GSE71118_exp_processed,g,1000,2,brewer.pal(11,"Spectral")[c(4,10,11)])->tmp_res;#
#------------
#---heatmap:
list("TCGA"=CNVCor_METCor_df_genes.shared.coxph_df,"GSE21050"=GSE21050.coxph_df,"GSE71118"=GSE71118.coxph_df)->coxph_list;
prepare_superheatmap<-function(coxphd,cutp){
	c()->datasets;
	list()->genes;
	list()->pvalues;
	c()->x_length;
	for(x in coxphd){
		x[x$Pvalue<cutp,"gName"]->x_genes;
		x[x$Pvalue<cutp,"HR"]->x_pvalues;
		c(x_length,length(x_genes))->x_length;
		c(genes,list(x_genes))->genes;
		c(pvalues,list(x_pvalues))->pvalues;
	}
	as.character(unique(unlist(genes)))->uniq_genes;
	names(coxphd)->coxphd_names;
	c()->uniq_pvalues;
	coxphd_names->names(genes);coxphd_names->names(pvalues);
	for(i in coxphd_names){
		genes[[i]]->i_genes;
		pvalues[[i]]->i_pvalues;
		setdiff(uniq_genes,i_genes)->i_diff_genes;
		c(i_pvalues,rep(0,length(i_diff_genes)))->i_pvalues;
		c(as.character(i_genes),i_diff_genes)->names(i_pvalues);
		i_pvalues[uniq_genes]->i_pvalues;
		c(uniq_pvalues,i_pvalues)->uniq_pvalues;
	}
	matrix(uniq_pvalues,ncol=length(uniq_genes),byrow=T)->res;
	uniq_genes->colnames(res);
	coxphd_names->rownames(res);
	return(res);
}
prepare_superheatmap(coxph_list,0.05)->coxph_list.merged_df;
#--
unlist(lapply(1:127,function(x){length(which(coxph_list.merged_df[,x]!=0))}))->coxph_list.count;
coxph_list.merged_df[,which(coxph_list.count>1)]->coxph_list.merged_df_filter;
superheat(coxph_list.merged_df_filter,X.text=round(coxph_list.merged_df_filter,3),X.text.angle=90,bottom.label.text.angle = 90,heat.pal = c("blue","white","red"))
###################################################################################################################
#-----------SNV:
read.table("multi-omics/SARC.mutect2.gene_vcf2matrix.res",header=T,sep="\t",stringsAsFactors=F)->gene_snv_merged##
CNVCor_METCor_iC_cluster->snv_nmf_3_cluster
gsub("-",".",snv_nmf_3_cluster$A0_Samples)->snv_nmf_3_cluster$A0_Samples;
snv_nmf_3_cluster$A0_Samples->rownames(snv_nmf_3_cluster);
intersect(colnames(gene_snv_merged),snv_nmf_3_cluster$A0_Samples)->exp_snv_shared_samples;
snv_nmf_3_cluster[exp_snv_shared_samples,]->snv_nmf_3_cluster;
cal_gene_snv_group_ftest<-function(snvd,groups){
	snvd->tmp_matrix;
	#--------------------------------
	c()->gName;
	c()->group1_Mut;
	c()->group1_Nor;
	c()->group2_Mut;
	c()->group2_Nor;
	c()->FisherP;
	#c()->padj;
	groups[groups$iCluster=="iC2",1]->group1;
	groups[groups$iCluster!="iC2",1]->group2;
	for(i in 1:nrow(tmp_matrix)){
		tmp_matrix[i,as.character(group1)]->group1_values;
		tmp_matrix[i,as.character(group2)]->group2_values;
		c(length(group1_values[group1_values!=0]),length(group1_values[group1_values==0]))->group1_count;
		c(length(group2_values[group2_values!=0]),length(group2_values[group2_values==0]))->group2_count;
		matrix(c(group1_count,group2_count),nrow=2,byrow=F)->group1_group2_table;
		#print(group1_group2_table);
		#flush.console();
		fisher.test(group1_group2_table)->group1_group2_table_test;
		#p.adjust(group1_group2_table_test$p.value,n=sum(group1_count,group2_count))->test_padj;
		c(FisherP,group1_group2_table_test$p.value)->FisherP;
		c(group1_Mut,group1_count[1])->group1_Mut;
		c(group1_Nor,group1_count[2])->group1_Nor;
		c(group2_Mut,group2_count[1])->group2_Mut;
		c(group2_Nor,group2_count[2])->group2_Nor;
		#c(padj,test_padj)->padj;
	}
	data.frame("gName"=tmp_matrix$gName,"FisherP"=FisherP,"Mut_C1"=group1_Mut,"Nor_C1"=group1_Nor,"Mut_C234"=group2_Mut,"Nor_C234"=group2_Nor)->res;
	p.adjust(res$FisherP,n=nrow(res))->res$padj;
	res[order(res$FisherP),]->res;
	return(res);
}
cal_gene_snv_group_ftest(gene_snv_merged,snv_nmf_3_cluster)->gene_snv_count;
write.table(gene_snv_count,"multi-omics/gene_snv_count.txt",quote=F,sep="\t",row.names=F)

#--top50 mutated genes:
gene_snv_count$Mut_C1+gene_snv_count$Mut_C234->gene_snv_count$MutCount;
gene_snv_count[order(gene_snv_count$MutCount,decreasing=T),]->gene_snv_count;
as.character(gene_snv_count[1:50,1])->genes_snv_count_top50;
unlist(lapply(as.character(genes_snv_count_top50),function(x){which(gene_snv_merged$gName==x)->i;i}))->tmp_index;
gene_snv_merged[tmp_index,]->gene_snv_merged_filter;
#----------
preoder_snv_gene_samples<-function(snvd,groups,clusters){
	data.frame("gName"=snvd$gName)->res;
	for(g in groups){
		unlist(lapply(g,function(x){which(clusters$iCluster==x)->i;i}))->tmp_index;
		gsub("-",".",as.character(clusters[tmp_index,1]))->g_samples
		snvd[,g_samples]->snvd_filter
		unlist(lapply(1:ncol(snvd_filter),function(x){which(snvd_filter[,x]!=0)->x_index;length(snvd_filter[x_index,x])}))->tmp_index_count;
		cbind(res,snvd_filter[,order(tmp_index_count)])->res;
	}
	res[,-1]->res;
	return(res);
}
preoder_snv_gene_samples(gene_snv_merged_filter,c("iC1","iC2","iC3","iC4"),iCluster_merged)->gene_snv_merged_filter_sort;
unlist(lapply(c("iC1","iC2","iC3","iC4"),function(x){which(iCluster_merged$iCluster==x)->i;i}))->tmp_index;
iCluster_merged[,c(1:(g_start_column-1))]->tmp_annotation_col;
iCluster_merged$iCluster->tmp_annotation_col$iCluster;
gsub("-","\\.",as.character(tmp_annotation_col$A0_Samples))->rownames(tmp_annotation_col)
#tmp_annotation_col[rownames(tmp_filter_merged)[tmp_index],]->tmp_annotation_col
generate_gaps(as.numeric(substr(tmp_annotation_col$iCluster,3,3)))->gaps_col;

for(i in 1:nrow(gene_snv_merged_filter_sort)){
	which(gene_snv_merged_filter_sort[i,]!=0)->tmp_index;
	1->gene_snv_merged_filter_sort[i,tmp_index];
}
list("HistologicalType"=generate_color(iCluster_merged$HistologicalType,myd_colors),"iCluster"=generate_color(iCluster_merged$iCluster,myd_colors),"A2_Event"=generate_color(tmp_annotation_col$A2_Event,brewer.pal(9,"Greys")[c(3,9)]),"A1_OS"=colorRampPalette(brewer.pal(9,"Purples"))(200))->col_colors;
pheatmap(gene_snv_merged_filter_sort,cluster_rows=F,cluster_cols=F,labels_row=genes_snv_count_top50,cellheight=7.5,color=brewer.pal(9,"Greys")[c(2,6)],gaps_col=gaps_col,show_colnames=F,annotation_col=tmp_annotation_col[,c("HistologicalType","iCluster","A2_Event","A1_OS")],annotation_colors=col_colors,border_color="white")->gene_snv.p1
#as.grob(~barplot(gene_snv_count_filter$Mut_C3,horiz=T,))->gene_snv.p2;
#grid.arrange(gene_snv.p1$gtable,gene_snv.p2,ncol=2)
#data.frame("Symbol"=c(as.character(gene_snv_count_filter$gName),as.character(gene_snv_count_filter$gName)),"MutCount"=c(gene_snv_count_filter$Mut_C3,gene_snv_count_filter$Mut_C12),"iCluster"=rep(c("C3","C12"),each=nrow(gene_snv_count_filter)),"Index"=c(seq(nrow(gene_snv_count_filter),1,-1),seq(nrow(gene_snv_count_filter),1,-1)))->gene_snv_count.df;
#ggplot(gene_snv_count.df,aes(x=Index,y=MutCount))->p
#p+geom_bar(aes(fill=IRGCluster),stat="identity")+coord_flip()
prepare_snv_cluster_barplot<-function(snvd,clusters,f){
	names(table(clusters[,f]))->f_names;
	c()->g_snv_count;
	for(fn in f_names){
		which(clusters[,f]==fn)->fn_index;
		as.character(clusters[fn_index,"A0_Samples"])->fn_samples;
		#which(colnames(snvd)%in%fn_samples)->snvd_fn_index;
		for(i in 1:nrow(snvd)){
			length(which(snvd[i,fn_samples]!=0))->i_length;#print(i);
			c(g_snv_count,i_length)->g_snv_count;#print(i_length);
		}
		#flush.console();
	}
	matrix(g_snv_count,nrow=length(f_names),byrow=T)->g_snv_count.table
	t(g_snv_count.table)->g_snv_count.table;
	f_names->colnames(g_snv_count.table)
	snvd$gName->rownames(g_snv_count.table)
	unlist(lapply(1:nrow(g_snv_count.table),function(j){
		sum(g_snv_count.table[j,]);
	}))->g_snv_count.sum;
	g_snv_count.table[order(g_snv_count.sum,decreasing=T),]->g_snv_count.table;
	return(g_snv_count.table);
}
draw_snv_cluster_barplot<-function(snvd,top_genes,myd_colors){
	snvd->g_snv_count.table;
	colnames(snvd)->f_names
	#----plot barplot:
	t(g_snv_count.table)->g_snv_count.table;
	g_snv_count.table[,1:top_genes]->g_snv_count_top.table;
	barplot(g_snv_count_top.table[,seq(top_genes,1,-1)],col=myd_colors[1:length(f_names)],beside=F,las=2,horiz=T,border=NA)
	legend("bottomright",legend=f_names,fill=myd_colors[1:length(f_names)],border=NA,bty='n')
}
prepare_snv_cluster_barplot(gene_snv_merged,snv_nmf_3_cluster,"iCluster")->snv_nmf_3_cluster.snv_sum;
layout(matrix(c(1,2,1,3),ncol=2,byrow=T))
draw_snv_cluster_barplot(snv_nmf_3_cluster.snv_sum,50,myd_colors)
#------------------------------------------TMB: 
calculate_immunescore_pvalue<-function(myd,columns,clusters,cluster_col){
	list()->c_groups;
	for(i in clusters){
		which(myd[,cluster_col]==i)->tmp_index;
		c(c_groups,list(tmp_name=tmp_index))->c_groups;
	}
	clusters->names(c_groups);
	c()->compare_clusters;
	c()->rank_test_pvalues;
	for(i in 1:(length(clusters)-1)){
		for(j in (i+1):length(clusters)){
			paste(clusters[i],clusters[j],sep="~")->tmp_name;
			c(compare_clusters,tmp_name)->compare_clusters;
			for(k in columns){
				as.numeric(unlist(c_groups[clusters[i]]))->c1_index;
				as.numeric(unlist(c_groups[clusters[j]]))->c2_index;
				wilcox.test(myd[c1_index,k],myd[c2_index,k])$p.value->tmp_pvalue;
				c(rank_test_pvalues,tmp_pvalue)->rank_test_pvalues;
			}
		}
	}

	length(c_groups)->c_groups_len;
	matrix(rank_test_pvalues,c_groups_len*(c_groups_len-1)/2,length(columns),byrow=T)->res_df;	
	#print(c_groups_len*(c_groups_len-1)/2);
	#flush.console();
	compare_clusters->rownames(res_df);
	colnames(myd)[columns]->colnames(res_df);
	return(res_df);
}
retrive_TCGA_immuneValues(CNVCor_METCor_iC_cluster$A0_Samples,c("SNV.Neoantigens","Silent.Mutation.Rate","Nonsilent.Mutation.Rate"))->sample_immune_landscape;
merge(sample_immune_landscape,CNVCor_METCor_iC_cluster,by.x="Sample",by.y="A0_Samples")->iC_cluster_immune_landscape;
#write.table(iC_cluster_immune_landscape,"multi-omics/iC_cluster_immune_landscape.txt",quote=F,sep="\t",row.names=F)
draw_boxplot_genes_by_factors_v2(iC_cluster_immune_landscape,c("Silent.Mutation.Rate","Nonsilent.Mutation.Rate","SNV.Neoantigens"),"iCluster",myd_colors,"yes")
#---for CNV/MetHyper
merge(myd.CNV_MET_abnormal_frequency,CNVCor_METCor_iC_cluster,by.x="Samples",by.y="A0_Samples")->myd.CNV_MET_abnormal_frequency.cluster;
draw_boxplot_genes_by_factors_v2(myd.CNV_MET_abnormal_frequency.cluster,c("Gain","Loss","MetHyper","MetHypo"),"iCluster",myd_colors,"yes")->x
cnv_met_count_summary<-function(freqd,genes,f,sub_clusters){
	c()->sum_values;
	for(sc in sub_clusters){
		which(freqd[,f]==sc)->sc_index;
		for(g in genes){
			sum(freqd[sc_index,g])->g_sum;
			c(sum_values,g_sum)->sum_values;
		}
	}
	matrix(sum_values,nrow=length(sub_clusters),byrow=T)->sum_values.matrix;
	genes->colnames(sum_values.matrix);
	sub_clusters->rownames(sum_values.matrix);
	return(sum_values.matrix);
}
cnv_met_count_summary(myd.CNV_MET_abnormal_frequency.cluster,c("Gain","Loss","MetHyper","MetHypo"),"iCluster",c("iC1","iC2","iC3","iC4"))->cnv_met_count.summary
#-----------count the SNVs of iC1,iC2,iC3,iC4: 
cal_gene_snv_group_ftest_v2<-function(snvd,groups,sub_group){
	snvd->tmp_matrix;
	#--------------------------------
	c()->gName;
	c()->group1_Mut;
	c()->group1_Nor;
	c()->group2_Mut;
	c()->group2_Nor;
	c()->FisherP;
	#c()->padj;
	groups[groups$iCluster==sub_group,1]->group1;
	groups[groups$iCluster!=sub_group,1]->group2;
	for(i in 1:nrow(tmp_matrix)){
		tmp_matrix[i,as.character(group1)]->group1_values;
		tmp_matrix[i,as.character(group2)]->group2_values;
		c(length(group1_values[group1_values!=0]),length(group1_values[group1_values==0]))->group1_count;
		c(length(group2_values[group2_values!=0]),length(group2_values[group2_values==0]))->group2_count;
		matrix(c(group1_count,group2_count),nrow=2,byrow=F)->group1_group2_table;
		#print(group1_group2_table);
		#flush.console();
		fisher.test(group1_group2_table)->group1_group2_table_test;
		#p.adjust(group1_group2_table_test$p.value,n=sum(group1_count,group2_count))->test_padj;
		c(FisherP,group1_group2_table_test$p.value)->FisherP;
		c(group1_Mut,group1_count[1])->group1_Mut;
		c(group1_Nor,group1_count[2])->group1_Nor;
		c(group2_Mut,group2_count[1])->group2_Mut;
		c(group2_Nor,group2_count[2])->group2_Nor;
		#c(padj,test_padj)->padj;
	}
	data.frame("gName"=tmp_matrix$gName,"FisherP"=FisherP,"Mut_C1"=group1_Mut,"Nor_C1"=group1_Nor,"Mut_C234"=group2_Mut,"Nor_C234"=group2_Nor)->res;
	p.adjust(res$FisherP,n=nrow(res))->res$padj;
	res[order(res$FisherP),]->res;
	return(res);
}
cal_gene_snv_group_ftest_v2(gene_snv_merged,snv_nmf_3_cluster,"iC1")->gene_snv_count;
gene_snv_count->x;
matrix(c(sum(x[,3]),sum(x[,4]),sum(x[,5]),sum(x[,6])),nrow=2,byrow=F)->x_matrix;
chisq.test(x_matrix)#p-value < 2.2e-16










