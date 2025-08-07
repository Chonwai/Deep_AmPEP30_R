require(seqinr)
library(keras)
library(kerasR)
kerasR::keras_init()
require(caret)
require(randomForest)
library(ROCR)
require(protr)
type8.cluster17=c('AT','C','DE','F','G','H','IV','K','L','M','N','P','Q','R','S','V','W')
type3a.cluster19=c('FA','P','G','S','T','D','E','Q','N','K','R','H','W','Y','M','L','I','V','C')
type12.cluster18=c('TVLI','M','F','W','Y','C','A','H','G','N','Q','P','R','K','S','T','D','E')
type7.cluster15=c('C','K','R','W','Y','A','FILV','M','D','E','Q','H','TP','GS','N')
type12.cluster17=c('TVLI','M','F','W','Y','C','A','H','G','N','Q','P','R','K','S','T','DE')
best_types=c('Type_8_17','Type_3A_19','Type_12_18','Type_7_15','Type_12_17')
final_paper2_test_path=list(po='/home/jieluyan/yan_ic50/seq/paper2_testing/lee_po_0.8',ne='/home/jieluyan/yan_ic50/seq/paper2_testing/lee_ne_0.8')
read_seq <- function(test_name='~/yan_ic50/seq/chris_seq.fasta',check_len=T){
	protdata = read.fasta(test_name,seqtype="AA",as.string=TRUE)
	name=getName(protdata)
	test_seqs = unlist(getSequence(protdata,as.string=TRUE),recursive=FALSE)
	test <- unlist(test_seqs)
	l=nchar(test)
	index=which(l %in% 5:30)
	if (check_len){
		if (length(index)==0){
			return(NULL)
		}else{
			return(list(seq=test[index],name=name[index]))
		}
	}else{
		return(list(seq=test,name=name))
	}
}
print_error_seqs <- function(path){
        protdata = read.fasta(path,seqtype="AA",as.string=TRUE)
        name=getName(protdata)
        test_seqs = unlist(getSequence(protdata,as.string=TRUE),recursive=FALSE)
        test <- unlist(test_seqs)
        l=nchar(test)
        index=which(l %in% 5:30)
	if (length(index)==length(test)){return(NULL)}else{return(list(seq=test[-index],name=name[-index]))}
}
read_pn <- function(active_name='/home/axpep/AxPEP_WWW_RUNS/Deep-AmPEP30/dataset/model_po.fasta',unactive_name='/home/axpep/AxPEP_WWW_RUNS/Deep-AmPEP30/dataset/model_ne.fasta',check_len=T){
    protdata = read.fasta(active_name,seqtype="AA",as.string=TRUE)
    po.name=getName(protdata)
    active_seqs = unlist(getSequence(protdata,as.string=TRUE),recursive=FALSE)
    po <- unlist(active_seqs)
    lp=nchar(po)
    p.index=which(lp %in% 5:30)
    protdata = read.fasta(unactive_name,seqtype="AA",as.string=TRUE)
    ne.name=getName(protdata)
    active_seqs = unlist(getSequence(protdata,as.string=TRUE),recursive=FALSE)
    ne <- unlist(active_seqs)
    ln=nchar(ne)
    n.index=which(ln %in% 5:30)
    if (check_len){
        return(list(po=po[p.index],ne=ne[n.index],nname=ne.name[n.index],pname=po.name[p.index]))
        }else{
        return(list(po=po,ne=ne,nname=ne.name,pname=po.name))
    }
}
gene_ftdata <- function(pn){
	po=pn$po
	ne=pn$ne
	pdata=gene_aac(list(seq=po,name=pn$pname),label=1)
	ndata=gene_aac(list(seq=ne,name=pn$nname),label=0)
	data<-rbind(pdata,ndata)
	data=data.frame(data)
	return(data)
}
gene_aac <- function(seq_name,label=NULL){
	ft=NULL
	seqs=seq_name$seq
	name=seq_name$name
	n=length(seqs)
	for (seq in seqs){
		t=extractAAC(seq)
		t=t*nchar(seq)
		ft=rbind(ft,t)
	}
	rownames(ft)=name
	if (!is.null(label)){
		class=rep(label,n)
		ft=cbind(ft,class)
	}
	ft=data.frame(ft)
	return(ft)
}
totally_psekraac_best5 <- function(test_path=NULL,test=F,po_ne.path=NULL,without_fasta_name=F,check_len=T){
	if (test) {
		if (without_fasta_name){
			seq=readLines(test_path)
			name=paste0('test_seq_',1:length(seq))
			seqs=list(seq=seq,name=name)
		}else{
			seqs=read_seq(test_path,check_len=check_len)
		}
		data.aac=gene_aac(seqs)
	} else{
		if (is.null(po_ne.path)){pn=read_pn(check_len=check_len)}else{pn=read_pn(po_ne.path$po,po_ne.path$ne,check_len=check_len)}
		data_aac=gene_ftdata(pn)
		data.aac=data_aac[,-21]
		class=data_aac$class
	}
	data.temp=NULL
	descriptors=c(type8.cluster17,type3a.cluster19,type12.cluster18,type7.cluster15,type12.cluster17)
	for (x in descriptors){	
		if (nchar(x)>1){
			char_multi=substring(x, seq(1,nchar(x),1), seq(1,nchar(x),1))
			t=apply(data.aac[,char_multi],1,sum)
			data.temp=cbind(data.temp,t)
		}else{
			t=data.aac[,x]
			data.temp=cbind(data.temp,t)
		}
	}
	n=length(descriptors)
	name=paste0(descriptors,'_',1:n)
	if (test) {
		data=data.temp
		whole_name=name
	} else {
		data=cbind(data.temp,class)
		whole_name=c(name,'class')
	}
	colnames(data)=whole_name
	data=data.frame(data)
	return(data)
}
# save_model_hdf5(mdl,path)
cnn_yan_predict <- function(test_path,mdl_path='../Deep-AmPEP30/AmPEP30-CNN.mdl',without_fasta_name=F,check_len=T){
	cnn.mdl=load_model_hdf5(mdl_path)
	test_data=totally_psekraac_best5(test_path,test=T,without_fasta_name=without_fasta_name,check_len=check_len)
	temp_data=test_data
	row_col=c(86,1)
	temp_data=data.matrix(temp_data)
	n_test=nrow(temp_data)
	test_x=array(temp_data,dim=c(n_test, row_col,1))
	sco=keras_predict(cnn.mdl, test_x)
	pre=rep(1,n_test)
	pre[which(sco<0.5)]=0
	error_seq=print_error_seqs(test_path)
	name.error=NULL
	sco.error=NULL
	if (!is.null(error_seq)){
		name.error=error_seq$name
		len_error=length(name.error)
		sco.error=rep(-1,len_error)
	}
	return(list(sco=c(sco,sco.error),pre=c(pre,sco.error),seq.name=c(rownames(test_data),name.error)))
}
rf_yan_predict <- function(test_path,mdl_path='/home/axpep/AxPEP_WWW_RUNS/Deep-AmPEP30/AmPEP30-RF-1200tree.mdl',without_fasta_name=F,check_len=T){
        rf.mdl=readRDS(mdl_path)
        test_data=totally_psekraac_best5(test_path,test=T,without_fasta_name=without_fasta_name,check_len=check_len)
	sco <- predict(rf.mdl, newdata = test_data, type = "vote")
	pre=rep(1,nrow(sco))
	pre[which(sco[,2]<0.5)]=0
        error_seq=print_error_seqs(test_path)
        name.error=NULL
        sco.error=NULL
        if (!is.null(error_seq)){
                name.error=error_seq$name
                len_error=length(name.error)
                sco.error=rep(-1,len_error)
        }
        #sco = c(sco,sco.error)
        #pre = c(pre,sco.error)
        #seq.name = c(rownames(test_data),name.error)
	return(list(sco=c(sco[,2],sco.error),pre=c(pre,sco.error),seq.name=c(rownames(test_data),name.error)))
}
predict_c_grabarta <- function(){
#	path='./C_glabarta_protein/remove_duplicated_sequence_5aato30aa_cgrabarate.fasta'
	path='./C_glabarta_protein/remove_first_char_M_and_duplicated_seq_with_model_dataset.fasta'
	rs=cnn_yan_predict(path)
	protdata=read.fasta(path,seqtype="AA",as.string=TRUE)
	seq=unlist(unlist(getSequence(protdata,as.string=TRUE),recursive=FALSE))
	seq=read_seq(path)
	prob=data.frame('seq_name'=rs$seq.name,'class'=rs$pre,'AMP_probablity'=rs$sco,'seq'=seq)
	prob=prob[order(prob$AMP_probablity,decreasing=T),]
	write.table(prob, "./rs_chris_seq_remove_initial_M_predict_by_ampep30-cnn.out", sep=",",col.names = T, row.names = F, quote = FALSE)
}
develop_cnn_mdl_yan <- function(rf=F){
	pn=read_pn()
	data=totally_psekraac_best5(pn)
	test_path='/home/axpep/AxPEP_WWW_RUNS/Deep-AmPEP30/dataset/test_ne.fasta'
	test_data=totally_psekraac_best5(test_path,test=T)
	if (!rf){
		cb=cnn_yan_binary(data,cnn_filters=c(128,128),kerner_size=c(3,3),dense_layer=c(10),batch=64,epoch=100,row_col=c(ncol(data)-1,1),pad='same',test_data=test_data)
		save_model_hdf5(cb$mdl,'/home/axpep/AxPEP_WWW_RUNS/Deep-AmPEP30/AmPEP30-CNN-time-test.mdl')
	}else{
		rf=rf_yan(data,test_data=data[1:10,],mtry=1,ntree=1200)
		saveRDS(rf$mdl,'/home/axpep/AxPEP_WWW_RUNS/Deep-AmPEP30/AmPEP30-RF-time-test.mdl')
	}
}
rf_yan <- function(data=data_d,data_rest=NULL,mtry=1,ntree=500,test_data=NULL,def=F){
        col=ncol(data)
        data$class=factor(data$class)
        if (!is.null(data_rest)) {data_rest$class=factor(data_rest$class)}
        flds <- createFolds(data$class, k = 10, list = TRUE, returnTrain = FALSE)
        pre=NULL
        ori=NULL
        if (is.null(test_data)){
                for (i in 1:10){
                        train=data[-flds[[i]],]
                        test=data[flds[[i]],]
                        if (def){
                                rf <- randomForest(class ~., rbind(train,data_rest),ntree=100,proximity=TRUE)
                        }else{
                                rf <- randomForest(class ~., rbind(train,data_rest),proximity=TRUE,mtry=mtry,ntree=ntree)
                        }
                        prediction=predict(rf,test[,-col],type='prob')
                        pre=c(pre,prediction[,2])
                        ori=c(ori,(test[,col]))
                }
                ori=ori-1
        } else{
                if (def){
                        rf <- randomForest(class ~., rbind(data,data_rest),ntree=100,proximity=TRUE)
                }else{
                        rf <- randomForest(class ~., rbind(data,data_rest),proximity=TRUE,mtry=mtry,ntree=ntree)
                }
                if (ncol(test_data)==col){
                        prediction=predict(rf,test_data[,-col],type='prob')
                        ori=test_data[,col]
                } else{
                        prediction=predict(rf,test_data,type='prob')
                }
                pre=prediction[,2]
        }
        return(list(mdl=rf,sco=pre,ori=ori))
}

cnn_yan_binary <- function(data,epoch=10,batch=128,r=0.2,dense_layer=c(100,50),cnn_filters=c(32,16),kerner_size=c(5,3),max_pool=NULL,row_col=NULL,stride=NULL,cv=10,test_data=NULL,pad='same',maxpool=T){
	col=ncol(data)
	if (is.null(row_col)){row_col=c(col-1,1)}
	data_y=array(data$class)
	data_x=data.matrix(data[,-col])
	if (is.null(test_data)){
		row=nrow(data)
		set.seed(1)
		flds <- createFolds(1:row, k = cv, list = TRUE, returnTrain = FALSE)
	}else{
		cv=1
	}
	pre=NULL
	ori=NULL
	for (i in 1:cv){
		if (is.null(test_data)){
			n_test=length(flds[[i]])
			n_train=nrow(data_x)-n_test
			test_x=array(data_x[flds[[i]],],dim=c(n_test, row_col,1))
			train_x=array(data_x[-flds[[i]],],dim=c(n_train, row_col,1))
			train_y=data_y[-flds[[i]]]
			test_y=data_y[flds[[i]]]
		} else {
			temp_data=test_data
			if (ncol(test_data)==col){
				temp_data=test_data[,-col]
				ori=test_data$class
			}
			temp_data=data.matrix(temp_data)
			n_test=nrow(temp_data)
			n_train=nrow(data_x)
			test_x=array(temp_data,dim=c(n_test, row_col,1))
			train_x=array(data_x,dim=c(n_train, row_col,1))
			train_y=data_y
		}
		mod <- Sequential()

		for (i in 1:length(cnn_filters)){
			if (row_col[2]==1){
				kerner=c(kerner_size[i],1)
				if (is.null(max_pool)){ pool=c(2,1) } else{pool=max_pool}
			}else{
				kerner=rep(kerner_size[i],2)
				if (is.null(max_pool)){ pool=c(2,2) } else{pool=max_pool}
			}
			if (i==1){ 
				mod$add(BatchNormalization(input_shape=c(row_col,1)))
				mod$add( Conv2D(cnn_filters[i], kerner, input_shape=c(row_col,1), activation="relu",padding=pad)) 
			} 
			else{ 
				mod$add(Conv2D(cnn_filters[i], kerner, activation="relu",padding=pad)) 
			}
			if (is.null(stride)){
				if (maxpool){mod$add(MaxPooling2D(pool_size=pool))}
			} else{ if (maxpool){mod$add(MaxPooling2D(pool_size=pool,strides=stride,padding=pad))} }
			if (maxpool){mod$add(Dropout(r))}
		}
		if (is.null(stride)){
			if (!maxpool){mod$add(MaxPooling2D(pool_size=pool,padding=pad))}
		} else{ if (!maxpool){mod$add(MaxPooling2D(pool_size=pool,strides=stride,padding=pad))} }
		if (!maxpool){mod$add(Dropout(r))}
		mod$add(Flatten())
		for (i in dense_layer){
			mod$add(Dense(i,activation="relu"))
			mod$add(Dropout(r))
		}
		mod$add(Dense(1))
		mod$add(Activation('sigmoid'))

		keras_compile(mod,  loss = 'binary_crossentropy', optimizer = RMSprop(lr = 0.00025))
		keras_fit(mod, train_x, train_y, batch_size = batch, epochs = epoch, verbose = 1)
		prediction=keras_predict(mod, test_x)
		#prediction=predict_on_batch(model, test_x)
		pre=c(pre,prediction)
		if (is.null(test_data)){ ori=c(ori,test_y)}	
	}
	return(list(mdl=mod,sco=pre,ori=ori))
}
cal_mat <- function(sco,ori,threshold=0.5,cal='auc'){
#library(ROCR)
#library(SDMTools)	
	conf_mat=confusion.matrix(ori,sco)
	tp=as.double(conf_mat[2,2])
	tn=as.double(conf_mat[1,1])
	fp=as.double(conf_mat[2,1])
	fn=as.double(conf_mat[1,2])
	acc=(tp+tn)/(tp+tn+fp+fn)
	sn=tp/(tp+fn)# ture positive rate, sensitivity or recall
	sp=tn/(tn+fp)# ture negative rate or specificity
	po=acc
	pe=((tp+fp)*(tp+fn)+(tn+fp)*(tn+fn))/(tp+tn+fp+fn)^2
	kappa=(po-pe)/(1-pe)
	mcc=(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
	if (cal=='auc'){
		mod=prediction(sco,ori)
		mod.auc=performance(mod,'auc')
        	auc.roc=mod.auc@y.values[[1]]
		pr.curve=performance(mod,'prec','rec')
		x=pr.curve@x.values[[1]]
		y=pr.curve@y.values[[1]]
		auc.pr=sum(diff(x)*y[-1])
		all_rs=c(acc,auc.roc,auc.pr,kappa,sn,sp,mcc)
		return(list(acc=acc,auc.roc=auc.roc,auc.pr=auc.pr,kappa=kappa,sn=sn,sp=sp,mcc=mcc,all_rs=unlist(all_rs)))
	}else{
		all_rs=c(acc,kappa,sn,sp,mcc)
		return(list(acc=acc,kappa=kappa,sn=sn,sp=sp,mcc=mcc,all_rs=unlist(all_rs)))
	}
}
# this codes are time test codes for table 5 in our manuscript ---- run time performance:

	# feature generateion run time based on seq_10000.fasta:
	#feature_generation_run_time_performance = totally_psekraac_best5(test_path='seq_10000.fasta',test=T,without_fasta_name=F,check_len=F)

        #training rf model based on 3298 sequences:
        #develop_cnn_mdl_yan(rf=T)

	#training cnn model based on 3298 sequences:
	#develop_cnn_mdl_yan(rf=F)

	#testing rf model based on seq_10000.fasta:
	#rs = rf_yan_predict('seq_10000.fasta',check_len=F)

	#testing rf model based on seq_10000.fasta:
	#rs = cnn_yan_predict('seq_10000.fasta',check_len=F)

# all codes above are main backend codes for our AxPEP server:
	args = commandArgs(trailingOnly=TRUE)
	seq.name = read_seq(test_name=args[1],check_len=T)
	if (is.null(seq.name)){
		seq.name = read_seq(test_name=args[1],check_len=F)
		name = seq.name$name
		pre = rep(-1,length(name))
		prob = data.frame('seq_name'= name, 'class' = pre, 'AMP_probablity' = pre)
	}else{
		rs=cnn_yan_predict(args[1],check_len=T)
		prob=data.frame('seq_name'=rs$seq.name,'class'=rs$pre,'AMP_probablity'=sprintf("%1.6f",rs$sco))
	}
	write.table(prob, args[2], sep=" ",col.names = F, row.names = F, quote = FALSE)

