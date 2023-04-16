
library("pbmcapply")
library("msgps")


likelihood.pca=function(splt,dep,ind,min.leaf){
  one=dep[ind>=splt]
  two=dep[ind<splt]
  test1=length(one)
  test2=length(two)
  if((test1<min.leaf)|(test2<min.leaf)){
    return(Inf)
  }
  else{
    return(var(one)*(test1-1)+var(two)*(test2-1))
  }
}

likelihood.cat=function(splt,dep,ind,min.leaf){
  one=dep[ind==splt]
  two=dep[ind!=splt]
  test1=length(one)
  test2=length(two)
  if((test1<min.leaf)|(test2<min.leaf)){
    return(Inf)
  }
  else{
    return(var(one)*(test1-1)+var(two)*(test2-1))
  }
}

one.split=function(data,dep.name,ind.names.pca,ind.names.cat,min.leaf){
  optimal.splits=data.frame()
  dep=data[,dep.name]
  for(col in ind.names.cat){
    ind=data[,col]
    if(length(unique(ind))<2){next}
    temp=c()
    splits=unique(ind)
    for(i in splits){
        temp=c(temp,likelihood.cat(splt=i,dep=dep,ind=ind,min.leaf = min.leaf))
    }
    opt=which(temp==min(temp))
    rows=data.frame(Variable=col,splt.test=splits[opt],output=temp[opt])
    rows=rows[1,]
    optimal.splits=rbind(optimal.splits,rows)
  }
  for(col in ind.names.pca){
    ind=data[,col]
    ordered=sort(ind)
    ordered=unique(ordered)
    len=length(ordered)
    splits=((ordered[-1]+ordered[-len])/2)
    if(length(unique(splits))<2){next}
    temp=rep(Inf,len-1)
    for(i in 1:(len-1)){
      temp[i]=likelihood.pca(splt = splits[i],dep=dep,ind=ind,min.leaf=min.leaf)
    }
    opt=which(temp==min(temp))
    rows=data.frame(Variable=col,splt.test=splits[opt],output=temp[opt])
    rows=rows[1,]
    optimal.splits=rbind(optimal.splits,rows)
  }
  print(optimal.splits)
  if(nrow(optimal.splits)==0|min(optimal.splits$output)==Inf){
    return(0)
  }else if(min(optimal.splits$output)>(var(dep)*(length(dep)-1))){
    return("retry")
  }else{
    best.split=optimal.splits[which(optimal.splits$output==min(optimal.splits$output)),][1,]

    best.split=rbind(best.split,best.split)
    best.split=best.split[,c("Variable","splt.test")]
    if(best.split[1,"Variable"] %in% ind.names.cat){
      best.split["Direction"]=c("Equals","Not Equals")
    }else{
      best.split["Direction"]=c("Greater","Less")
    }
    return(best.split)
  }
}



rrtree=function(data,dep.name,ind.names.pca,independent.variables.categorical,num.of.ind.vars.to.sample,min.leaf){
  counter=1
  tree=data.frame("Variable"="No split","splt.test"=0,"Direction"="No split","Status"="Unknown","Parent"=0,"Node"=1)

  indices.per.leaf=list()
  indices.per.leaf[[counter]]=1:nrow(data)
  while(!all(tree$Status=="End")){
    iter=which(tree[,"Status"]!="End")
    for(i in iter){
      if(tree[i,"Status"]=="End"){next}
      temp.leaf=data[indices.per.leaf[[tree[i,"Node"]]],]
      tree[i,"Terminal"]=1
      ind.names=sample(c(ind.names.pca,independent.variables.categorical),num.of.ind.vars.to.sample,replace=FALSE)
      splt.criteria=one.split(data=temp.leaf,dep.name=dep.name,ind.names.pca=ind.names[(ind.names %in% ind.names.pca)],ind.names.cat=ind.names[(ind.names %in% independent.variables.categorical)],min.leaf=min.leaf)
      if(is.null(nrow(splt.criteria))){
        if(splt.criteria=="retry"){
          splt.criteria=one.split(data=temp.leaf,dep.name=dep.name,ind.names.pca=ind.names.pca[!(ind.names.pca %in% ind.names)],ind.names.cat=independent.variables.categorical[!(independent.variables.categorical %in% ind.names)],min.leaf=min.leaf)
        }
      }
      if(is.null(nrow(splt.criteria))){
        tree[i,"Status"]="End"
        tree[i,"Terminal"]=-1
        print("Leaf ended early 1")
        next
      }
      splt.criteria["Parent"]=tree[i,"Node"]
      splt.criteria["Status"]="Unknown"
      splt.criteria["Terminal"]=-1

      counter=counter+1
      splt.criteria[1,"Node"]=counter
      one.leaf=leaf.data(data=data,temp.tree=tree,leaf.criteria=splt.criteria[1,],indices.per.leaf=indices.per.leaf)
      indices.per.leaf[[counter]]=one.leaf

      counter=counter+1
      splt.criteria[2,"Node"]=counter
      two.leaf=leaf.data(data=data,temp.tree=tree,leaf.criteria=splt.criteria[2,],indices.per.leaf=indices.per.leaf)
      indices.per.leaf[[counter]]=two.leaf

      if(length(one.leaf)<min.leaf){
        splt.criteria[1,"Status"]="End"
      }
      if(length(two.leaf)<min.leaf){
        splt.criteria[2,"Status"]="End"
      }
      tree[i,"Status"]="End"
      tree=rbind(tree,splt.criteria)
      row.names(tree)=NULL
      print(tail(tree,2))
      print(length(one.leaf))
      print(length(two.leaf))
    }
  }
  tree.info=list(tree.table=tree,indices.per.leaf=indices.per.leaf)
  return(tree.info)
}

leaf.data=function(data, temp.tree, leaf.criteria,indices.per.leaf){
  node=leaf.criteria[1,"Parent"]
  if(node==0){
    return(indices.per.leaf[[1]])
  }else{
    aux=indices.per.leaf[[node]]
    dir=leaf.criteria[1,"Direction"]
    if(dir=="Greater"){
      temp.data=data[aux,leaf.criteria[1,"Variable"]]
      aux=aux[which(temp.data>=(leaf.criteria[1,"splt.test"]))]
    }
    if(dir=="Less"){
      temp.data=data[aux,leaf.criteria[1,"Variable"]]
      aux=aux[which(temp.data<(leaf.criteria[1,"splt.test"]))]
    }
    if(dir=="Equals"){
      temp.data=data[aux,leaf.criteria[1,"Variable"]]
      aux=aux[which(temp.data==leaf.criteria[1,"splt.test"])]
    }
    if(dir=="Not Equals"){
      temp.data=data[aux,leaf.criteria[1,"Variable"]]
      aux=aux[which(temp.data!=leaf.criteria[1,"splt.test"])]
    }
    return(aux)
  }
}



tree.coordinator=function(i,training.data,dependent.variables,independent.variables.pca,independent.variables.categorical,samples.per.tree,number.of.independent.variables.to.sample,cat.levels,indices,minimum.leaf.size){
  sample.row.numbers=sample(nrow(training.data),samples.per.tree,replace=TRUE)
  samp=training.data[sample.row.numbers,]
  min.leaf=minimum.leaf.size
  to.output=(rrtree(data=samp,dep.name=dependent.variables,ind.names.pca = independent.variables.pca,independent.variables.categorical=independent.variables.categorical,num.of.ind.vars.to.sample=number.of.independent.variables.to.sample,min.leaf = min.leaf))
  if(is.vector(indices)){
    to.output[["index"]]=indices[sample.row.numbers]
  }else{
    to.output[["index"]]=indices[sample.row.numbers,]
  }
  return(to.output)
}

out.of.bag.tree.indices=function(bagged.model){
  all.indices=bagged.model$data.as.matrix[,"index"]
  output=list()
  all.ind=1:length(bagged.model$trees)
  for(i in all.indices){
    indices=c()
    for(t in all.ind){
      samp=bagged.model$trees[[t]]$index[,"idx"]
      if(!(i %in% samp)){
        indices=c(indices,t)
      }
    }
    output[[i]]=indices
  }
  return(output)
}

build.tree=function(training.data,dependent.variable,independent.variables,samples.per.tree,minimum.leaf.size=2,number.of.independent.variables.to.sample,number.of.cores,number.of.trees,index.columns=NULL){
  tree.index=1:number.of.trees
  if(is.null(index.columns)){
    training.data[,"index"]=1:nrow(training.data)
    index.columns="index"
  }
  independent.variables.numerical=c()
  independent.variables.categorical=c()
  for(variable in independent.variables){
    if(is.numeric(training.data[,variable])){
      independent.variables.numerical=c(independent.variables.numerical,variable)
    }else{
      independent.variables.categorical=c(independent.variables.categorical,variable)
    }
  }
  print("Numerical variables:")
  print(independent.variables.numerical)
  print("Categorical variables:")
  print(independent.variables.categorical)
  cat.levels=list()
  train=training.data
  for(cat in independent.variables.categorical){
    cat.levels[[cat]]=unique(train[,cat])
    train[,cat]=as.numeric(factor(train[,cat],levels=cat.levels[[cat]]))
    train[is.na(train[,cat]),cat]=-1
  }
  train=as.matrix(train[,c(dependent.variable,independent.variables.numerical,independent.variables.categorical)])
  idx=1:nrow(train)
  train=cbind(idx,train)
  dimnames(train)[[2]][1]="index"
  index=training.data[,index.columns]
  index=cbind(idx,index)
  trees=pbmclapply(tree.index,tree.coordinator,mc.cores=number.of.cores,training.data=train,dependent.variables=dependent.variable,independent.variables.pca=independent.variables.numerical,independent.variables.categorical=independent.variables.categorical,samples.per.tree=samples.per.tree,number.of.independent.variables.to.sample=number.of.independent.variables.to.sample,indices=index,minimum.leaf.size=minimum.leaf.size)
  model=list(trees=trees,data.as.matrix=train,original.data=training.data,dependent.variable=dependent.variable,independent.variables.numerical=independent.variables.numerical,independent.variables.categorical=independent.variables.categorical,categorical.levels=cat.levels,index.columns=index.columns)
  oob.indices=out.of.bag.tree.indices(bagged.model = model)
  model[["oob.indices"]]=oob.indices
  return(model)
}


subset.from.case=function(tree.object,data,case){
  tree.from.fit=tree.object$tree.table
  current.node=tree.from.fit[1,]
  while(current.node[1,"Terminal"]!=-1){
    if(nrow(tree.from.fit)==1){
      return(data)
      break
    }
    node=current.node[1,"Node"]
    child=tree.from.fit[tree.from.fit$Parent==node,]
    if(child[1,"Direction"]=="Greater"){
      v=case[child[1,"Variable"]]
      if(v>=child[1,"splt.test"]){
        current.node=child[1,]
      }else{
        current.node=child[2,]
      }
    }
    if(current.node[1,"Terminal"]==-1){break}
    if(child[1,"Direction"]=="Equals"){
      v=case[child[1,"Variable"]]
      if(v==child[1,"splt.test"]){
        current.node=child[1,]
      }else{
        current.node=child[2,]
      }
    }
  }
  return(data[tree.object$index[,"idx"],][tree.object$indices.per.leaf[[current.node[1,"Node"]]],])
}



subset.from.ensemble=function(bagged.model,case,oob.indices){
  return=matrix(ncol=dim(bagged.model$data.as.matrix)[2])
  if(is.null(oob.indices)){
    for(i in 1:length(bagged.model$trees)){
      sub=subset.from.case(tree.object = bagged.model$trees[[i]],data=bagged.model$data.as.matrix,case=case)
      return=rbind(return,sub)
    }
    return(return[-1,])
  }else{
    for(i in oob.indices){
      sub=subset.from.case(tree.object = bagged.model$trees[[i]],data=bagged.model$data.as.matrix,case=case)
      return=rbind(return,sub)
    }
    return(return[-1,])
  }
}


predictor=function(training.data,newdata,categorical.variables,pca.variables,dependent.variable,model,predictions.only,categorical.levels){
  n.cat=c()
  d.cat=c()
  newdata=t(newdata)
  y=training.data[,dependent.variable]
  for(q in categorical.variables){
    if(newdata[,q] %in% training.data[,q]){
      n.cat=c(n.cat,q)
    }
  }
  if(length(n.cat)==1){
    cat=rbind(newdata[,n.cat],training.data[,n.cat])
    cat=data.frame(cat)
    cat[is.na(cat)]=0
    names(cat)[1]=n.cat
    cat[,1]=as.factor(cat[,1])
    levels(cat[,1])=categorical.levels[[n.cat]]
    f=as.formula(paste("~",paste(n.cat,collapse=" + ")))
    d.cat=model.matrix(f,cat)
    n.cat=d.cat[1,-1]
    xs=cbind(d.cat[-1,-1],training.data[,pca.variables])
  }else if(length(n.cat)==0){
    xs=training.data[,pca.variables]
  }else{
    cat=rbind(newdata[,n.cat],training.data[,n.cat])
    cat=data.frame(cat)
    cat[is.na(cat)]=0
    names(cat)=n.cat
    for(i in n.cat){
      cat[,i]=as.factor(cat[,i])
      levels(cat[,i])=categorical.levels[[i]]
    }
    f=as.formula(paste("~",paste(n.cat,collapse=" + ")))
    d.cat=model.matrix(f,cat)
    n.cat=d.cat[1,-1]
    xs=cbind(d.cat[-1,-1],training.data[,pca.variables])
  }
  new=as.matrix(c(1,n.cat,newdata[,pca.variables]))
  names(new)=c(n.cat,pca.variables)
  mod=msgps(X=xs,y=y,penalty = "enet",alpha=0)
  pred=t(new) %*% mod$dfgcv_result$coef
  if(predictions.only==TRUE){
    return(list(prediction=pred))
  }else{
    resid=cbind(rep(1,nrow(xs)),xs) %*% mod$dfgcv_result$coef - y
    rmse=mean((resid)^2)^(1/2)
    coef=data.frame(t(mod$dfgcv_result$coef))
    names(coef)=c("Intercept",dimnames(d.cat)[[2]][-1],pca.variables)
    return(list(prediction=pred,coefs=t(coef),rmse.error.of.lasso.model=rmse))
  }
}


prediction.coordinator=function(i,bagged.model,validation.data,labels=NULL,model,predictions.only,oob.indices=NULL){
  observation=validation.data[i,]
  if(!is.null(oob.indices)){
    indices=oob.indices[[i]]
  }else{indices=NULL}
  t.data=subset.from.ensemble(bagged.model = bagged.model, case = observation,oob.indices=indices)
  cat=c()
  for(z in bagged.model$independent.variables.categorical){
    if(length(unique(t.data[,z]))!=1){
      cat=c(cat,z)
    }
  }
  if(model=="Average"){
    ret=list(prediction=mean(t.data[,bagged.model$dependent.variable]))
  }else{
    ret=predictor(training.data=t.data,newdata=observation,categorical.variables=cat,pca.variables=bagged.model$independent.variables.numerical,dependent.variable=bagged.model$dependent.variable,model=model,predictions.only=predictions.only,categorical.levels=bagged.model$categorical.levels)
  }
  if(!is.null(labels)){
    for(p in names(labels)){
      ret[[names(labels)[p]]]=labels[i,p]
    }
  }
  return(ret)
}


predict.lm.from.rf=function(bagged.model,newdata,model="lasso",predictions.only=FALSE,number.of.cores){
  iter=1:nrow(newdata)
  cat.levels=bagged.model$categorical.levels
  new=newdata
  for(cat in bagged.model$independent.variables.categorical){
    new[,cat]=as.numeric(factor(new[,cat],levels=cat.levels[[cat]]))
    new[is.na(new[,cat]),cat]=-1
  }
  new=as.matrix(new[,c(bagged.model$dependent.variables,bagged.model$independent.variables.numerical,bagged.model$independent.variables.categorical)])
  if(all(bagged.model$index.columns %in% names(newdata))){
    indices=newdata[,bagged.model$index.columns]
  }else{indices=NULL}
  details=pbmclapply(iter,prediction.coordinator,mc.cores=number.of.cores,validation.data=new,bagged.model=bagged.model,labels=indices,model=model,predictions.only=predictions.only)
  return(details)
}

leaf.data.pruning=function(data, temp.tree, leaf.criteria){
  aux=1:nrow(data)
  if(leaf.criteria[1,"Parent"]==0){
    return(data)
    break
  }
  current.criteria=leaf.criteria
  rel.tree.subset=leaf.criteria
  while(current.criteria[1,"Parent"]!=0){
    current.criteria=temp.tree[temp.tree$Node==current.criteria[1,"Parent"],]
    rel.tree.subset=rbind(rel.tree.subset,current.criteria)
  }
  current.criteria=rel.tree.subset[rel.tree.subset$Parent==1,]
  while(nrow(current.criteria)!=0){
    node=current.criteria[1,"Parent"]
    dir=current.criteria[1,"Direction"]
    if(dir=="Greater"){
      temp.data=data[aux,current.criteria[1,"Variable"]]
      aux=aux[which(temp.data>=(current.criteria[1,"splt.test"]))]
    }
    if(dir=="Less"){
      temp.data=data[aux,current.criteria[1,"Variable"]]
      aux=aux[which(temp.data<(current.criteria[1,"splt.test"]))]
    }
    if(dir=="Equals"){
      temp.data=data[aux,current.criteria[1,"Variable"]]
      aux=aux[which(temp.data==current.criteria[1,"splt.test"])]
    }
    if(dir=="Not Equals"){
      temp.data=data[aux,current.criteria[1,"Variable"]]
      aux=aux[which(temp.data!=current.criteria[1,"splt.test"])]
    }
    current.criteria=rel.tree.subset[rel.tree.subset$Parent==current.criteria[1,"Node"],]
  }
  return(data[aux,])
}

prune.tree=function(tree.object,out.of.samp.data,sign.level,dep.name){
  tree=tree.object$tree.table
  tree["Optimized"]="Unknown"
  agg=aggregate(tree$Terminal,by=list(Parent=tree$Parent),FUN=sum)
  to.test=agg[agg$x==-2,"Parent"]
  while(!is.null(to.test)&!all(tree[which(tree$Parent %in% to.test),"Optimized"]=="End")){
    for(i in to.test){
      parent=tree[tree$Node==i,]
      children=tree[tree$Parent==i,]
      if(all(children$Optimized=="End")){next}
      one.leaf.os=leaf.data.pruning(data=out.of.samp.data,temp.tree=tree,leaf.criteria = children[1,])
      two.leaf.os=leaf.data.pruning(data=out.of.samp.data,temp.tree=tree,leaf.criteria = children[2,])

      if(is.null(nrow(one.leaf.os))|is.null(nrow(two.leaf.os))){
        tree=tree[tree$Parent!=i,]
        next
      }

      if((nrow(one.leaf.os)<2)|(nrow(two.leaf.os)<2)){
        tree=tree[tree$Parent!=i,]
        next
      }
      one=one.leaf.os[,dep.name]
      two=two.leaf.os[,dep.name]
      combined=c(one,two)


      if((var(one)*(length(one)-1)+var(two)*(length(two)-1))*sign.level<(var(combined)*(length(combined)-1))){
        tree[tree$Parent==i,"Optimized"]="End"
      }else{
        tree=tree[tree$Parent!=i,]
      }
    }
    for(t in unique(tree$Parent)){
      nodes=tree[tree$Parent==t,"Node"]
      if(nrow(tree[tree$Parent==nodes[1],])==0){
        tree[tree$Node==nodes[1],"Terminal"]=-1
      }
      if(nrow(tree[tree$Parent==nodes[2],])==0){
        tree[tree$Node==nodes[2],"Terminal"]=-1
      }
    }
    agg=aggregate(tree$Terminal,by=list(Parent=tree$Parent),FUN=sum)
    to.test=agg[agg$x==-2,"Parent"]
    print(to.test)
  }
  tree.info=list(tree.table=tree,indices.per.leaf=tree.object$indices.per.leaf,index=tree.object$index)
  print(nrow(tree))
  return(tree.info)
}

pruning.coordinator=function(i,bagged.model,sig.level){
  out.samp=bagged.model$data.as.matrix
  pruned.tree=prune.tree(tree.object=bagged.model$trees[[i]],out.of.samp.data=out.samp,sign.level=sig.level,dep.name=bagged.model$dependent.variable)
  return(pruned.tree)
}

pruning.parallel=function(bagged.model,sig.level,number.of.cores){
  iter=1:length(bagged.model$trees)
  p.tree=pbmclapply(iter,pruning.coordinator,mc.cores=number.of.cores,bagged.model=bagged.model,sig.level=sig.level)
  bagged.model$trees=p.tree
  return(bagged.model)
}

out.of.bag.error=function(bagged.model,model="lasso",number.of.cores){
  iter=bagged.model$data.as.matrix[,"index"]
  preds=pbmclapply(iter,prediction.coordinator,mc.cores=number.of.cores,validation.data=bagged.model$data.as.matrix,bagged.model=bagged.model,labels=NULL,model=model,predictions.only=TRUE,oob.indices=bagged.model$oob.indices)
  var=c()
  for(n in iter){
    var=c(var,(preds[[n]]$prediction-bagged.model$data.as.matrix[n,bagged.model$dependent.variable])^2)
  }
  rmse=mean(var)^(1/2)
  return(rmse)
}

auto.optimize=function(bagged.model,step.size=.025,guess=1.1,number.of.cores){
  error.recording=data.frame()
  for(sig in c(guess-step.size,guess,guess+step.size)){
    temp.model=pruning.parallel(bagged.model=bagged.model,sig.level=sig,number.of.cores=number.of.cores)
    tree.size=c()
    for(i in 1:length(temp.model$trees)){
      tree.size=c(tree.size,nrow(temp.model$trees[[i]]$tree.table))
    }
    size=mean(tree.size)
    error=out.of.bag.error(bagged.model=temp.model,model="lasso",number.of.cores=number.of.cores)
    to.add=data.frame("Significance"=sig,"RMSE"=error,"Tree Size"=size)
    error.recording=rbind(error.recording,to.add)
    print(error.recording)
  }
  last.error=Inf
  if(error.recording[2,"RMSE"]==min(error.recording[,"RMSE"])){
    print(paste("Optimal significance level: ",as.character(error.recording[2,"Significance"]),sep=""))
    return(pruning.parallel(bagged.model=bagged.model,sig.level=error.recording[2,"Significance"],number.of.cores=number.of.cores))
  }
  most.recent=min(error.recording[,"RMSE"])
  if(error.recording[1,"RMSE"]==min(error.recording[,"RMSE"])){
    direction="Down"
    sign.level=error.recording[1,"Significance"]
  }else if(error.recording[3,"RMSE"]==min(error.recording[,"RMSE"])){
    direction="Up"
    sign.level=error.recording[3,"Significance"]
  }else{break}

  while(last.error>most.recent){
    last.error=most.recent
    if(direction=="Down"){
      sign.level=sign.level-step.size
    }else{
      sign.level=sign.level+step.size
    }
    temp.model=pruning.parallel(bagged.model=bagged.model,sig.level=sign.level,number.of.cores=number.of.cores)
    tree.size=c()
    for(i in 1:length(temp.model$trees)){
      tree.size=c(tree.size,nrow(temp.model$trees[[i]]$tree.table))
    }
    size=mean(tree.size)
    error=out.of.bag.error(bagged.model=temp.model,model="lasso",number.of.cores=number.of.cores)
    to.add=data.frame("Significance"=sign.level,"RMSE"=error,"Tree Size"=size)
    error.recording=rbind(error.recording,to.add)
    print(error.recording)
    most.recent=error
    if(sign.level<=1){
      break
    }
  }
  error=out.of.bag.error(bagged.model=bagged.model,model="lasso",number.of.cores=number.of.cores)
  tree.size=c()
  for(i in 1:length(bagged.model$trees)){
    tree.size=c(tree.size,nrow(bagged.model$trees[[i]]$tree.table))
  }
  size=mean(tree.size)
  to.add=data.frame("Significance"=-1,"RMSE"=error,"Tree Size"=size)
  print(paste("Optimal significance level: ",error.recording[error.recording[,"RMSE"]==min(error.recording[,"RMSE"]),"Significance"],sep=""))
  return(pruning.parallel(bagged.model=bagged.model,sig.level=error.recording[error.recording[,"RMSE"]==min(error.recording[,"RMSE"]),"Significance"],number.of.cores=number.of.cores))
}


