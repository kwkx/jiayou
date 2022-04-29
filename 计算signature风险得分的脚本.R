###################下面开始CXCRscore的计算，是12417_96########################################
setwd("E:\\10分文章新训练集\\表达谱汇总")
cox_gene_list=read.table("4genelist.txt", header=T, sep="\t", check.names=F)
fix(cox_gene_list)
rawexpression=read.table("GSE12417_GPL96.txt", header=T, sep="\t", check.names=F)
fix(rawexpression)
intersect(rawexpression[,1],cox_gene_list[,1])
out=data.frame()
for(i in 1:nrow(cox_gene_list))
{
  num=which(rawexpression[,1]==cox_gene_list[i,1])
  out=rbind(out,rawexpression[num,])
  print(i)
  
}

fix(out)
tcgaOut=out
row.names(tcgaOut)=tcgaOut[,1]
tcgaOut=tcgaOut[,2:ncol(tcgaOut)]
fix(tcgaOut)
###下面开始添加临床信息并计算风险值
#colnames(tcgaOut)=substr(colnames(tcgaOut), 1, 12) 
tcgaOut=t(tcgaOut)
fix(tcgaOut)
Targetcli=read.table("GSE12417_96_cli.txt", header=T, sep="\t", check.names=F)
fix(Targetcli)
fustat=vector(length = nrow(tcgaOut))
futime=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  
  try(
    {
      num=which(Targetcli[,1]==rownames(tcgaOut)[i])
      fustat[i]=Targetcli[num,3]
      futime[i]=Targetcli[num,2]
    },silent = T
  )
  print(i)
}
tcgaOut=cbind(tcgaOut,fustat,futime)
fix(tcgaOut)
tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365
tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
tcgaOut=as.data.frame(tcgaOut)
#####下面开始添加riskScore
riskScore=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  sumvalue=0
  for(j in 1:4)
  {
    sumvalue=sumvalue+tcgaOut[i,j]*cox_gene_list[j,2]
  }
  riskScore[i]=sumvalue
  print(i)
  
}
median(riskScore)
############分组的中位数为42.30105
medianTrainRisk=42.30105
tcgaOut=cbind(tcgaOut,riskScore)

risk=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  if(tcgaOut[i,7]<= 42.30105)
  {
    risk[i]="low"
  }
  else
  {
    risk[i]="high"
  }
  print(i)
}
tcgaOut=cbind(tcgaOut,risk)
fix(tcgaOut)
tcgaOut[1,ncol(tcgaOut)]
data.class(tcgaOut)
write.table(tcgaOut,
            file="E:\\10分文章新训练集\\预后模型和列线图汇总\\CXCR模型\\09.survival\\riskTest.txt",
            sep="\t",
            quote = F,
            row.names=T)


##################下面开始CXCRscore的计算，是TCGA########################################
setwd("E:\\10分文章新训练集\\表达谱汇总")
cox_gene_list=read.table("4genelist.txt", header=T, sep="\t", check.names=F)
fix(cox_gene_list)
rawexpression=read.table("TCGA.txt", header=T, sep="\t", check.names=F)
fix(rawexpression)
intersect(rawexpression[,1],cox_gene_list[,1])
out=data.frame()
for(i in 1:nrow(cox_gene_list))
{
  num=which(rawexpression[,1]==cox_gene_list[i,1])
  out=rbind(out,rawexpression[num,])
  print(i)
  
}

fix(out)
tcgaOut=out
row.names(tcgaOut)=tcgaOut[,1]
tcgaOut=tcgaOut[,2:ncol(tcgaOut)]
fix(tcgaOut)
###下面开始添加临床信息并计算风险值
colnames(tcgaOut)=substr(colnames(tcgaOut), 1, 12) 
tcgaOut=t(tcgaOut)
fix(tcgaOut)
Targetcli=read.table("TCGA_cli.txt", header=T, sep="\t", check.names=F)
fix(Targetcli)
fustat=vector(length = nrow(tcgaOut))
futime=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  
  try(
    {
      num=which(Targetcli[,1]==rownames(tcgaOut)[i])
      fustat[i]=Targetcli[num,3]
      futime[i]=Targetcli[num,2]
    },silent = T
  )
  print(i)
}
tcgaOut=cbind(tcgaOut,fustat,futime)
fix(tcgaOut)
tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365
tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
tcgaOut=as.data.frame(tcgaOut)
#####下面开始添加riskScore
riskScore=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  sumvalue=0
  for(j in 1:4)
  {
    sumvalue=sumvalue+tcgaOut[i,j]*cox_gene_list[j,2]
  }
  riskScore[i]=sumvalue
  print(i)
  
}
median(riskScore)
############8.471892
medianTrainRisk=8.471892
tcgaOut=cbind(tcgaOut,riskScore)

risk=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  if(tcgaOut[i,7]<=8.471892)
  {
    risk[i]="low"
  }
  else
  {
    risk[i]="high"
  }
  print(i)
}
tcgaOut=cbind(tcgaOut,risk)
fix(tcgaOut)
tcgaOut[1,ncol(tcgaOut)]
data.class(tcgaOut)
write.table(tcgaOut,
            file="E:\\10分文章新训练集\\预后模型和列线图汇总\\CXCR模型\\09.survival\\riskTest.txt",
            sep="\t",
            quote = F,
            row.names=T)

##################下面开始CXCRscore的计算，是106291########################################
setwd("E:\\10分文章新训练集\\表达谱汇总")
cox_gene_list=read.table("4genelist.txt", header=T, sep="\t", check.names=F)
fix(cox_gene_list)
rawexpression=read.table("GSE106291.txt", header=T, sep="\t", check.names=F)
fix(rawexpression)
intersect(rawexpression[,1],cox_gene_list[,1])
out=data.frame()
for(i in 1:nrow(cox_gene_list))
{
  num=which(rawexpression[,1]==cox_gene_list[i,1])
  out=rbind(out,rawexpression[num,])
  print(i)
  
}

fix(out)
tcgaOut=out
row.names(tcgaOut)=tcgaOut[,1]
tcgaOut=tcgaOut[,2:ncol(tcgaOut)]
fix(tcgaOut)
###下面开始添加临床信息并计算风险值
#colnames(tcgaOut)=substr(colnames(tcgaOut), 1, 12) 
tcgaOut=t(tcgaOut)
fix(tcgaOut)
Targetcli=read.table("GSE106291_cli.txt", header=T, sep="\t", check.names=F)
fix(Targetcli)
fustat=vector(length = nrow(tcgaOut))
futime=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  
  try(
    {
      num=which(Targetcli[,1]==rownames(tcgaOut)[i])
      fustat[i]=Targetcli[num,3]
      futime[i]=Targetcli[num,2]
    },silent = T
  )
  print(i)
}
tcgaOut=cbind(tcgaOut,fustat,futime)
fix(tcgaOut)
tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365
tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
tcgaOut=as.data.frame(tcgaOut)
#####下面开始添加riskScore
riskScore=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  sumvalue=0
  for(j in 1:4)
  {
    sumvalue=sumvalue+tcgaOut[i,j]*cox_gene_list[j,2]
  }
  riskScore[i]=sumvalue
  print(i)
  
}
median(riskScore)
############0.3090131
medianTrainRisk=0.3090131
tcgaOut=cbind(tcgaOut,riskScore)

risk=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  if(tcgaOut[i,7]<=0.3090131)
  {
    risk[i]="low"
  }
  else
  {
    risk[i]="high"
  }
  print(i)
}
tcgaOut=cbind(tcgaOut,risk)
fix(tcgaOut)
tcgaOut[1,ncol(tcgaOut)]
data.class(tcgaOut)
write.table(tcgaOut,
            file="E:\\10分文章新训练集\\预后模型和列线图汇总\\CXCR模型\\09.survival\\riskTest.txt",
            sep="\t",
            quote = F,
            row.names=T)


##################下面开始CXCRscore的计算，是Beat########################################
setwd("E:\\10分文章新训练集\\表达谱汇总")
cox_gene_list=read.table("4genelist.txt", header=T, sep="\t", check.names=F)
fix(cox_gene_list)
rawexpression=read.table("beat.txt", header=T, sep="\t", check.names=F)
fix(rawexpression)
intersect(rawexpression[,1],cox_gene_list[,1])
out=data.frame()
for(i in 1:nrow(cox_gene_list))
{
  num=which(rawexpression[,1]==cox_gene_list[i,1])
  out=rbind(out,rawexpression[num,])
  print(i)
  
}

fix(out)
tcgaOut=out
row.names(tcgaOut)=tcgaOut[,1]
tcgaOut=tcgaOut[,2:ncol(tcgaOut)]
fix(tcgaOut)
###下面开始添加临床信息并计算风险值
#colnames(tcgaOut)=substr(colnames(tcgaOut), 1, 12) 
tcgaOut=t(tcgaOut)
fix(tcgaOut)
Targetcli=read.table("beat_cli.txt", header=T, sep="\t", check.names=F)
fix(Targetcli)
fustat=vector(length = nrow(tcgaOut))
futime=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  
  try(
    {
      num=which(Targetcli[,1]==rownames(tcgaOut)[i])
      fustat[i]=Targetcli[num,2]
      futime[i]=Targetcli[num,3]
    },silent = T
  )
  print(i)
}
tcgaOut=cbind(tcgaOut,fustat,futime)
fix(tcgaOut)
tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365
tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
tcgaOut=as.data.frame(tcgaOut)
#####下面开始添加riskScore
riskScore=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  sumvalue=0
  for(j in 1:4)
  {
    sumvalue=sumvalue+tcgaOut[i,j]*cox_gene_list[j,2]
  }
  riskScore[i]=sumvalue
  print(i)
  
}
median(riskScore)
############11.75191
medianTrainRisk=11.75191
tcgaOut=cbind(tcgaOut,riskScore)

risk=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  if(tcgaOut[i,7]<=11.75191)
  {
    risk[i]="low"
  }
  else
  {
    risk[i]="high"
  }
  print(i)
}
tcgaOut=cbind(tcgaOut,risk)
fix(tcgaOut)
tcgaOut[1,ncol(tcgaOut)]
data.class(tcgaOut)
write.table(tcgaOut,
            file="E:\\10分文章新训练集\\预后模型和列线图汇总\\CXCR模型\\09.survival\\riskTest.txt",
            sep="\t",
            quote = F,
            row.names=T)
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################


###################下面开始CIN25的计算，是12417_96########################################
setwd("E:\\10分文章新训练集\\表达谱汇总")
cox_gene_list=read.table("CIN25list.txt", header=T, sep="\t", check.names=F)
fix(cox_gene_list)
rawexpression=read.table("GSE12417_GPL96.txt", header=T, sep="\t", check.names=F)
fix(rawexpression)
intersect(rawexpression[,1],cox_gene_list[,1])
cox_gene_list=cox_gene_list[which(cox_gene_list[,1]%in%intersect(rawexpression[,1],cox_gene_list[,1])),]
out=data.frame()
for(i in 1:nrow(cox_gene_list))
{
  num=which(rawexpression[,1]==cox_gene_list[i,1])
  out=rbind(out,rawexpression[num,])
  print(i)
  
}

fix(out)
tcgaOut=out
row.names(tcgaOut)=tcgaOut[,1]
tcgaOut=tcgaOut[,2:ncol(tcgaOut)]
fix(tcgaOut)
###下面开始添加临床信息并计算风险值
#colnames(tcgaOut)=substr(colnames(tcgaOut), 1, 12) 
tcgaOut=t(tcgaOut)
fix(tcgaOut)
Targetcli=read.table("GSE12417_96_cli.txt", header=T, sep="\t", check.names=F)
fix(Targetcli)
fustat=vector(length = nrow(tcgaOut))
futime=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  
  try(
    {
      num=which(Targetcli[,1]==rownames(tcgaOut)[i])
      fustat[i]=Targetcli[num,3]
      futime[i]=Targetcli[num,2]
    },silent = T
  )
  print(i)
}
tcgaOut=cbind(tcgaOut,fustat,futime)
fix(tcgaOut)
tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365
tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
tcgaOut=as.data.frame(tcgaOut)
#####下面开始添加riskScore
riskScore=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  sumvalue=0
  for(j in 1:22)
  {
    sumvalue=sumvalue+tcgaOut[i,j]*cox_gene_list[j,2]
  }
  riskScore[i]=sumvalue
  print(i)
  
}
median(riskScore)
############分组的中位数为 212.2699
medianTrainRisk= 212.2699
tcgaOut=cbind(tcgaOut,riskScore)

risk=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  if(tcgaOut[i,25]<= 212.2699)
  {
    risk[i]="low"
  }
  else
  {
    risk[i]="high"
  }
  print(i)
}
tcgaOut=cbind(tcgaOut,risk)
fix(tcgaOut)
tcgaOut[1,ncol(tcgaOut)]
data.class(tcgaOut)
write.table(tcgaOut,
            file="E:\\10分文章新训练集\\预后模型和列线图汇总\\CIN25模型\\09.survival\\riskTest.txt",
            sep="\t",
            quote = F,
            row.names=T)


##################下面开始CIN25的计算，是TCGA########################################
setwd("E:\\10分文章新训练集\\表达谱汇总")
cox_gene_list=read.table("CIN25list.txt", header=T, sep="\t", check.names=F)
fix(cox_gene_list)
rawexpression=read.table("TCGA.txt", header=T, sep="\t", check.names=F)
fix(rawexpression)
intersect(rawexpression[,1],cox_gene_list[,1])
cox_gene_list=cox_gene_list[which(cox_gene_list[,1]%in%intersect(rawexpression[,1],cox_gene_list[,1])),]
out=data.frame()
for(i in 1:nrow(cox_gene_list))
{
  num=which(rawexpression[,1]==cox_gene_list[i,1])
  out=rbind(out,rawexpression[num,])
  print(i)
  
}

fix(out)
tcgaOut=out
row.names(tcgaOut)=tcgaOut[,1]
tcgaOut=tcgaOut[,2:ncol(tcgaOut)]
fix(tcgaOut)
###下面开始添加临床信息并计算风险值
colnames(tcgaOut)=substr(colnames(tcgaOut), 1, 12) 
tcgaOut=t(tcgaOut)
fix(tcgaOut)
Targetcli=read.table("TCGA_cli.txt", header=T, sep="\t", check.names=F)
fix(Targetcli)
fustat=vector(length = nrow(tcgaOut))
futime=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  
  try(
    {
      num=which(Targetcli[,1]==rownames(tcgaOut)[i])
      fustat[i]=Targetcli[num,3]
      futime[i]=Targetcli[num,2]
    },silent = T
  )
  print(i)
}
tcgaOut=cbind(tcgaOut,fustat,futime)
fix(tcgaOut)
tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365
tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
tcgaOut=as.data.frame(tcgaOut)
#####下面开始添加riskScore
riskScore=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  sumvalue=0
  for(j in 1:21)
  {
    sumvalue=sumvalue+tcgaOut[i,j]*cox_gene_list[j,2]
  }
  riskScore[i]=sumvalue
  print(i)
  
}
median(riskScore)
############分组的中位数为352.7999
medianTrainRisk=352.7999
tcgaOut=cbind(tcgaOut,riskScore)

risk=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  if(tcgaOut[i,24]<=352.7999)
  {
    risk[i]="low"
  }
  else
  {
    risk[i]="high"
  }
  print(i)
}
tcgaOut=cbind(tcgaOut,risk)
fix(tcgaOut)
tcgaOut[1,ncol(tcgaOut)]
data.class(tcgaOut)
write.table(tcgaOut,
            file="E:\\10分文章新训练集\\预后模型和列线图汇总\\CIN25模型\\09.survival\\riskTest.txt",
            sep="\t",
            quote = F,
            row.names=T)


#################下面开始CIN25的计算，是GSE10629#######################################
setwd("E:\\10分文章新训练集\\表达谱汇总")
cox_gene_list=read.table("CIN25list.txt", header=T, sep="\t", check.names=F)
fix(cox_gene_list)
rawexpression=read.table("GSE106291.txt", header=T, sep="\t", check.names=F)
fix(rawexpression)
intersect(rawexpression[,1],cox_gene_list[,1])
cox_gene_list=cox_gene_list[which(cox_gene_list[,1]%in%intersect(rawexpression[,1],cox_gene_list[,1])),]
out=data.frame()
for(i in 1:nrow(cox_gene_list))
{
  num=which(rawexpression[,1]==cox_gene_list[i,1])
  out=rbind(out,rawexpression[num,])
  print(i)
  
}

fix(out)
tcgaOut=out
row.names(tcgaOut)=tcgaOut[,1]
tcgaOut=tcgaOut[,2:ncol(tcgaOut)]
fix(tcgaOut)
###下面开始添加临床信息并计算风险值
#colnames(tcgaOut)=substr(colnames(tcgaOut), 1, 12) 
tcgaOut=t(tcgaOut)
fix(tcgaOut)
Targetcli=read.table("GSE106291_cli.txt", header=T, sep="\t", check.names=F)
fix(Targetcli)
fustat=vector(length = nrow(tcgaOut))
futime=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  
  try(
    {
      num=which(Targetcli[,1]==rownames(tcgaOut)[i])
      fustat[i]=Targetcli[num,3]
      futime[i]=Targetcli[num,2]
    },silent = T
  )
  print(i)
}
tcgaOut=cbind(tcgaOut,fustat,futime)
fix(tcgaOut)
tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365
tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
tcgaOut=as.data.frame(tcgaOut)
#####下面开始添加riskScore
riskScore=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  sumvalue=0
  for(j in 1:22)
  {
    sumvalue=sumvalue+tcgaOut[i,j]*cox_gene_list[j,2]
  }
  riskScore[i]=sumvalue
  print(i)
  
}
median(riskScore)
############分组的中位数为 1.060002
medianTrainRisk= 1.060002
tcgaOut=cbind(tcgaOut,riskScore)

risk=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  if(tcgaOut[i,25]<= 1.060002)
  {
    risk[i]="low"
  }
  else
  {
    risk[i]="high"
  }
  print(i)
}
tcgaOut=cbind(tcgaOut,risk)
fix(tcgaOut)
tcgaOut[1,ncol(tcgaOut)]
data.class(tcgaOut)
write.table(tcgaOut,
            file="E:\\10分文章新训练集\\预后模型和列线图汇总\\CIN25模型\\09.survival\\riskTest.txt",
            sep="\t",
            quote = F,
            row.names=T)
##################下面开始CIN25score的计算，是Beat########################################
setwd("E:\\10分文章新训练集\\表达谱汇总")
cox_gene_list=read.table("CIN25list.txt", header=T, sep="\t", check.names=F)
fix(cox_gene_list)
rawexpression=read.table("beat.txt", header=T, sep="\t", check.names=F)
fix(rawexpression)
intersect(rawexpression[,1],cox_gene_list[,1])
out=data.frame()
for(i in 1:nrow(cox_gene_list))
{
  num=which(rawexpression[,1]==cox_gene_list[i,1])
  out=rbind(out,rawexpression[num,])
  print(i)
  
}

fix(out)
tcgaOut=out
row.names(tcgaOut)=tcgaOut[,1]
tcgaOut=tcgaOut[,2:ncol(tcgaOut)]
fix(tcgaOut)
###下面开始添加临床信息并计算风险值
#colnames(tcgaOut)=substr(colnames(tcgaOut), 1, 12) 
tcgaOut=t(tcgaOut)
fix(tcgaOut)
Targetcli=read.table("beat_cli.txt", header=T, sep="\t", check.names=F)
fix(Targetcli)
fustat=vector(length = nrow(tcgaOut))
futime=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  
  try(
    {
      num=which(Targetcli[,1]==rownames(tcgaOut)[i])
      fustat[i]=Targetcli[num,2]
      futime[i]=Targetcli[num,3]
    },silent = T
  )
  print(i)
}
tcgaOut=cbind(tcgaOut,fustat,futime)
fix(tcgaOut)
tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365
tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
tcgaOut=as.data.frame(tcgaOut)
#####下面开始添加riskScore
riskScore=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  sumvalue=0
  for(j in 1:22)
  {
    sumvalue=sumvalue+tcgaOut[i,j]*cox_gene_list[j,2]
  }
  riskScore[i]=sumvalue
  print(i)
  
}
median(riskScore)
###########114.1442
medianTrainRisk=114.1442
tcgaOut=cbind(tcgaOut,riskScore)

risk=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  if(tcgaOut[i,25]<=114.1442)
  {
    risk[i]="low"
  }
  else
  {
    risk[i]="high"
  }
  print(i)
}
tcgaOut=cbind(tcgaOut,risk)
fix(tcgaOut)
tcgaOut[1,ncol(tcgaOut)]
data.class(tcgaOut)
write.table(tcgaOut,
            file="E:\\10分文章新训练集\\预后模型和列线图汇总\\CIN25模型\\09.survival\\riskTest.txt",
            sep="\t",
            quote = F,
            row.names=T)
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################
##########################################################################################################


###################下面开始hypoxia的计算，是12417_96########################################
setwd("E:\\10分文章新训练集\\表达谱汇总")
cox_gene_list=read.table("hypoxialist.txt", header=T, sep="\t", check.names=F)
fix(cox_gene_list)
rawexpression=read.table("GSE12417_GPL96.txt", header=T, sep="\t", check.names=F)
fix(rawexpression)
intersect(rawexpression[,1],cox_gene_list[,1])
cox_gene_list=cox_gene_list[which(cox_gene_list[,1]%in%intersect(rawexpression[,1],cox_gene_list[,1])),]
out=data.frame()
for(i in 1:nrow(cox_gene_list))
{
  num=which(rawexpression[,1]==cox_gene_list[i,1])
  out=rbind(out,rawexpression[num,])
  print(i)
  
}

fix(out)
tcgaOut=out
row.names(tcgaOut)=tcgaOut[,1]
tcgaOut=tcgaOut[,2:ncol(tcgaOut)]
fix(tcgaOut)
###下面开始添加临床信息并计算风险值
#colnames(tcgaOut)=substr(colnames(tcgaOut), 1, 12) 
tcgaOut=t(tcgaOut)
fix(tcgaOut)
Targetcli=read.table("GSE12417_96_cli.txt", header=T, sep="\t", check.names=F)
fix(Targetcli)
fustat=vector(length = nrow(tcgaOut))
futime=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  
  try(
    {
      num=which(Targetcli[,1]==rownames(tcgaOut)[i])
      fustat[i]=Targetcli[num,3]
      futime[i]=Targetcli[num,2]
    },silent = T
  )
  print(i)
}
tcgaOut=cbind(tcgaOut,fustat,futime)
fix(tcgaOut)
tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365
tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
tcgaOut=as.data.frame(tcgaOut)
#####下面开始添加riskScore
riskScore=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  sumvalue=0
  for(j in 1:4)
  {
    sumvalue=sumvalue+tcgaOut[i,j]*cox_gene_list[j,2]
  }
  riskScore[i]=sumvalue
  print(i)
  
}
median(riskScore)
############分组的中位数为 19.19479
medianTrainRisk= 19.19479
tcgaOut=cbind(tcgaOut,riskScore)

risk=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  if(tcgaOut[i,7]<= 19.19479)
  {
    risk[i]="low"
  }
  else
  {
    risk[i]="high"
  }
  print(i)
}
tcgaOut=cbind(tcgaOut,risk)
fix(tcgaOut)
tcgaOut[1,ncol(tcgaOut)]
data.class(tcgaOut)
write.table(tcgaOut,
            file="E:\\10分文章新训练集\\预后模型和列线图汇总\\hypoxia模型\\09.survival\\riskTest.txt",
            sep="\t",
            quote = F,
            row.names=T)


##################下面开始Hypoxia的计算，是TCGA########################################
setwd("E:\\10分文章新训练集\\表达谱汇总")
cox_gene_list=read.table("hypoxialist.txt", header=T, sep="\t", check.names=F)
fix(cox_gene_list)
rawexpression=read.table("TCGA.txt", header=T, sep="\t", check.names=F)
fix(rawexpression)
intersect(rawexpression[,1],cox_gene_list[,1])
cox_gene_list=cox_gene_list[which(cox_gene_list[,1]%in%intersect(rawexpression[,1],cox_gene_list[,1])),]
out=data.frame()
for(i in 1:nrow(cox_gene_list))
{
  num=which(rawexpression[,1]==cox_gene_list[i,1])
  out=rbind(out,rawexpression[num,])
  print(i)
  
}

fix(out)
tcgaOut=out
row.names(tcgaOut)=tcgaOut[,1]
tcgaOut=tcgaOut[,2:ncol(tcgaOut)]
fix(tcgaOut)
###下面开始添加临床信息并计算风险值
colnames(tcgaOut)=substr(colnames(tcgaOut), 1, 12) 
tcgaOut=t(tcgaOut)
fix(tcgaOut)
Targetcli=read.table("TCGA_cli.txt", header=T, sep="\t", check.names=F)
fix(Targetcli)
fustat=vector(length = nrow(tcgaOut))
futime=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  
  try(
    {
      num=which(Targetcli[,1]==rownames(tcgaOut)[i])
      fustat[i]=Targetcli[num,3]
      futime[i]=Targetcli[num,2]
    },silent = T
  )
  print(i)
}
tcgaOut=cbind(tcgaOut,fustat,futime)
fix(tcgaOut)
tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365
tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
tcgaOut=as.data.frame(tcgaOut)
#####下面开始添加riskScore
riskScore=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  sumvalue=0
  for(j in 1:4)
  {
    sumvalue=sumvalue+tcgaOut[i,j]*cox_gene_list[j,2]
  }
  riskScore[i]=sumvalue
  print(i)
  
}
median(riskScore)
############分组的中位数为12.97726
medianTrainRisk=12.97726
tcgaOut=cbind(tcgaOut,riskScore)

risk=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  if(tcgaOut[i,7]<=12.97726)
  {
    risk[i]="low"
  }
  else
  {
    risk[i]="high"
  }
  print(i)
}
tcgaOut=cbind(tcgaOut,risk)
fix(tcgaOut)
tcgaOut[1,ncol(tcgaOut)]
data.class(tcgaOut)
write.table(tcgaOut,
            file="E:\\10分文章新训练集\\预后模型和列线图汇总\\hypoxia模型\\09.survival\\riskTest.txt",
            sep="\t",
            quote = F,
            row.names=T)


#################下面开始Hypoxia的计算，是GSE10629#######################################
setwd("E:\\10分文章新训练集\\表达谱汇总")
cox_gene_list=read.table("hypoxialist.txt", header=T, sep="\t", check.names=F)
fix(cox_gene_list)
rawexpression=read.table("GSE106291.txt", header=T, sep="\t", check.names=F)
fix(rawexpression)
intersect(rawexpression[,1],cox_gene_list[,1])
cox_gene_list=cox_gene_list[which(cox_gene_list[,1]%in%intersect(rawexpression[,1],cox_gene_list[,1])),]
out=data.frame()
for(i in 1:nrow(cox_gene_list))
{
  num=which(rawexpression[,1]==cox_gene_list[i,1])
  out=rbind(out,rawexpression[num,])
  print(i)
  
}

fix(out)
tcgaOut=out
row.names(tcgaOut)=tcgaOut[,1]
tcgaOut=tcgaOut[,2:ncol(tcgaOut)]
fix(tcgaOut)
###下面开始添加临床信息并计算风险值
#colnames(tcgaOut)=substr(colnames(tcgaOut), 1, 12) 
tcgaOut=t(tcgaOut)
fix(tcgaOut)
Targetcli=read.table("GSE106291_cli.txt", header=T, sep="\t", check.names=F)
fix(Targetcli)
fustat=vector(length = nrow(tcgaOut))
futime=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  
  try(
    {
      num=which(Targetcli[,1]==rownames(tcgaOut)[i])
      fustat[i]=Targetcli[num,3]
      futime[i]=Targetcli[num,2]
    },silent = T
  )
  print(i)
}
tcgaOut=cbind(tcgaOut,fustat,futime)
fix(tcgaOut)
tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365
tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
tcgaOut=as.data.frame(tcgaOut)
#####下面开始添加riskScore
riskScore=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  sumvalue=0
  for(j in 1:4)
  {
    sumvalue=sumvalue+tcgaOut[i,j]*cox_gene_list[j,2]
  }
  riskScore[i]=sumvalue
  print(i)
  
}
median(riskScore)
############分组的中位数为 0.0006625542
medianTrainRisk= 0.0006625542
tcgaOut=cbind(tcgaOut,riskScore)

risk=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  if(tcgaOut[i,7]<= 0.0006625542)
  {
    risk[i]="low"
  }
  else
  {
    risk[i]="high"
  }
  print(i)
}
tcgaOut=cbind(tcgaOut,risk)
fix(tcgaOut)
tcgaOut[1,ncol(tcgaOut)]
data.class(tcgaOut)
write.table(tcgaOut,
            file="E:\\10分文章新训练集\\预后模型和列线图汇总\\hypoxia模型\\09.survival\\riskTest.txt",
            sep="\t",
            quote = F,
            row.names=T)
##################下面开始Hypoxia score的计算，是Beat########################################
setwd("E:\\10分文章新训练集\\表达谱汇总")
cox_gene_list=read.table("hypoxialist.txt", header=T, sep="\t", check.names=F)
fix(cox_gene_list)
rawexpression=read.table("beat.txt", header=T, sep="\t", check.names=F)
fix(rawexpression)
intersect(rawexpression[,1],cox_gene_list[,1])
out=data.frame()
for(i in 1:nrow(cox_gene_list))
{
  num=which(rawexpression[,1]==cox_gene_list[i,1])
  out=rbind(out,rawexpression[num,])
  print(i)
  
}

fix(out)
tcgaOut=out
row.names(tcgaOut)=tcgaOut[,1]
tcgaOut=tcgaOut[,2:ncol(tcgaOut)]
fix(tcgaOut)
###下面开始添加临床信息并计算风险值
#colnames(tcgaOut)=substr(colnames(tcgaOut), 1, 12) 
tcgaOut=t(tcgaOut)
fix(tcgaOut)
Targetcli=read.table("beat_cli.txt", header=T, sep="\t", check.names=F)
fix(Targetcli)
fustat=vector(length = nrow(tcgaOut))
futime=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  
  try(
    {
      num=which(Targetcli[,1]==rownames(tcgaOut)[i])
      fustat[i]=Targetcli[num,2]
      futime[i]=Targetcli[num,3]
    },silent = T
  )
  print(i)
}
tcgaOut=cbind(tcgaOut,fustat,futime)
fix(tcgaOut)
tcgaOut[,"futime"]=as.numeric(tcgaOut[,"futime"])/365
tcgaOut[,"fustat"]=as.numeric(tcgaOut[,"fustat"])
tcgaOut=as.data.frame(tcgaOut)
#####下面开始添加riskScore
riskScore=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  sumvalue=0
  for(j in 1:4)
  {
    sumvalue=sumvalue+tcgaOut[i,j]*cox_gene_list[j,2]
  }
  riskScore[i]=sumvalue
  print(i)
  
}
median(riskScore)
###########8.302912
medianTrainRisk=8.302912
tcgaOut=cbind(tcgaOut,riskScore)

risk=vector(length = nrow(tcgaOut))
for(i in 1:nrow(tcgaOut))
{
  if(tcgaOut[i,7]<=8.302912)
  {
    risk[i]="low"
  }
  else
  {
    risk[i]="high"
  }
  print(i)
}
tcgaOut=cbind(tcgaOut,risk)
fix(tcgaOut)
tcgaOut[1,ncol(tcgaOut)]
data.class(tcgaOut)
write.table(tcgaOut,
            file="E:\\10分文章新训练集\\预后模型和列线图汇总\\hypoxia模型\\09.survival\\riskTest.txt",
            sep="\t",
            quote = F,
            row.names=T)
