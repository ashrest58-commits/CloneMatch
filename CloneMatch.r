library(ape)
library(ggtree)
library(phangorn)
library(igraph)

#for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15;do Rscript ConcorTest_clean1.r Tree="/home/esl/Desktop/scRNAphyloData/Final/AML1per/LE${i}.all-cell-types.iqtree.tree" CloID="/home/esl/Desktop/scRNAphyloData/Final/AML1per/Initial/LE${i}/LE${i}/initial/clusters.tsv.gz"; done
##########
test_pair <- function(clone_A,clone_B,states,D,N_PERM) {
    keep <- states %in% c(clone_A, clone_B)
    labs <- states[keep]
    tips <- names(labs)
    nA <- sum(labs == clone_A)
    nB <- sum(labs == clone_B)

    A_tips <- tips[labs == clone_A]
    B_tips <- tips[labs == clone_B]
    observed <- mean(D[A_tips, B_tips])

    null <- numeric(N_PERM)
    for (i in seq_len(N_PERM)) {
        perm <- sample(labs)
        A_perm <- tips[perm == clone_A]
        B_perm <- tips[perm == clone_B]
        null[i] <- mean(D[A_perm, B_perm])
    }
    p_value <- (1 + sum(null >= observed)) /
               (N_PERM + 1)

    data.frame(
        Clone1 = clone_A,
        Clone2 = clone_B,
        N1 = nA,
        N2 = nB,
        ObservedDistance = observed,
        NullMean = mean(null),
        NullSD = sd(null),
        P = p_value,
        stringsAsFactors = FALSE
    )
}

test_cell<-function(Clade=Clade, D=D){
    CloneLs=unique(Clade$clone_id)
    Head=c('cell_id','clone_id',paste0(rep(CloneLs, each = 2), c("_mean", "_p")),'filter')
    Start='y'
    for (i in rownames(Clade)){
        Cell=Clade[i,1]
        TarClone=Clade[i,2]
        TarCloTips=Clade[Clade$clone_id==TarClone & Clade$cell_id!=Cell,]$cell_id
        TarMeanDis=mean(D[Cell,TarCloTips])
        Res=c(Cell,TarClone)        
        Fil='n'
        for (Clone in CloneLs){
            if (Clone==TarClone){Res=c(Res,TarMeanDis,NA)}else{
                CloTips=Clade[Clade$clone_id==Clone,]$cell_id
                MeanDis=mean(D[Cell,CloTips])
                if (MeanDis<=TarMeanDis){
                    Res=c(Res,MeanDis,1)
                    Fil='y'
                } else {
                    ResT=t.test(D[Cell,TarCloTips],D[Cell,CloTips],alternative = "less")
                    P=ResT$p.value
                    Res=c(Res,MeanDis,P)
                    if (P>0.05){Fil='y'}
                }                
            }
        }    
        Res=c(Res,Fil)    
        names(Res)=Head
     #   Add=data.frame(Res)
        if (Start=='y'){
                All=Res
                Start='n'
        }else{All=rbind(All,Res)}
    }
    return(
      list(
        table = All
      )
    )  
}

test_within<-function(Clade=FillInfo, D=D){
    CloneLs=unique(Clade$clone_id)
    
    Start='y'
    Head0=c('clone_id','mean','sd','N_afterFilter')
    for (i in CloneLs){
        CellLs=Clade[Clade$clone_id==i & Clade$filter=='n',]$cell_id
        DisLs <- D[CellLs, CellLs]
        DisLs[lower.tri(DisLs, diag = TRUE)] <- NA
        DisLs <- DisLs[!is.na(DisLs)]
        MeanDis=mean(DisLs)
        SDDis=sd(DisLs)
        Res=c(i,MeanDis,SDDis,length(CellLs))     
        names(Res)=Head0
        if (Start =='y'){
            All=Res
            Start='n'
        }else{All=rbind(All,Res)}
    }    
    Start='y'
    All=data.frame(All)
    CloneLs=All[All$N_afterFilter>100,]$clone_id
    Head=c('clone_id','mean','sd','N_afterFilter',paste0(rep(CloneLs, each = 1), c("_p")),'spread')
    for (TarClone in CloneLs){
        In=c(TarClone,All[All$clone_id==TarClone,]$mean,All[All$clone_id==TarClone,]$sd,All[All$clone_id==TarClone,]$N_afterFilter)
        Spread='y'
        for (AnoClone in CloneLs){
            if (TarClone==AnoClone){ In=c(In,NA)}else{
                if (All[All$clone_id==TarClone,]$mean<=All[All$clone_id==AnoClone,]$mean){
                   In=c(In,1)
                   Spread='n'
                   }else{
                       TtestP=Ttest(mean1=All[All$clone_id==TarClone,]$mean,sd1=All[All$clone_id==TarClone,]$sd,n1=All[All$clone_id==TarClone,]$N_afterFilter,mean2=All[All$clone_id==AnoClone,]$mean,sd2=All[All$clone_id==AnoClone,]$sd,n2=All[All$clone_id==AnoClone,]$N_afterFilter)
                       In=c(In,TtestP)
                       print (TtestP)
                       if (TtestP>=0.05){Spread='n'}
                }
            }
        }
        In=c(In,Spread)
        names(In)=Head
        if (Start =='y'){
            All1=In
            Start='n'
        }else{All1=rbind(All1,In)}    
    }    
    return(
      list(
        table = All1
      )
    )  
}

Ttest<-function(mean1,sd1,n1,mean2,sd2,n2){
    mean1 <- as.numeric(mean1)
    sd1   <- as.numeric(sd1)
    n1    <- as.numeric(n1)
    mean2 <- as.numeric(mean2)
    sd2   <- as.numeric(sd2)
    n2    <- as.numeric(n2)
    se <- sqrt(sd1^2 / n1 + sd2^2 / n2)
    t_stat <- (mean1 - mean2) / se
    df <- (sd1^2 / n1 + sd2^2 / n2)^2 /
      ((sd1^2 / n1)^2 / (n1 - 1) +
       (sd2^2 / n2)^2 / (n2 - 1))
    p_value <- 2 * pt(-abs(t_stat), df)
    return (p_value)
}
###main
###Data
# Create a 10-tip tree
#set.seed(123)

#TreeR2 <- rtree(10)

# Give the tips meaningful IDs
#TreeR2$tip.label <- paste0("Tip_", 1:10)

# Define two distinct states
#tip_states2 <- c(
#  Tip_1 = 0,
#  Tip_2 = 0,
#  Tip_3 = 0,
#  Tip_4 = 1,
#  Tip_5 = 1,
#  Tip_6 = 1,
#  Tip_7 = 0,
#  Tip_8 = 1,
#  Tip_9 = 0,
#  Tip_10 = 1
#)


#Dir='/home/esl/Desktop/scRNAphyloData/SPECTRUM/OV003/'
#Tree=paste0(Dir,'MERGED.Ovarian-cancer.fasttree.tree')
#CloID=paste0(Dir,'MERGED.Ovarian-cancer.clones.manifold.tsv.gz')
#####
#Tree='/home/esl/Desktop/scRNAphyloData/Final/AML1per/LE9.all-cell-types.iqtree.tree'
#CloID='/home/esl/Desktop/scRNAphyloData/Final/AML1per/Initial/LE9/LE9/initial/clusters.tsv.gz'
args <- commandArgs(trailingOnly = TRUE)

Tree <- sub("^Tree=", "", args[grep("^Tree=", args)])
CloID <- sub("^CloID=", "", args[grep("^CloID=", args)])

cat("Tree:", Tree, "\n")
cat("CloID:", CloID, "\n")
####
Clade=read.delim(CloID,sep='\t')
#Clade$clone_id[Clade$clone_id %in% c(1, 4)] <- 14 #need to manually merge clones
tip_states=Clade$clone_id
names(tip_states)=Clade$cell_id
TreeR=read.tree(Tree)
TreeR <- root(TreeR, outgroup = "REFERENCE", resolve.root = TRUE)
plot(TreeR, show.tip.label = FALSE)

length(TreeR$tip.label)
length(tip_states)

setdiff(TreeR$tip.label, names(tip_states))   # tips missing states
setdiff(names(tip_states), TreeR$tip.label)   # extra states not in tree
tips_to_drop <- setdiff(TreeR$tip.label, names(tip_states))
TreeR2 <- drop.tip(TreeR, tips_to_drop)

tip_states2 <- tip_states[TreeR2$tip.label]
#OriginalTree=TreeR2


#start
Clade <- data.frame(
  cell_id = TreeR2$tip.label,
  clone_id = as.numeric(tip_states2[TreeR2$tip.label]),
  stringsAsFactors = FALSE
)

#Clade
p=ggtree(TreeR2)
p$data=merge(p$data,Clade ,by.x='label',by.y='cell_id',all.x=TRUE)
head(p$data)
A=p+geom_tippoint(aes(color = as.factor(clone_id)), size = 1) 
ggsave(
  paste0(Tree, ".png"),
  plot = A,
  width = 10,
  height = 10,
  dpi = 300
)


###Concordance test and merge clones
tr <- TreeR2

## Match clone labels to tree tips
states <- tip_states2[tr$tip.label]

stopifnot(all(!is.na(states)))

clones <- sort(unique(states))

cat("Number of tips:", Ntip(tr), "\n")
cat("Number of clones:", length(clones), "\n")
print(table(states))

tr_unweighted <- tr
tr_unweighted$edge.length <- rep(1, nrow(tr_unweighted$edge))

## Number of edges between every pair of tips
D <- cophenetic(tr_unweighted)

pairs <- combn(clones, 2, simplify = FALSE)

results <- do.call(
    rbind,
    lapply(
        pairs,
        function(x) {
            test_pair(
                clone_A = x[1],
                clone_B = x[2],
                states = states,
                D = D,
                N_PERM = 100
            )
        }
    )
)


## Multiple-testing correction

results$FDR <- p.adjust(results$P, method = "BH")

results$Difference=results$ObservedDistance- results$NullMean 
results$Classification <- ifelse(
  results$Difference < 1, "Overlapping",
  ifelse(results$Difference <= 2, "Weak overlap", "Separated")
)
write.table(results,paste0(Tree,'Concordance.txt'),sep='\t',quote=FALSE,row.names = FALSE)
##merge and cleaning
OverPair=results[results$Classification== "Overlapping",]
pairs <- OverPair[OverPair$Classification == "Overlapping",
                  c("Clone1", "Clone2")]

# Create graph
g <- graph_from_data_frame(
  pairs,
  directed = FALSE
)

# Find connected components
comp <- components(g)

# Create group table
groups <- data.frame(
  Clone = names(comp$membership),
  Group = paste0("Group_", comp$membership)
)

Gls=unique(groups$Group)
for (G in Gls){
   Cls=paste0(groups[groups$Group == G, ]$Clone, collapse = "")
   groups[groups$Group == G, "Group"] <- Cls

}
groups <- groups[order(groups$Group, groups$Clone), ]

groups
LineLs=rownames(groups)
for (Line in LineLs){
 C=groups[Line,]$Clone
 G=groups[Line,]$Group
 Clade$clone_id[Clade$clone_id %in% c(C)] <- G #need to manually merge clones
}

p=ggtree(TreeR2)
p$data=merge(p$data,Clade ,by.x='label',by.y='cell_id',all.x=TRUE)
head(p$data)
A=p+geom_tippoint(aes(color = as.factor(clone_id)), size = 1) 
ggsave(
  paste0(Tree, "Merge.png"),
  plot = A,
  width = 10,
  height = 10,
  dpi = 300
)

###Find root clone (before filtering). do not need?
#depth <- node.depth.edgelength(tr_unweighted)
#top30 <- tr_unweighted$tip.label[
#    order(depth[seq_len(Ntip(tr_unweighted))])[1:30]
#]
#RootCells <- Clade[Clade$cell_id %in% top30, ]
#MajorClo=names(which.max(table(RootCells$clone_id)))
#MajorPro=max(table(RootCells$clone_id))/sum(table(RootCells$clone_id))
#if (MajorPro<=0.5){MajorClo='None'}
#print (MajorClo)
#print (MajorPro)

###Clean tree using cell position
resultsC <- test_cell(Clade=Clade, D=D)
FillInfo=data.frame(resultsC$table)
write.table(FillInfo,paste0(Tree,'Concordance.txt'),sep='\t',quote=FALSE,row.names = FALSE)

tips_to_drop=FillInfo[FillInfo$filter=='y',]$cell_id
TreeR3 <- drop.tip(TreeR2, tips_to_drop)
p=ggtree(TreeR3)
p$data=merge(p$data,Clade ,by.x='label',by.y='cell_id',all.x=TRUE)
head(p$data)
A=p+geom_tippoint(aes(color = as.factor(clone_id)), size = 1) 
ggsave(
  paste0(Tree, "Filter.png"),
  plot = A,
  width = 10,
  height = 10,
  dpi = 300
)
Out=data.frame(cell_id=p$data[p$data$isTip==TRUE,]$label,clone_id=p$data[p$data$isTip==TRUE,]$clone_id)
write.table(Out,paste0(Tree,'ConcordanceCellAnno.txt'),sep='\t',quote=FALSE,row.names = FALSE)
write.tree(TreeR3, file = paste0(Tree,'Concordance.nwk'))

##make summary
Original=data.frame(tip_states)
Original$cell_id=row.names(Original)
All=merge(Original,Out,by='cell_id',all.x=TRUE)
Count <- aggregate(cell_id ~ tip_states + clone_id,
                   data = All,
                   FUN = function(x) length(unique(x)))
IncNA <- aggregate(cell_id ~ tip_states,
                data = All,
                FUN = function(x) length(unique(x)))
names(IncNA)[2] <- "Include_removed"
Count <- merge(Count, IncNA, by = "tip_states", all.y=TRUE)
names(Count)=c('Original_clone_id','Merged_clone_id','cell_count_FinalTree','cell_count_Original')
write.table(Count,paste0(Tree,'CellCount.txt'),sep='\t',quote=FALSE,row.names = FALSE)

##Do not use below
###check some clone widely distributed than others
#resultsWithin=test_within(Clade=FillInfo, D=D)
#DistTa=data.frame(resultsWithin$table)
#write.table(DistTa,paste0(Tree,'Spread.txt'),sep='\t',quote=FALSE,row.names = FALSE)
#SpreadClone=DistTa[DistTa$spread=='y',]$clone_id
#SpreadCloneTips=FillInfo[FillInfo$clone_id==SpreadClone,]$cell_id
#if (length(SpreadClone)==1){ #there is a spread clone
#    tips_to_drop=FillInfo[FillInfo$filter=='y' & FillInfo$clone_id!=SpreadClone,]$cell_id
#    TreeR3 <- drop.tip(TreeR2, tips_to_drop)
#    p=ggtree(TreeR3)
#    p$data=merge(p$data,Clade ,by.x='label',by.y='cell_id',all.x=TRUE)
#    head(p$data)
#    A=p+geom_tippoint(aes(color = as.factor(clone_id)), size = 1) 
#    ggsave(
#      paste0(Tree, "Filter_nonSpread.png"),
#      plot = A,
#      width = 10,
#      height = 10,
#      dpi = 300
#    )
#    CloneLs=unique(DistTa$clone_id)
#    RmLs=c()
#    for (Clone in CloneLs){
#        if (Clone!=SpreadClone){
#            Tips=FillInfo[FillInfo$filter=='n' & FillInfo$clone_id==Clone,]$cell_id
#            mrca_node <- getMRCA(TreeR3, Tips)
#            CladeTips <- TreeR3$tip.label[
#               Descendants(TreeR3, mrca_node, type = "tips")[[1]]
#            ]
#            SpreadTipIn=intersect(CladeTips,SpreadCloneTips)
#            if (length(SpreadTipIn)>length(Tips)){
#                RmLs=c(RmLs,Tips)
#            }else{
#                RmLs=c(RmLs,SpreadTipIn)
#            }
#        }
#    }
#    RmLs=unique(RmLs)
#    TreeR3 <- drop.tip(TreeR3, RmLs)
#    p=ggtree(TreeR3)
#    p$data=merge(p$data,Clade ,by.x='label',by.y='cell_id',all.x=TRUE)
#    head(p$data)
#    A=p+geom_tippoint(aes(color = as.factor(clone_id)), size = 1) 
#    ggsave(
#      paste0(Tree, "Filter_Final.png"),
#      plot = A,
#      width = 10,
#      height = 10,
#      dpi = 300
#    )
#    FillInfo$filter[FillInfo$cell_id %in% RmLs] <- "y"
#    write.table(FillInfo,paste0(Tree,'Concordance_Final.txt'),sep='\t',quote=FALSE,row.names = FALSE)
#    Out=data.frame(cell_id=p$data[p$data$isTip==TRUE,]$label,clone_id=p$data[p$data$isTip==TRUE,]$clone_id)
#    write.table(Out,paste0(Tree,'ConcordanceCellAnno_Final.txt'),sep='\t',quote=FALSE,row.names = FALSE)
#    write.tree(TreeR3, file = paste0(Tree,'Concordance_Final.nwk'))

#}




