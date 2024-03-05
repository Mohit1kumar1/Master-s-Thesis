################################################################################
# Recurrent Genomic Selection: wheat data set
# Strategy GDI-L
################################################################################

##        110   111     113   112   
                             
##        RND   CB80-L GDI-L  MAX-L 
## ---------------------------------
## SE-IL  GEV   <-     <-     <-    
## CR-IL  RND   CB80   GDI    MAX   
## CR-F1  RND   <-     <-     <-    
## SE-C1  GEV   <-     <-     <-    
## CR-C1  RND   <-     <-     <-    
## SE-C2  GEV   <-     <-     <-    
## CR-C2  RND   <-     <-     <-    
## SE-C3  GEV   <-     <-     <-           

sel.crs <- function (d.sorted,ncp,nct,old.crs=NULL)
                                        # ncp Number of crosses per parent  
                                        # nct Total number of crosses       
    {
        indx <- rep(FALSE,nrow(d.sorted))           

        if (!is.null(old.crs))
          for (i in 1:nrow(old.crs))
            for (j in 1:nrow(d.sorted))
              if ((old.crs[i,1]==d.sorted[j,1])&&(old.crs[i,2]==d.sorted[j,2]))
                indx[j] <- TRUE
        
        new.crs <- 0
        for (i in 1:nrow(d.sorted))  {
            crosses.so.far <- c(d.sorted$OTU1[indx],d.sorted$OTU2[indx])

            search1 <- paste("\\b",d.sorted$OTU1[i],"\\b",sep="")
            search2 <- paste("\\b",d.sorted$OTU2[i],"\\b",sep="")
            
            if ( (new.crs<nct)                               &&
                (ncp > length(grep(search1,crosses.so.far))) &&
                (ncp > length(grep(search2,crosses.so.far))) &&
                (FALSE == indx[i])                                 )
                {
                    indx[i] <- TRUE
                    new.crs = 1 + new.crs
                }
        }
        crosses <- d.sorted[indx,]

        return (crosses)
    }


library ("ggplot2"); library ("gridExtra")

library ("SelectionTools")
#st.script.dir <- "~/1-daten/p111-selection-tools/SelectionTools-dev/src/"
#source(paste(st.script.dir,"SelectionTools.R",sep=""))

st.input.dir  <- "input"
st.output.dir <- "output1"
st.set.info.level (-2)
gs.set.num.threads(4)

st.read.marker.data ("marker.data.mpo",format="m", data.set ="PP")
st.read.map         ("map.file.txt",m.stretch=1, 
                                    skip=1, format="mcp", data.set = "PP")
st.read.performance.data ("yld.txt", data.set="PP")

st.restrict.marker.data   ( NoAll.MAX = 2 , data.set ="PP"  )  
st.restrict.marker.data   ( MaMis.MAX = 0.1, data.set ="PP" )  
st.restrict.marker.data   ( InMis.MAX = 0.1, data.set ="PP" )  
st.restrict.marker.data   ( ExHet.MIN = 0.1, data.set ="PP" )  


st.copy.marker.data("PA", "PP" )

gs.esteff.rr("BLUP",data.set = "PA")
yld.eff <- gs.return.effects(data.set="PA")

st.set.simpop ( pop.name="PP", data.set="PP" ) 
load.effmap("yld","data/effect.file.eff")
################################################################################

genotype.population("PP")
evaluate.population("PP", "yld")

# Select the best 144 P lines / GEGV

copy.population("PA","PP")
genotype.population("PA")
evaluate.population("PA", "yld")
population.sort("PA", decreasing=TRUE) 
population.divide("Psel", "PA", 144)                # SE-L: GEGV

###########################################
NREP <- 1000
###########################################

e    <- NULL

for (REP in 1:NREP) {
    
cat (sprintf("%05i/%05i\r",REP,NREP))

# Cross selected partental lines                      # CR-IL: GDI

st.get.simpop("Psel","Psel")
gs.set.effects(eff=yld.eff,data.set="Psel")
#gs.cross.eval.ma(data.set="Psel")
gs.cross.eval.gd(data.set="Psel")
crs.Psel  <- gs.cross.info(data.set="Psel",sortby ="gd")
d.sorted <- data.frame(OTU1     =crs.Psel$P1No,
                       OTU2     =crs.Psel$P2No,
                       Measure  =crs.Psel$gd,stringsAsFactors=FALSE)
crs     <- sel.crs (d.sorted,ncp=1,nct=72)  

population.copy("tmp","Psel")
for (ii in 1:144) {
  nme <- sprintf("p%03i",ii)
  population.divide (nme,"tmp",1)
}

for (ii in 1:72) {
  nme1 <- sprintf("p%03i",crs[ii,1])
  nme2 <- sprintf("p%03i",crs[ii,2])
  nme3 <- sprintf("f1%02i",ii)
  cross(nme3,nme1,nme2,1)
}

idx <- sample(1:72,72)

for (ii in 1:36) {                                  # CR-1: Random
  nme1 <- sprintf("f1%02i",idx[ii]) 
  nme2 <- sprintf("f1%02i",idx[36+ii])
  nme3 <- sprintf("SYN1%02i",ii)
  cross(nme3,nme1,nme2,10)
}

remove.population("SYN1")
for (ii in 1:36) {
  nme3 <- sprintf("SYN1%02i",ii)
  append.population("SYN1",nme3)
}

# Select the best 72 SYN1 plants / GEGV

copy.population("SYN","SYN1")

for (C in 2:3 ) {

genotype.population("SYN")                 # Starts with a population SYN
evaluate.population("SYN", "yld")
population.sort("SYN", decreasing=TRUE)    
population.divide("SYNsel", "SYN", 72)     # SE-1 / SE-2: GEGV

nme <- sprintf("SYN%1isel",C-1)
copy.population(nme,"SYNsel")

population.copy("split","SYNsel")
for (ii in 1:72) {
  nme <- sprintf("s%03i",ii)
  population.divide (nme,"split",1)
}

idx <- sample(1:72,72)                     # CR-2 / CR-3: Random
                                           
for (ii in 1:36) {
  nme1 <- sprintf("s%03i",idx[ii])
  nme2 <- sprintf("s%03i",idx[36+ii])
  nme3 <- sprintf("SYNn%02i",ii)
  cross(nme3,nme1,nme2,10)                 
}

remove.population("SYN")
for (ii in 1:36) {
  nme3 <- sprintf("SYNn%02i",ii)
  append.population("SYN",nme3)
}

nme <- sprintf("SYN%1i",C)
copy.population (nme,"SYN")                # Generates the new SYN population

} # for (C in 2:3 )

# Select 36 SYN3 plants

genotype.population("SYN")                 
evaluate.population("SYN", "yld")
population.sort("SYN", decreasing=TRUE)    
population.divide("SYNsel", "SYN", 36)     
copy.population("SYN3sel","SYNsel")       # SE-3

# 6 DH Lines per selected SYN3 plant

remove.population("DH")
population.copy("split","SYNsel")
for (ii in 1:36) {
    population.divide ("ss","split",1)
    dh ("dd","ss",6)
    append.population("DH","dd")
}

eval <- c("PP","SYN1","SYN2","SYN3","DH")

v <- d <- NULL
for (i in 1:length(eval))
{
    genotype.population ( eval[i] )                 
    evaluate.population ( eval[i], "yld")
    P <- get.population.gvalue(eval[i])$gvalue
    d <- rbind(d,data.frame (gen=eval[i],y=P))
#
    D <- eval[i]
    st.get.simpop(D,D)
    st.def.hblocks ( hap=5, hap.unit=2, data.set=D ) 
    st.recode.hil  (data.set=D)
    div <- mean(st.marker.data.statistics(data.set=D)$marker.list$ExHet)
    v <- rbind(v,data.frame (gen=eval[i],d=div))
}
m <- tapply(d$y,d$gen,mean);
e <- rbind(e,merge(data.frame(gen=names(m),y=m),v))
rownames(e) <- c()

} # for (REP in 1:NREP)

e
 ########################################################################
#write.table (e,"data/c002-113-e001.dta" ) #  200 REPS
write.table (e,"data/c002-113-e010.dta" ) #  1000 REPS