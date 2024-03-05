rm(list=ls())

################################################################################
# Recurrent Genomic Selection: Wheat data set
# Strategy Random
################################################################################

##        RND   CB80-L MAX-L 
## --------------------------
## SE-IL  GEV   <-     <-
## CR-IL  RND   CB80   MAX
## CR-F1  RND   <-     <- 
## SE-C1  GEV   <-     <- 
## CR-C1  RND   <-     <- 
## SE-C2  GEV   <-     <- 
## CR-C2  RND   <-     <-    
## SE-C3  GEV   <-     <- 
rm(list = ls())
library ("SelectionTools")
library ("ggplot2");library ("gridExtra")
#st.script.dir <- "~/1-daten/p111-selection-tools/SelectionTools-dev/src/"
#source(paste(st.script.dir,"SelectionTools.R",sep=""))

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
# Random intermating of the selected P 

population.copy("tmp","Psel")
for (ii in 1:144) {
  nme <- sprintf("p%03i",ii)
  population.divide (nme,"tmp",1)
}

idx <- sample(1:144,144)

for (ii in 1:72) {                                # CR-L: Random
  nme1 <- sprintf("p%03i",idx[ii])
  nme2 <- sprintf("p%03i",idx[72+ii])
  nme3 <- sprintf("f1%02i",ii)
  cross(nme3,nme1,nme2,1)
}

for (ii in 1:36) {                                # CR-1: Random
  nme1 <- sprintf("f1%02i",ii)
  nme2 <- sprintf("f1%02i",36+ii)
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

idx <- sample(1:72,72)                     
                                           
for (ii in 1:36) {
  nme1 <- sprintf("s%03i",idx[ii])
  nme2 <- sprintf("s%03i",idx[36+ii])
  nme3 <- sprintf("SYNn%02i",ii)
  cross(nme3,nme1,nme2,10)                 # CR-2 / CR-3: Random
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
# Save simulation results
########################################################################

# write.table (e,"data/c002-110-e001.dta" ) #  200 REPS
# write.table (e,"data/c002-110-e002.dta" ) #  200 REPS
write.table (e,"data/c002-110-e010.dta" ) # 1000 REPS

########################################################################