################################################################################
# Recurrent Genomic Selection: Wheat data set
# Strategy EMBV
################################################################################

##              222    223    224    227    226
##        RND   MAX-F  MSD-1  MSD-2  UC-OCS EMBV
## ------------------------------------------------- 
## SE-IL  GEV   <-     <-     <-     <-     <-
## CR-IL  RND   MAX    MAX    MAX    UC-OCS <-
## CR-F1  RND   MAX    MAX    MAX    UC-OCS <-  
## SE-C1  GEV   <-      |      |     |      <-
## CR-C1  RND   <-     MSD    MSD    UC-OCS EMBV 
## SE-C2  GEV   <-     <-      |     |      <-
## CR-C2  RND   <-     <-     MSD    UC-OCS EMBV
## SE-C3  GEV   <-     <-     <-     <-     <-

calc.embv <- function(pop,N,G,REP){
    l <- list.populations()
    n <- l[ pop == l[,1],2]
    embv <- rep(0,n)
    gv   <- rep(0,REP)
    population.copy("EMBVsplit",pop)
    for (ii in 1:n) {
        population.divide ("EMBVSC","EMBVsplit",1)
        for (r in 1:REP) { 
            dh ("EMBVDH","EMBVSC",G)
            genotype.population("EMBVDH")
            evaluate.population("EMBVDH", "yld")
            population.sort("EMBVDH", decreasing=TRUE) 
            population.divide("EMBVDHsel", "EMBVDH", N)
            gval <-get.population.gvalue("EMBVDHsel")
            gv[r] <- mean(gval$gvalue)
        }
        embv[ii] <- mean(gv)
    }
    return(embv)
}

library ("SelectionTools")
#library("sqldf")

################################################################################
st.input.dir  <- "input3"
st.output.dir <- "output3"
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

gs.esteff.rr ( method="BLUP", data.set="PP")  
yld.eff <- gs.return.effects (data.set="PP" ) 

gs.write.pseff( yld.eff, file="data/effect.file1.eff")    

#gs.plot.validation("PP","PP")

st.set.simpop ( pop.name="PP", data.set="PP" )
load.effmap("yld","data/effect.file1.eff")

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
# Loop for the simulation
###########################################
NREP <-200
e    <- NULL

for (REP in 1:NREP) {

cat (sprintf("%05i\r",REP))
# Random intermating of the selected P 

population.copy("tmp","Psel")
for (ii in 1:144) {
  nme <- sprintf("p%03i",ii)
  population.divide (nme,"tmp",1)
}

idx <- sample(1:144,144)

for (ii in 1:72) {
  nme1 <- sprintf("p%03i",idx[ii])
  nme2 <- sprintf("p%03i",idx[72+ii])
  nme3 <- sprintf("f1%02i",ii)
  cross(nme3,nme1,nme2,1)
}


for (ii in 1:36) {
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

#genotype.population("SYN")                 # Instead of selection accoriding
#evaluate.population("SYN", "yld")          # to the GEBV

embv <- calc.embv (pop="SYN", N=2, G=6, REP=10)  # We use the EMBV
set.population.gvalue("SYN",embv) 

population.sort("SYN", decreasing=TRUE)    # Selection accoriding to
population.divide("SYNsel", "SYN", 72)     # GEGV

nme <- sprintf("SYN%1isel",C-1)
copy.population(nme,"SYNsel")

population.copy("split","SYNsel")
for (ii in 1:72) {
  nme <- sprintf("s%03i",ii)
  population.divide (nme,"split",1)
}

idx <- sample(1:72,72)                     # Random mating of the selected
                                           # fraction
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

#genotype.population("SYN")                 
#evaluate.population("SYN", "yld")

embv <- calc.embv (pop="SYN", N=2, G=6, REP=20)  
set.population.gvalue("SYN",embv) 

population.sort("SYN", decreasing=TRUE)    
population.divide("SYNsel", "SYN", 36)     
copy.population("SYN3sel","SYNsel")

# 6 DH Lines per selected SYN3 plant


remove.population("DH")
population.copy("split","SYNsel")
for (ii in 1:36) {
    population.divide ("ss","split",1)
    dh ("dd","ss",6)
    append.population("DH","dd")
}

eval <- c("PP","Psel","SYN1","SYN1sel","SYN2","SYN2sel","SYN3","SYN3sel","DH")

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


#m <- tapply(d$y,d$gen,mean);
#e <- merge(data.frame(gen=names(m),y=m),v)
#rownames(e) <- c()
#e$Strategy <- strategy
#e$Crop <- crop
#e$Date <- date()

# Save simulation results in data base
#conn <- dbConnect(RSQLite::SQLite(), dbfile)
#dbWriteTable(conn, rfile, e, append=TRUE)
#dbDisconnect(conn)

} # for (REP in 1:NREP)
write.table(e,"data/228best.dta")
