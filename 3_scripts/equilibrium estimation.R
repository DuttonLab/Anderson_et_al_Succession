#libraries
library(here)
library(matlib)
library(tidyverse)

theme.mine=theme(plot.title = element_text(face="bold", size=24),
                 text = element_text(size=rel(5)),
                 axis.text.x = element_text(size=rel(5)),
                 axis.text.y = element_text(size=rel(5)),
                 legend.text=element_text(size=rel(5)),
                 panel.background = element_blank(),
                 axis.line=element_line(linetype=1, color="black"))

## matrix approach to solving the systems:: https://www.jstor.org/stable/1935355?seq=1#metadata_info_tab_contents
specvec = c("BC9", "BC10", "JB5", "JB7", "X135E", "JBC", "JB370")

## Ph 5 ----
raw = read.csv(here("4_res", "coefficient-estimates-pH5.csv"))

# make vector of Ks, ph5
ktemp = raw[raw$name=="K",]
kvec5 = ktemp$est
names(kvec5)=ktemp$spec

#make interaction matrix
A = diag(7)
colnames(A) = specvec
rownames(A) = specvec

# We know that K = Nstar %*% A, so 
# Nstar = K %*% inv(A)

for(cur.focal in specvec){
  dat.cur = raw %>% 
    filter(spec == cur.focal) %>% 
    filter(name != "K")
  A[dat.cur$name, cur.focal] = dat.cur$est
}
#saving a version for later, with pH identity
A.5=A
#remember, the part that matters is the interactions (as NEGATIVE values) (interactions)
# aug = cbind((kvec), -A.5)
# aug = rbind(c(1, rep(0,ncol(aug)-1)),
#             aug)
# eigen(aug)$values
# eigen(cbind((kvec), -A.5))$values
temp=A.5
diag(temp)=-kvec5
kmat = matrix((kvec5), ncol=1)
kmat=kmat[,rep(1,nrow(kmat))]
temp = -temp/kmat
eigen(-A.5)$values
eigen(temp)$values

# det(A)
equilib.5 = as.vector(kvec5 %*% inv(A))
names(equilib.5)=specvec


## sanity test: how does that compare to our estimates from the prediction approach
pred.5 = read.csv(here("4_res","predictions-from-pH5.csv"))

equilib.predtemp = pred.5 %>% 
  group_by(spec) %>% 
  summarize(dens = mean(fit))
equilib.pred5 = equilib.predtemp$dens
names(equilib.pred5)=equilib.predtemp$spec
equilib.pred5=equilib.pred5[names(equilib.5)]
rbind(equilib.5, equilib.pred5)

## How is stability in system missing the zero-density species?
A.mod = A.5[names(equilib.5)[equilib.5>0], names(equilib.5)[equilib.5>0]]
eigen(-A.mod)$values

## Repeat with pH 7 data ----

raw = read.csv(here("4_res", "coefficient-estimates-pH7.csv"))

# make vector of Ks
ktemp = raw[raw$name=="K",]
kvec7 = ktemp$est 
names(kvec7)=ktemp$spec

#make interaction matrix
A = diag(7)
colnames(A) = specvec
rownames(A) = specvec

# We know that K = Nstar %*% A, so 
# Nstar = K %*% inv(A)

for(cur.focal in specvec){
  dat.cur = raw %>% 
    filter(spec == cur.focal) %>% 
    filter(name != "K")
  A[dat.cur$name, cur.focal] = dat.cur$est
}
#saving a version for later, with pH info
A.7 = A
# we have a problem: NA for effect of X135E on JBC. BUT! We also know there is 
# no X135E in the final community observations, so as a quick fix let's cut it out. 

Aminus = A[-which(rownames(A)=="X135E"), -which(colnames(A)=="X135E")]
kminus = kvec7[-which(names(kvec7)=="X135E")]

#stability:
eigen(-Aminus)$values


##STABLE!!

# det(A)
equilib.7 = as.vector(kminus %*% inv(Aminus))
names(equilib.7)=names(kminus)
# stability:
A.mod = Aminus[names(equilib.7)[equilib.7>0], names(equilib.7)[equilib.7>0]]
eigen(-A.mod)$values
#stable!!

## sanity test: how does that compare to our estimates from the prediction approach
pred.7 = read.csv(here("4_res","predictions-from-pH7.csv"))

equilib.predtemp = pred.7 %>% 
  na.omit() %>% 
  group_by(spec) %>% 
  summarize(dens = mean(fit))
equilib.pred7 = equilib.predtemp$dens
names(equilib.pred7)=equilib.predtemp$spec
equilib.pred7=equilib.pred7[names(equilib.7)]
rbind(equilib.7, equilib.pred7)


obs.raw = read.csv(here("1_raw_data","2018_08_21_Bayley_allPW.csv"))
obs.com = obs.raw %>% 
  filter(cond == "community") %>% 
  filter(day == 21) %>% 
  group_by(spec) %>% 
  na.omit() %>% 
  summarise(cfus = median((cfus)))

obs.com$cond="obs"

#put our equilibrium calculations to data frames
equilib.5.df=data.frame(cfus = equilib.5, spec = names(equilib.5))
equilib.7.df=data.frame(cfus = equilib.7, spec = names(equilib.7))
# set negatives to 0s
equilib.5.df$cfus=pmax(equilib.5.df$cfus, 1)
equilib.5.df$cond = "LV pH5"
equilib.7.df$cfus=pmax(equilib.7.df$cfus, 1)
equilib.7.df$cond = "LV pH7"
res = rbind(obs.com, equilib.5.df, equilib.7.df)


## Now some reduced version:
## (a) what if only diutina (X135E) and penicillium (JBC) interactions mattered (pH5)
## ## pH 5 P + D ----
A.pd5=A.5
A.pd5[,!(colnames(A.pd5) %in% c("X135E", "JBC"))]=0
#make sure the diagonals are all exactly 0
A.pd5=A.pd5*(1-diag(7))+diag(7)
A.pd5
equilib.pd5 = as.vector(kvec5 %*% inv(A.pd5))
names(equilib.pd5)=names(kvec5)
cur.df = data.frame(spec = names(equilib.pd5),
                    cfus=equilib.pd5,
                    cond="penicillium + diutina pH5"
                    )
res=rbind(res, cur.df)


## Second version of A: p only
## pH 5 P ----
A.p5=A.5
A.p5[,!(colnames(A.p5) %in% c("JBC"))]=0
#make sure the diagonals are all exactly 0
A.p5=A.p5*(1-diag(7))+diag(7)
A.p5
equilib.p5 = as.vector(kvec5 %*% inv(A.p5))
names(equilib.p5)=names(kvec5)
cur.df = data.frame(spec = names(equilib.pd5),
                    cfus=equilib.pd5,
                    cond="penicillium pH5"
)
res=rbind(res, cur.df)

## pH 5 no Penicillium ----

A.np5=A.5
A.np5[,(colnames(A.np5) %in% c("JBC"))]=0
#make sure the diagonals are all exactly 0
diag(A.np5)=1
A.np5
equilib.np5 = as.vector(kvec5 %*% inv(A.np5))
names(equilib.np5)=names(kvec5)
cur.df = data.frame(spec = names(equilib.np5),
                    cfus=equilib.np5,
                    cond="no penicillium pH5")
res=rbind(res, cur.df)

## pH 5 no Diutina -----
A.nd5=A.5
A.nd5[,(colnames(A.nd5) %in% c("X135E"))]=0
#make sure the diagonals are all exactly 0
diag(A.nd5)=1
A.nd5
equilib.nd5 = as.vector(kvec5 %*% inv(A.nd5))
names(equilib.nd5)=names(kvec5)
cur.df = data.frame(spec = names(equilib.np5),
                    cfus=equilib.np5,
                    cond="no diutina pH5")
res=rbind(res, cur.df)

## pH 7 P ----
A.p7 = A.7
A.p7[,!(colnames(A.p7) %in% c("JBC"))]=0
A.p7=A.p7*(1-diag(7))+diag(7)
Aminus = A.p7[-which(rownames(A.p7)=="X135E"), -which(colnames(A.p7)=="X135E")]
kminus = kvec7[-which(names(kvec7)=="X135E")]

equilib.p7 = as.vector(kminus %*% inv(Aminus))
names(equilib.p7)=names(kminus)
cur.df = data.frame(spec = names(equilib.p7),
                    cfus=equilib.p7,
                    cond="penicillium pH7"
)


## pH 7 no penicillium ----

A.np7 = A.7
A.np7[,(colnames(A.np7) %in% c("JBC"))]=0
diag(A.np7)=1
Aminus = A.p7[-which(rownames(A.p7)=="X135E"), -which(colnames(A.p7)=="X135E")]
kminus = kvec7[-which(names(kvec7)=="X135E")]


equilib.np7 = as.vector(kminus %*% inv(Aminus))
names(equilib.p7)=names(kminus)
cur.df = data.frame(spec = names(equilib.p7),
                    cfus=equilib.p7,
                    cond="no penicillium pH7"
)
res=rbind(res, cur.df)

#quich cleans: make sure naming is consistent, all values are 1 or higher.
res$spec[res$spec=="135E"] = "X135E"
res$cfus = pmax(res$cfus,1)

gp = ggplot(data=res, aes(x=spec, y = cfus, color = cond)) + 
  geom_point(data = res , 
             size=4,
             position = position_dodge2(width=.5))+
  # geom_point(data=equilib.5.df, color="red", position = position_nudge(-.1), size=4)+
  # geom_point(data=equilib.7.df, color="blue", position = position_nudge(.1), size=4)+
  scale_y_log10()+
  ylab("cfus")+
  xlab("")+
  ggtitle("Actual data vs estimated equilibria")+
  theme.mine
gp

ggsave(here("5_figs", "equilib-vs-dat-all.jpg"),
       gp,
       width=12,
       height=8)

write.csv(res,
          file = here("4_res","LV-equilib-calcs.CSV"),
          row.names = FALSE)

