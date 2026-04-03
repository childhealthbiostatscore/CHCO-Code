#### explore msm package for TODAY study ####
library(msm)
options(scipen = 999)
home_dir = paste0("/Users/",Sys.info()["user"],"/Library/CloudStorage/OneDrive-UW/Laura Pyle's files - Biostatistics Core Shared Drive/TODAY subaward")
setwd(home_dir)

# dataset with CKM stges, but 1 row pp in a summary fashion:
dat.ckm<-read.csv("./CKM/Results/CKM stage summary.csv")

# long dataset with CKM stages:
dat<-read.csv("./Data clean/comorb_ckm2_for_hce.csv")
nrow(dat)
length(unique(dat$RELEASEID))

#CKM stage variable:
#CKM_syn??

table(dat$CKM_syn)

ckm_sub<-dat[,c("RELEASEID","days","CKM_syn","CKM_syn2","CKM_syn_base","CKM_syn_base2",
                "progress_CKM","progress_CKM2")]

# it seems like CKM_syn is the main variable?? ask laura
#is stage 4 = death? how is death incorporated in this long dataset?

dat$ckm_num<-NA
dat$ckm_num[dat$CKM_syn=="Stage 2"]<-1
dat$ckm_num[dat$CKM_syn=="Stage 2+"]<-2
dat$ckm_num[dat$CKM_syn=="Stage 3"]<-3
dat$ckm_num[dat$CKM_syn=="Stage 4"]<-4

### create a statetable:
statetable.msm(CKM_syn, RELEASEID, data=dat)

#patients with only one observation?
ckm_1<- data.frame(dat %>%
  group_by(RELEASEID) %>%
  summarise(n_rows = n()
  ))
table(dat$n_rows)
dat<-merge(dat,ckm_1,by="RELEASEID",all.x=T)

#remove those with only one observation for now:
dat<-subset(dat,dat$n_rows>1)

#initial values? this doesn't converge
# twoway4.q <- rbind(c(0, 0.25, 0, 0.25), c(0.166, 0, 0.166, 0.166), + c(0, 0.25, 0, 0.25), c(0, 0, 0, 0))
# rownames(twoway4.q) <- colnames(twoway4.q) <- c("Stage 2", "Stage 2+", "Stage 3", "Stage 4")
Q<-crudeinits.msm(ckm_num ~ days, subject = RELEASEID, data = dat, qmatrix = twoway4.q)
test_msm<-msm(ckm_num ~ days, subject = RELEASEID, data = dat, qmatrix = Q)
#probability of transitioning in the next year:
pmatrix.msm(test_msm, t = 365, ci = "normal")

#add a covariate as an example:
dat<-merge(dat[,-which(colnames(dat)=="sex")],dat.ckm[,c("RELEASEID","sex")],by="RELEASEID",all.x=T)
sex_msm<-msm(ckm_num ~ days, subject = RELEASEID, data = dat, qmatrix = Q,
             covariates = ~ sex)

hazard.msm(sex_msm)
