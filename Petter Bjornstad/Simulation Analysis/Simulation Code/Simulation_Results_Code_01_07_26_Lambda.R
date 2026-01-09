#Simulation Results ----
library(reticulate)
use_python("/home/hhampson/miniconda3/bin/python", required = TRUE)

#1. Set Directories ----
# dir.dat <- c("/Users/hhampson/Library/Mobile Documents/com~apple~CloudDocs/BHRM_microbiome/Simulation Code 11_04/Simularion Results/Final")
# dir.code <- c("/Users/hhampson/Library/Mobile Documents/com~apple~CloudDocs/BHRM_microbiome/Simulation Code 11_04/HPC Code")
# dir.lib <- c("/Users/hhampson/Library/Mobile Documents/com~apple~CloudDocs/BHRM_microbiome/Simulation Code 11_04")
# dir.results <- c("/Users/hhampson/Library/CloudStorage/OneDrive-UW/Biostatistics Core Shared Drive/Hailey Hampson/1_Ongoing Projects/Microbiome Manuscript/Simulation Results/Processed Simulation Results")


#2. Load Libraries ----
#Load Libraries
source("Libraries.R")

#3. Set up to load data ----
## Load boto3 and pandas
boto3 <- reticulate::import("boto3")
pd <- reticulate::import("pandas")

## Create an S3 client
# install.packages("jsonlite")  # Install if not already installed
library(jsonlite)  # Load the package

keys <- fromJSON("/home/hhampson/keys.json") # replace with your Lambda username
session <- boto3$session$Session(
  aws_access_key_id = keys$MY_ACCESS_KEY,
  aws_secret_access_key = keys$MY_SECRET_KEY
)

## Create an S3 client with the session
s3 <- session$client("s3", endpoint_url = "https://s3.kopah.uw.edu")

#Bucket
bucket <- "simulation" # bucket name in Kopah

#Load Taxa Z Matrices
#Genus
temp_file <- tempfile(fileext = ".RDS") 
s3$download_file(bucket,"Z.s.g.RDS", temp_file)
Z.s.g <- readRDS(temp_file)
# Z.s.g <- readRDS(fs::path(dir.code,"Z.s.g.RDS"))

#Family
temp_file <- tempfile(fileext = ".RDS") 
s3$download_file(bucket,"Z.g.f.RDS", temp_file)
Z.g.f <- readRDS(temp_file)
# Z.g.f <- readRDS(fs::path(dir.code,"Z.g.f.RDS"))

#Order
temp_file <- tempfile(fileext = ".RDS") 
s3$download_file(bucket,"Z.f.o.RDS", temp_file)
Z.f.o <- readRDS(temp_file)
# Z.f.o <- readRDS(fs::path(dir.code,"Z.f.o.RDS"))

#Class
temp_file <- tempfile(fileext = ".RDS") 
s3$download_file(bucket,"Z.o.c.RDS", temp_file)
Z.o.c <- readRDS(temp_file)
# Z.o.c <- readRDS(fs::path(dir.code,"Z.o.c.RDS"))

#Phylum
temp_file <- tempfile(fileext = ".RDS") 
s3$download_file(bucket,"Z.c.p.RDS", temp_file)
Z.c.p <- readRDS(temp_file)
# Z.c.p <- readRDS(fs::path(dir.code,"Z.c.p.RDS"))

#3. Load Results & Simulation Parameters ----
#Load simulation scenarios
temp_file <- tempfile(fileext = ".rds") 
s3$download_file(bucket,"Simulation_Scenarios_03_12_25.rds", temp_file)
simulation.parameters <- readRDS(temp_file)
# simulation.parameters <- readRDS(fs::path(dir.code,"Simulation_Scenarios_03_12_25.rds"))

#Load all simulatin parameters
temp_file <- tempfile(fileext = ".csv") 
s3$download_file(bucket,"Simulation_Parameters.csv", temp_file)
sim.par.all <- read.csv(temp_file)
# sim.par.all <- read.csv(fs::path(dir.code,"Simulation_Parameters.csv")) %>% 
sim.par.all <- sim.par.all %>% 
  dplyr::select(-X) %>% 
  dplyr::select(c("P.e.causal","P.s.causal.all","N","OR.exposure")) %>% 
  dplyr::rename(P.s.scenario=P.s.causal.all)

#4. Format Results ----
# results.data <- results
format.fxn <- function(results.data){
  results.data <- results.data %>% 
    dplyr::mutate(component=case_when(component=="Count model coefficients" ~ "Means",
                                      component=="Zero-inflation model coefficients" ~ "Probability",
                                      component=="Means" ~ "Means",
                                      component=="Probability" ~ "Probability")) %>% 
    mutate(sig=case_when(pdir > 0.975 ~ "*",
                         pdir <= 0.975 ~ "N.S.",
                         is.na(pdir) ~ sig))
  
  zing.data <- results.data %>%
    filter(model=="ZINB" | model =="Poisson") %>%
    filter(component=="Means") %>%
    mutate(pval = as.numeric(pval)) %>%
    dplyr::group_by(domain) %>%
    mutate(fdr = p.adjust(pval, method = "fdr")) %>%
    ungroup() %>%
    mutate(fdr.sig=ifelse(fdr<0.05,"*","N.S.")) %>%
    mutate(model="ZING") %>% 
    mutate(lcl = estimate - (1.96 * SE),ucl = estimate + (1.96 * SE)
    ) %>% 
    dplyr::select(c("taxa_full","exposure","estimate","lcl","ucl","component","model","domain","Scenario","fdr.sig","iteration"))
  
  # zing.data.prob <- results.data %>%
  #   filter(model=="ZINB" | model =="Poisson") %>%
  #   filter(component=="Probability") %>%
  #   mutate(pval = as.numeric(pval)) %>%
  #   dplyr::group_by(domain) %>%
  #   mutate(fdr = p.adjust(pval, method = "fdr")) %>%
  #   ungroup() %>%
  #   mutate(fdr.sig=ifelse(fdr<0.05,"*","N.S.")) %>%
  #   mutate(model="ZING") 
  
  bayesian <- results.data %>%
    filter(component=="Means") %>%
    filter(model=="BaHZING" | model =="RBaHZING") %>%
    mutate(taxa_full=str_replace_all(taxa_full,"k__species","species")) %>% 
    # mutate(fdr=p) %>% 
    # mutate(fdr.sig=ifelse(fdr<0.05,"*","N.S.")) 
    mutate(fdr.sig=sig) %>% 
    dplyr::rename(lcl=bci_lcl,
                  ucl=bci_ucl) %>% 
    dplyr::select(c("taxa_full","exposure","estimate","lcl","ucl","component","model","domain","Scenario","fdr.sig","iteration"))
  
  # formatted.data.prob <- results.data %>%
  #   filter(component=="Probability") %>%
  #   filter(model=="BaHZING" | model =="RBaHZING") %>%
  #   mutate(taxa_full=str_replace_all(taxa_full,"k__species","species")) %>% 
  #   mutate(fdr=pval) %>% 
  #   mutate(fdr.sig=ifelse(fdr<0.05,"*","N.S.")) 
  
  formatted.data <- rbind(zing.data,  bayesian)
  formatted.data <- formatted.data %>% 
    dplyr::rename(sig=fdr.sig)
  # formatted.data.prob <- rbind(zing.data.prob,formatted.data.prob)
  # rm(zing.data,zing.data.prob)
  rm(zing.data,bayesian)
  
  # formatted.data <- formatted.data %>%
  #   mutate(Taxa.Level=case_when(grepl("species",Taxa) ~ "Species",
  #                               grepl("genus",Taxa) ~ "Genus",
  #                               grepl("family",Taxa) ~ "Family",
  #                               grepl("order",Taxa) ~ "Order",
  #                               grepl("class",Taxa) ~ "Class",
  #                               grepl("phylum",Taxa) ~ "Phylum"))
  # #Remove columns not needed
  # formatted.data <- formatted.data %>%
  #   select(-Component)
  
  # # Define a function to assign replication numbers
  # assign_replication_numbers <- function(group) {
  #   group$Replication <- seq_len(nrow(group))
  #   return(group)
  # }
  
  # # Apply the function to each group
  # formatted.data <- formatted.data %>%
  #   tidylog::group_by(Taxa, Exposure, Model,iteration) %>%
  #   group_modify(~ assign_replication_numbers(.x))
  
  #Group number for plotting
  formatted.data <- formatted.data %>%
    tidylog::group_by(taxa_full,exposure,iteration) %>%
    tidylog::mutate(group_number = cur_group_id()) %>%
    ungroup()
  
  # formatted.data.prob <- formatted.data.prob %>%
  #   tidylog::group_by(taxa_full,exposure,iteration) %>%
  #   tidylog::mutate(group_number = cur_group_id()) %>%
  #   ungroup()
  
  sim.par <- sim.par.all[unique(formatted.data$Scenario),]
  par <- sim.par$P.s.scenario
  
  #Causal Exposures
  if (sim.par$P.e.causal==0){
    causal.e <- NULL
  } else {
    causal.e <- paste0("X.",rep(1:sim.par$P.e.causal))
  }
  
  #Species level indicator
  if(par==1){
    causal.s <- NULL
  }
  if(par!=1){
    causal.s <- paste0("species",simulation.parameters[[par]])
  }
  
  #Genus level indicator
  is_in_causal_s <- which(rownames(Z.s.g) %in% causal.s)
  subset_matrix <- Z.s.g[is_in_causal_s, ]
  genus_occurrences <- colSums(subset_matrix)
  causal.g <- names(which(genus_occurrences > 0))
  
  #Family level indicator
  is_in_causal_g <- which(rownames(Z.g.f) %in% causal.g)
  subset_matrix <- Z.g.f[is_in_causal_g, ]
  family_occurrences <- colSums(subset_matrix)
  causal.f <- names(which(family_occurrences > 0))
  
  #Order level indicator
  is_in_causal_f <- which(rownames(Z.f.o) %in% causal.f)
  subset_matrix <- Z.f.o[is_in_causal_f, ]
  order_occurrences <- colSums(subset_matrix)
  causal.o <- names(which(order_occurrences > 0))
  
  # Class level indicator
  is_in_causal_c <- which(rownames(Z.o.c) %in% causal.o)
  subset_matrix <- Z.o.c[is_in_causal_c, ]
  class_occurrences <- colSums(subset_matrix)
  causal.c <- names(which(class_occurrences > 0))
  
  #Phylum level indicator
  is_in_causal_p<- which(rownames(Z.c.p) %in% causal.c)
  subset_matrix <- Z.c.p[is_in_causal_p, ]
  phylum_occurrences <- colSums(subset_matrix)
  causal.p <- names(which(phylum_occurrences > 0))
  
  #Set the intended OR for the scenario
  setOR <- sim.par$OR.exposure
  null <- 1
  
  #Calculate expected values for Phyla
  ExpectedValues.p <- data.frame()
  phylum.names <- unique(formatted.data$taxa_full[which(grepl("phylum",formatted.data$taxa_full))])
  for (p in phylum.names) {
    causal <- length(which(paste0("class",which(Z.c.p[,p]==1)) %in% causal.c))
    total <- length(which(Z.c.p[,p] == 1))
    Expected.p <- data.frame(taxa_full=p,ExpectedLogOdds=log(setOR) *(causal / total),Index=as.numeric(str_remove(p,"phylum")))
    ExpectedValues.p <- rbind(ExpectedValues.p,Expected.p)
  }
  
  #Calculate expected values for Class
  ExpectedValues.c <- data.frame()
  class.names <- unique(formatted.data$taxa_full[which(grepl("class",formatted.data$taxa_full))])
  for (c in class.names) {
    causal <- length(which(paste0("order",which(Z.o.c[,c]==1)) %in% causal.o))
    total <- length(which(Z.o.c[,c] == 1))
    Expected.c <- data.frame(taxa_full=c,ExpectedLogOdds=log(setOR) *(causal / total),Index=as.numeric(str_remove(c,"class")))
    ExpectedValues.c <- rbind(ExpectedValues.c,Expected.c)
  }
  
  #Calculate expected values for Order
  ExpectedValues.o <- data.frame()
  order.names <- unique(formatted.data$taxa_full[which(grepl("order",formatted.data$taxa_full))])
  for (o in order.names) {
    causal <- length(which(paste0("family",which(Z.f.o[,o]==1)) %in% causal.f))
    total <- length(which(Z.f.o[,o] == 1))
    Expected.o <- data.frame(taxa_full=o,ExpectedLogOdds= log(setOR) *(causal / total),Index=as.numeric(str_remove(o,"order")))
    ExpectedValues.o <- rbind(ExpectedValues.o,Expected.o)
  }
  
  #Calculate expected values for Family
  ExpectedValues.f <- data.frame()
  family.names <- unique(formatted.data$taxa_full[which(grepl("family",formatted.data$taxa_full))])
  for (f in family.names) {
    causal <- length(which(paste0("genus",which(Z.g.f[,f]==1)) %in% causal.g))
    total <- length(which(Z.g.f[,f] == 1))
    Expected.f <- data.frame(taxa_full=f,ExpectedLogOdds=log(setOR) *(causal / total),Index=as.numeric(str_remove(f,"family")))
    ExpectedValues.f <- rbind(ExpectedValues.f,Expected.f)
  }
  
  #Calculate expected values for Genus
  ExpectedValues.g <- data.frame()
  genus.names <- unique(formatted.data$taxa_full[which(grepl("genus",formatted.data$taxa_full))])
  for (g in genus.names) {
    causal <- length(which(paste0("species",which(Z.s.g[,g]==1)) %in% causal.s))
    total <- length(which(Z.s.g[,g] == 1))
    Expected.g <- data.frame(taxa_full=g,ExpectedLogOdds=log(setOR) *(causal / total),Index=as.numeric(str_remove(g,"genus")))
    ExpectedValues.g <- rbind(ExpectedValues.g,Expected.g)
  }
  
  #Calculate expected values for Species - should this be unaffect by other levels?? Check this!!!
  ExpectedValues.s <- data.frame()
  species.names <- unique(formatted.data$taxa_full[which(grepl("species",formatted.data$taxa_full))])
  for (s in species.names) {
    # causal <- length(which(paste0("species",which(Z.s.g[,s]==1)) %in% causal.s))
    # total <- length(which(Z.s.g[,g] == 1))
    Expected.s <- data.frame(taxa_full=s,ExpectedLogOdds=ifelse(s %in% causal.s,log(setOR),0),Index=as.numeric(str_remove(s,"species")))
    ExpectedValues.s <- rbind(ExpectedValues.s,Expected.s)
  }
  
  ExpectedValues <- rbind(ExpectedValues.p,ExpectedValues.c)
  ExpectedValues <- rbind(ExpectedValues,ExpectedValues.o)
  ExpectedValues <- rbind(ExpectedValues,ExpectedValues.f)
  ExpectedValues <- rbind(ExpectedValues,ExpectedValues.g)
  ExpectedValues <- rbind(ExpectedValues,ExpectedValues.s)
  formatted.data <- tidylog::left_join(formatted.data,ExpectedValues,by="taxa_full")
  
  #Now make only the causal exposures causal
  formatted.data <- formatted.data %>%
    tidylog::mutate(expectedLogOdds = case_when(exposure=="Mixture" ~ NA,
                                                # exposure!="Mixture" ~ ExpectedLogOdds))
                                                exposure %in% causal.e ~ ExpectedLogOdds,
                                                !exposure %in% causal.e ~ log(null)))
  #Remove the intermediate ExpectedLogOdds variable
  formatted.data <- formatted.data %>% 
    dplyr::select(-c("ExpectedLogOdds"))
  
  
  # Expected value for mixture =sum of individual effects
  formatted.data <-  formatted.data %>%
    # filter(Exposure!="Mixture") %>% 
    tidylog::group_by(taxa_full,model,iteration) %>%
    tidylog::mutate(mixture_sum = sum(expectedLogOdds,na.rm=T)) %>%
    tidylog::mutate(expectedLogOdds=ifelse(is.na(expectedLogOdds),mixture_sum,expectedLogOdds)) %>%
    ungroup() %>% 
    select(-mixture_sum)
  
  #Create Causal/Non-causal indicator: If expected value > 0 C, else NC
  formatted.data <- formatted.data %>%
    mutate(indicator=ifelse(expectedLogOdds>0,"Associated","Not Associated"))
  
  associated.taxa <- c(causal.s,causal.g,causal.f,causal.o,causal.c,causal.p)
  
  return(formatted.data)
  
}

#Read in & Format all results
all_formatted_results <- data.frame()
temp_file <- tempfile(fileext = ".RDS") 
# for (scenario in 1:9) {
for (scenario in 4:9) {
  s3$download_file(bucket,paste0("scn",scenario,"_iters1_to_1000.RDS"), temp_file)
  results <- readRDS(temp_file)
  formatted.results <- format.fxn(results)
  all_formatted_results <- rbind(all_formatted_results,formatted.results)
}
rm(formatted.results)

#5. Plot Results ----
temp_file <- tempfile(fileext = ".RDS") 
s3$download_file(bucket,"taxa_structure.RDS", temp_file)
names <- readRDS(temp_file)
# names <- readRDS(fs::path(dir.dat,"taxa_structure.RDS"))

##a. Mixture Analysis ----
#Test on scenario 4
# formatted.data <- formatted.scenario4
# rm(formatted.scenario4)

#Select random sample of iterations to create and view plots
# Set seed for reproducibility
set.seed(123)

# Sample a small number of iterations
n_sample_iterations <- 25  # Adjust this number as needed

# Get unique iterations and sample from them
sampled_iterations <- sample(unique(all_formatted_results$iteration), n_sample_iterations)

# Filter your data to sampled iterations
all_formatted_sample <- all_formatted_results %>% 
  filter(iteration %in% sampled_iterations)

# plot.fxn <- fxn(data_ready) {}
# mixture <- formatted.data %>% 
# mixture <- all_formatted_results %>% 
mixture <- all_formatted_sample %>% 
  filter(exposure=="Mixture") 

mixture$model <- factor(mixture$model, levels = c("BaHZING","RBaHZING","ZING"))
mixture$Scenario <- paste0("Scenario ",mixture$Scenario)
# scenario_order <- c("Scenario 4")
# mixture$Scenario <- factor(mixture$Scenario, levels = scenario_order)
taxa_order <- c("Species","Genus","Family","Order","Class","Phylum")
mixture$domain <- factor(mixture$domain, levels = taxa_order)

NA.mix <- mixture %>%
  filter(indicator=="Not Associated")
A.mix <- mixture %>%
  filter(indicator=="Associated")

mix.plot <- ggplot()+
  geom_line(data = NA.mix, alpha = 0.5, aes(model,estimate,group = group_number,color=estimate)) +
  geom_point(data=NA.mix,aes(model,estimate,color=estimate),alpha = 0.5) +
  geom_line(data = A.mix, alpha = 0.5, aes(model,estimate,group = group_number,color=estimate)) +
  geom_point(data=A.mix,aes(model,estimate,color=estimate),alpha = 0.5) +
  geom_errorbar(data=NA.mix,
                aes(model,estimate,group = model),
                stat = "summary",
                fun.data = function(x) {
                  mean_val <- mean(x)
                  sd_val <- sd(x)
                  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                },
                width = 0.1,
                linewidth = 0.5,
                color = "black"
  ) +
  geom_errorbar(data=A.mix,
                aes(model,estimate,group = model),
                stat = "summary",
                fun.data = function(x) {
                  mean_val <- mean(x)
                  sd_val <- sd(x)
                  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                },
                width = 0.1,
                linewidth = 0.5,
                color = "black"
  ) +
  geom_point(data=NA.mix,
             aes(model,estimate,group = model),
             stat = "summary",
             fun = mean,   # Use the mean function to calculate the summary statistic
             size = 2,     # Size of the points
             shape = 16,   # Shape of the points
             color = "black",  # Color of the points
             fill = "black") +   # Fill color of the points
  geom_point(data=A.mix,
             aes(model,estimate,group = model),
             stat = "summary",
             fun = mean,   # Use the mean function to calculate the summary statistic
             size = 2,     # Size of the points
             shape = 16,   # Shape of the points
             color = "black",  # Color of the points
             fill = "black") +
# geom_errorbar(data=NA.mix,
#               aes(model,estimate,group = model),
#               stat = "summary",
#               fun.data = function(x) {
#                 estimate_val <- estimate(x)
#                 sd_val <- sd(x)
#                 data.frame(y = estimate_val, ymin = estimate_val - sd_val, ymax = estimate_val + sd_val)
#               },
#               width = 0.1,
#               linewidth = 0.5,
#               color = "black"
# ) +
# geom_errorbar(data=A.mix,
#               aes(model,estimate,group = model),
#               stat = "summary",
#               fun.data = function(x) {
#                 estimate_val <- estimate(x)
#                 sd_val <- sd(x)
#                 data.frame(y = estimate_val, ymin = estimate_val - sd_val, ymax = estimate_val + sd_val)
#               },
#               width = 0.1,
#               linewidth = 0.5,
#               color = "black"
# ) +
# geom_point(data=NA.mix,
#            aes(model,estimate,group = model),
#            stat = "summary",
#            fun = mean,   # Use the estimate function to calculate the summary statistic
#            size = 2,     # Size of the points
#            shape = 16,   # Shape of the points
#            color = "black",  # Color of the points
#            fill = "black") +   # Fill color of the points
# geom_point(data=A.mix,
#            aes(model,estimate,group = model),
#            stat = "summary",
#            fun = mean,   # Use the estimate function to calculate the summary statistic
#            size = 2,     # Size of the points
#            shape = 16,   # Shape of the points
#            color = "black",  # Color of the points
#            fill = "black") +
# facet_grid(Scenario ~ domain)+
facet_grid(Scenario ~ domain,scales = "free_y")+
theme_bw() +
theme(
  text = element_text(family = "Times", size = 20,color="black"),
  axis.text = element_text(family = "Times", size = 15),
  plot.title = element_text(family = "Times", face = "bold", size = 16),
  plot.subtitle = element_text(family = "Times", size = 15),
  axis.title.x = element_blank(),
  axis.title.y = element_text(family = "Times", size = 20,face="bold"),
  axis.text.x=element_text(angle=45,hjust=1),
  # strip.text.y = element_blank(),
  strip.text.x = element_text(size=15,face="bold"))+

scale_colour_gradientn(colours = c("#660000","#660000","#cc0000","#cc0000","#e69f00","#e69f00","#9a9a9a","#9a9a9a"),
                       # scale_colour_gradientn(colours = c("#9a9a9a","#9a9a9a","#9a9a9a","#9a9a9a","#e69f00","#cc0000","#660000"),
                       values = c(1.0,0.8,0.6,0.4,0.2,0)) +
# breaks = c(-0.2,-0.4,-0.6,0,0.2,0.4,0.6))+
theme(legend.position="left")

mix.plot
png("/home/hhampson/Results/Mixture_Plots_Test2.png",res=300,height=4000,width=4000)
plot(mix.plot)
dev.off()

mix.plot.edit <- ggplot()+
  geom_line(data = NA.mix, alpha = 0.3, 
            aes(model, estimate, group = group_number), 
            color = "gray60") +
  geom_line(data = A.mix, alpha = 0.5, 
            aes(model, estimate, group = group_number), 
            color = "firebrick3") +
  geom_errorbar(data = NA.mix,
                aes(model, estimate, group = model),
                stat = "summary",
                fun.data = function(x) {
                  mean_val <- mean(x)
                  sd_val <- sd(x)
                  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                },
                width = 0.2,
                linewidth = 0.8,
                color = "gray30") +
  geom_errorbar(data = A.mix,
                aes(model, estimate, group = model),
                stat = "summary",
                fun.data = function(x) {
                  mean_val <- mean(x)
                  sd_val <- sd(x)
                  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                },
                width = 0.2,
                linewidth = 0.8,
                color = "firebrick4") +
  geom_point(data = NA.mix,
             aes(model, estimate, group = model),
             stat = "summary",
             fun = mean,
             size = 3,
             shape = 21,
             color = "gray30",
             fill = "white",
             stroke = 1.5) +
  geom_point(data = A.mix,
             aes(model, estimate, group = model),
             stat = "summary",
             fun = mean,
             size = 3,
             shape = 21,
             color = "firebrick4",
             fill = "white",
             stroke = 1.5) +
  facet_grid(Scenario ~ domain, scales = "free_y") +
  theme_bw() +
  theme(
    text = element_text(family = "Times", size = 20, color = "black"),
    axis.text = element_text(family = "Times", size = 15),
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "Times", size = 20, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 15, face = "bold"),
    strip.text.y = element_text(size = 15, face = "bold")) +
  labs(y = "Estimate")
mix.plot.edit 

# ggsave(plot=mix.plot,fs::path(dir.results,"Figure 3.3 Scenario 15.jpeg"),height=4,width=15)
# png(fs::path(dir.results,"Scenario 4 Test.jpeg"),res=300,height=1000,width=3000)
png("/home/hhampson/Results/Mixture_Plots_Test.png",res=300,height=4000,width=4000)
plot(mix.plot.edit)
dev.off()


mix.plot.edit2 <- ggplot()+
  geom_line(data = NA.mix, alpha = 0.5, 
            aes(model, estimate, group = group_number, color = estimate)) +
  geom_point(data = NA.mix, aes(model, estimate), 
             alpha = 0.3, color = "gray50", size = 0.5) +  # Neutral color
  geom_line(data = A.mix, alpha = 0.5, 
            aes(model, estimate, group = group_number, color = estimate)) +
  geom_point(data = A.mix, aes(model, estimate), 
             alpha = 0.3, color = "gray50", size = 0.5) +  # Neutral color
  geom_errorbar(data = NA.mix,
                aes(model, estimate, group = model),
                stat = "summary",
                fun.data = function(x) {
                  mean_val <- mean(x)
                  sd_val <- sd(x)
                  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                },
                width = 0.1,
                linewidth = 0.5,
                color = "black") +
  geom_errorbar(data = A.mix,
                aes(model, estimate, group = model),
                stat = "summary",
                fun.data = function(x) {
                  mean_val <- mean(x)
                  sd_val <- sd(x)
                  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                },
                width = 0.1,
                linewidth = 0.5,
                color = "black") +
  geom_point(data = NA.mix,
             aes(model, estimate, group = model),
             stat = "summary",
             fun = mean,
             size = 3,
             shape = 21,
             color = "black",
             fill = "white",
             stroke = 1) +
  geom_point(data = A.mix,
             aes(model, estimate, group = model),
             stat = "summary",
             fun = mean,
             size = 3,
             shape = 21,
             color = "black",
             fill = "white",
             stroke = 1) +
  facet_grid(Scenario ~ domain, scales = "free_y") +
  theme_bw() +
  theme(
    text = element_text(family = "Times", size = 20, color = "black"),
    axis.text = element_text(family = "Times", size = 15),
    plot.title = element_text(family = "Times", face = "bold", size = 16),
    plot.subtitle = element_text(family = "Times", size = 15),
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "Times", size = 20, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 15, face = "bold"),
    legend.position = "right",  # Move legend to right for clarity
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)) +
  scale_colour_gradientn(
    colours = c("#660000", "#cc0000", "#e69f00", "#9a9a9a"),
    values = scales::rescale(c(0.75, 0.5, 0.25, 0, -0.25, -0.5)),
    name = "Effect\nEstimate",  # Clear legend title
    guide = guide_colorbar(barwidth = 1.5, barheight = 15)) +
  labs(y = "Estimate")

png("/home/hhampson/Results/Mixture_Plots_Test3.png",res=300,height=4000,width=4000)
plot(mix.plot.edit2)
dev.off()

# Add a categorical color variable
# Add a categorical color variable
NA.mix.cat <- NA.mix %>%
  mutate(effect_cat = case_when(
    estimate >= 0.5 ~ "Strong Positive",
    estimate >= 0.25 ~ "Moderate Positive",
    estimate >= 0 ~ "Weak Positive",
    estimate >= -0.25 ~ "Weak Negative",
    TRUE ~ "Negative"
  ))

A.mix.cat <- A.mix %>%
  mutate(effect_cat = case_when(
    estimate >= 0.5 ~ "Strong Positive",
    estimate >= 0.25 ~ "Moderate Positive",
    estimate >= 0 ~ "Weak Positive",
    estimate >= -0.25 ~ "Weak Negative",
    TRUE ~ "Negative"
  ))

# Create the plot with discrete colors
mix.plot.edit3 <- ggplot() +
  geom_line(data = NA.mix.cat, alpha = 0.5, 
            aes(model, estimate, group = group_number, color = effect_cat)) +
  geom_line(data = A.mix.cat, alpha = 0.5, 
            aes(model, estimate, group = group_number, color = effect_cat)) +
  geom_errorbar(data = NA.mix.cat,
                aes(model, estimate, group = model),
                stat = "summary",
                fun.data = function(x) {
                  mean_val <- mean(x)
                  sd_val <- sd(x)
                  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                },
                width = 0.1,
                linewidth = 0.5,
                color = "black") +
  geom_errorbar(data = A.mix.cat,
                aes(model, estimate, group = model),
                stat = "summary",
                fun.data = function(x) {
                  mean_val <- mean(x)
                  sd_val <- sd(x)
                  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                },
                width = 0.1,
                linewidth = 0.5,
                color = "black") +
  geom_point(data = NA.mix.cat,
             aes(model, estimate, group = model),
             stat = "summary",
             fun = mean,
             size = 2,
             shape = 16,
             color = "black",
             fill = "black") +
  geom_point(data = A.mix.cat,
             aes(model, estimate, group = model),
             stat = "summary",
             fun = mean,
             size = 2,
             shape = 16,
             color = "black",
             fill = "black") +
  facet_grid(Scenario ~ domain, scales = "free_y") +
  theme_bw() +
  theme(
    text = element_text(family = "Times", size = 20, color = "black"),
    axis.text = element_text(family = "Times", size = 15),
    plot.title = element_text(family = "Times", face = "bold", size = 16),
    plot.subtitle = element_text(family = "Times", size = 15),
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "Times", size = 20, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 15, face = "bold"),
    legend.position = "right") +
  scale_color_manual(
    name = "Effect Size",
    values = c("Strong Positive" = "#660000",
               "Moderate Positive" = "#cc0000",
               "Weak Positive" = "#e69f00",
               "Weak Negative" = "#cccccc",
               "Null" = "#9a9a9a"),
    breaks = c("Strong Positive", "Moderate Positive", "Weak Positive", 
               "Weak Negative", "Negative")) +
  labs(y = "Estimate")

mix.plot.edit3
png("/home/hhampson/Results/Mixture_Plots_Test4.png",res=300,height=4000,width=4000)
plot(mix.plot.edit3)
dev.off()

##b. Individual Analysis ----
# ind <- formatted.data %>%
#   tidylog::filter(Exposure!="Mixture") 
ind <- all_formatted_sample %>% 
  filter(exposure!="Mixture") 

ind$model <- factor(ind$model, levels = c("BaHZING","RBaHZING","ZING"))
ind$Scenario <- paste0("Scenario ",ind$Scenario)
# scenario_order <- c("Scenario 4")
# ind$Scenario <- factor(ind$Scenario, levels = scenario_order)
taxa_order <- c("Species","Genus","Family","Order","Class","Phylum")
ind$domain <- factor(ind$domain, levels = taxa_order)

NA.ind <- ind %>%
  filter(indicator=="Not Associated")
A.ind <- ind %>%
  filter(indicator=="Associated")


# ind.plot <- ggplot()+
#   geom_line(data = NC, alpha = 0.5, aes(Model,Mean,group = group_number,color=Mean)) +
#   geom_point(data=NC,aes(Model,Mean,color=Mean),alpha = 0.5) +
#   geom_line(data = C, alpha = 0.5, aes(Model,Mean,group = group_number,color=Mean)) +
#   geom_point(data=C,aes(Model,Mean,color=Mean),alpha = 0.5) +
#   geom_errorbar(data=NC,
#                 aes(Model,Mean,group = Model),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.5,
#                 color = "black"
#   ) +
#   geom_errorbar(data=C,
#                 aes(Model,Mean,group = Model),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.5,
#                 color = "black"
#   ) +
#   geom_point(data=NC,
#              aes(Model,Mean,group = Model),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 2,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +   # Fill color of the points
#   geom_point(data=C,
#              aes(Model,Mean,group = Model),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 2,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   facet_grid(Scenario ~ Taxa.Level)+
#   # facet_grid(Scenario ~ Taxa.Level,scales = "free_y", space = "free_y")+
#   theme_bw() +
#   theme(
#     text = element_text(family = "Times", size = 20,color="black"),
#     axis.text = element_text(family = "Times", size = 15),
#     plot.title = element_text(family = "Times", face = "bold", size = 16),
#     plot.subtitle = element_text(family = "Times", size = 15),
#     axis.title.x = element_blank(),
#     axis.title.y = element_text(family = "Times", size = 20,face="bold"),
#     axis.text.x=element_text(angle=45,hjust=1),
#     # strip.text.y = element_blank(),
#     strip.text.x = element_text(size=15,face="bold"))+
#   # labs(x="Log Odds Ratio")+
#   # scale_color_gradient2(low = "gray", mid = "orange", high = "#CC0066")
#   # scale_color_gradient2(low = "gray", mid = c("yellow","orange"), high = "#CC0066",
#   #                       midpoint = max/2, limits = c(min, max), na.value = NA)
#   scale_colour_gradientn(colours = c("#9a9a9a","#9a9a9a","#9a9a9a","#9a9a9a","#e69f00","#cc0000","#660000","#660000"),
#                          # values = c(1.0,0.8,0.6,0.4,0.2,0))+
#                          # c("#660000","#660000","#cc0000","#e69f00","#9a9a9a","#9a9a9a")
#                          # breaks = c(1,0.5,0,-0.5,-1)) +
#                          breaks = c(3,2,1,0,-1,-2,-3)) +
#   theme(legend.position="left")
# 
# ind.plot

# ggsave(plot=ind.plot,fs::path(dir.results,"Figure 3.2 Scenario 15.jpeg"),height=4,width=15)
ind.plot <- ggplot() +
  geom_line(data = NA.ind, alpha = 0.5, 
            aes(model, estimate, group = group_number, color = estimate)) +
  geom_line(data = A.ind, alpha = 0.5, 
            aes(model, estimate, group = group_number, color = estimate)) +
  geom_errorbar(data = NA.ind,
                aes(model, estimate, group = model),
                stat = "summary",
                fun.data = function(x) {
                  mean_val <- mean(x)
                  sd_val <- sd(x)
                  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                },
                width = 0.1,
                linewidth = 0.5,
                color = "black") +
  geom_errorbar(data = A.ind,
                aes(model, estimate, group = model),
                stat = "summary",
                fun.data = function(x) {
                  mean_val <- mean(x)
                  sd_val <- sd(x)
                  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                },
                width = 0.1,
                linewidth = 0.5,
                color = "black") +
  geom_point(data = NA.ind,
             aes(model, estimate, group = model),
             stat = "summary",
             fun = mean,
             size = 2,
             shape = 16,
             color = "black",
             fill = "black") +
  geom_point(data = A.ind,
             aes(model, estimate, group = model),
             stat = "summary",
             fun = mean,
             size = 2,
             shape = 16,
             color = "black",
             fill = "black") +
  facet_grid(Scenario ~ domain, scales = "free_y") +
  theme_bw() +
  theme(
    text = element_text(family = "Times", size = 20, color = "black"),
    axis.text = element_text(family = "Times", size = 15),
    plot.title = element_text(family = "Times", face = "bold", size = 16),
    plot.subtitle = element_text(family = "Times", size = 15),
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "Times", size = 20, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 15, face = "bold"),
    legend.position = "left") +
  scale_colour_gradientn(
    colours = c("#660000", "#660000", "#cc0000", "#cc0000", "#e69f00", "#e69f00", "#9a9a9a", "#9a9a9a"),
    values = c(1.0, 0.8, 0.6, 0.4, 0.2, 0)) +
  labs(y = "Estimate")

ind.plot
png("/home/hhampson/Results/Individual_Plots_Test1.png",res=300,height=4000,width=4000)
plot(ind.plot)
dev.off()

ind.plot.v2 <- ggplot() +
  geom_line(data = NA.ind, alpha = 0.3, 
            aes(model, estimate, group = group_number), 
            color = "gray60") +
  geom_line(data = A.ind, alpha = 0.5, 
            aes(model, estimate, group = group_number), 
            color = "firebrick3") +
  geom_errorbar(data = NA.ind,
                aes(model, estimate, group = model),
                stat = "summary",
                fun.data = function(x) {
                  mean_val <- mean(x)
                  sd_val <- sd(x)
                  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                },
                width = 0.2,
                linewidth = 0.8,
                color = "gray30") +
  geom_errorbar(data = A.ind,
                aes(model, estimate, group = model),
                stat = "summary",
                fun.data = function(x) {
                  mean_val <- mean(x)
                  sd_val <- sd(x)
                  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                },
                width = 0.2,
                linewidth = 0.8,
                color = "firebrick4") +
  geom_point(data = NA.ind,
             aes(model, estimate, group = model),
             stat = "summary",
             fun = mean,
             size = 3,
             shape = 21,
             color = "gray30",
             fill = "white",
             stroke = 1.5) +
  geom_point(data = A.ind,
             aes(model, estimate, group = model),
             stat = "summary",
             fun = mean,
             size = 3,
             shape = 21,
             color = "firebrick4",
             fill = "white",
             stroke = 1.5) +
  facet_grid(Scenario ~ domain, scales = "free_y") +
  theme_bw() +
  theme(
    text = element_text(family = "Times", size = 20, color = "black"),
    axis.text = element_text(family = "Times", size = 15),
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "Times", size = 20, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 15, face = "bold"),
    strip.text.y = element_text(size = 15, face = "bold")) +
  labs(y = "Estimate")

ind.plot.v2
png("/home/hhampson/Results/Individual_Plots_Test2.png",res=300,height=4000,width=4000)
plot(ind.plot.v2)
dev.off()

# Add categorical color variable
NA.ind.cat <- NA.ind %>%
  mutate(effect_cat = case_when(
    estimate >= 0.5 ~ "Strong Positive",
    estimate >= 0.25 ~ "Moderate Positive",
    estimate >= 0 ~ "Weak Positive",
    estimate >= -0.25 ~ "Weak Negative",
    TRUE ~ "Negative"
  )) %>%
  mutate(effect_cat = factor(effect_cat, 
                             levels = c("Strong Positive", 
                                        "Moderate Positive", 
                                        "Weak Positive", 
                                        "Weak Negative", 
                                        "Negative")))

A.ind.cat <- A.ind %>%
  mutate(effect_cat = case_when(
    estimate >= 0.5 ~ "Strong Positive",
    estimate >= 0.25 ~ "Moderate Positive",
    estimate >= 0 ~ "Weak Positive",
    estimate >= -0.25 ~ "Weak Negative",
    TRUE ~ "Negative"
  )) %>%
  mutate(effect_cat = factor(effect_cat, 
                             levels = c("Strong Positive", 
                                        "Moderate Positive", 
                                        "Weak Positive", 
                                        "Weak Negative", 
                                        "Negative")))

ind.plot.v3 <- ggplot() +
  geom_line(data = NA.ind.cat, alpha = 0.5, 
            aes(model, estimate, group = group_number, color = effect_cat)) +
  geom_line(data = A.ind.cat, alpha = 0.5, 
            aes(model, estimate, group = group_number, color = effect_cat)) +
  geom_errorbar(data = NA.ind.cat,
                aes(model, estimate, group = model),
                stat = "summary",
                fun.data = function(x) {
                  mean_val <- mean(x)
                  sd_val <- sd(x)
                  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                },
                width = 0.1,
                linewidth = 0.5,
                color = "black") +
  geom_errorbar(data = A.ind.cat,
                aes(model, estimate, group = model),
                stat = "summary",
                fun.data = function(x) {
                  mean_val <- mean(x)
                  sd_val <- sd(x)
                  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                },
                width = 0.1,
                linewidth = 0.5,
                color = "black") +
  geom_point(data = NA.ind.cat,
             aes(model, estimate, group = model),
             stat = "summary",
             fun = mean,
             size = 2,
             shape = 16,
             color = "black",
             fill = "black") +
  geom_point(data = A.ind.cat,
             aes(model, estimate, group = model),
             stat = "summary",
             fun = mean,
             size = 2,
             shape = 16,
             color = "black",
             fill = "black") +
  facet_grid(Scenario ~ domain, scales = "free_y") +
  theme_bw() +
  theme(
    text = element_text(family = "Times", size = 20, color = "black"),
    axis.text = element_text(family = "Times", size = 15),
    plot.title = element_text(family = "Times", face = "bold", size = 16),
    plot.subtitle = element_text(family = "Times", size = 15),
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "Times", size = 20, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 15, face = "bold"),
    legend.position = "right") +
  scale_color_manual(
    name = "Effect Size",
    values = c("Strong Positive" = "#660000",
               "Moderate Positive" = "#cc0000",
               "Weak Positive" = "#e69f00",
               "Weak Negative" = "#cccccc",
               "Negative" = "#9a9a9a"),
    breaks = c("Strong Positive", "Moderate Positive", "Weak Positive", 
               "Weak Negative", "Negative")) +
  labs(y = "Estimate")

ind.plot.v3
png("/home/hhampson/Results/Individual_Plots_Test3.png",res=300,height=4000,width=4000)
plot(ind.plot.v3)
dev.off()

ind.plot.v4 <- ggplot() +
  geom_line(data = NA.ind, alpha = 0.5, 
            aes(model, estimate, group = group_number, color = estimate)) +
  geom_point(data = NA.ind, aes(model, estimate), 
             alpha = 0.3, color = "gray50", size = 0.5) +
  geom_line(data = A.ind, alpha = 0.5, 
            aes(model, estimate, group = group_number, color = estimate)) +
  geom_point(data = A.ind, aes(model, estimate), 
             alpha = 0.3, color = "gray50", size = 0.5) +
  geom_errorbar(data = NA.ind,
                aes(model, estimate, group = model),
                stat = "summary",
                fun.data = function(x) {
                  mean_val <- mean(x)
                  sd_val <- sd(x)
                  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                },
                width = 0.1,
                linewidth = 0.5,
                color = "black") +
  geom_errorbar(data = A.ind,
                aes(model, estimate, group = model),
                stat = "summary",
                fun.data = function(x) {
                  mean_val <- mean(x)
                  sd_val <- sd(x)
                  data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
                },
                width = 0.1,
                linewidth = 0.5,
                color = "black") +
  geom_point(data = NA.ind,
             aes(model, estimate, group = model),
             stat = "summary",
             fun = mean,
             size = 3,
             shape = 21,
             color = "black",
             fill = "white",
             stroke = 1.5) +
  geom_point(data = A.ind,
             aes(model, estimate, group = model),
             stat = "summary",
             fun = mean,
             size = 3,
             shape = 21,
             color = "black",
             fill = "white",
             stroke = 1.5) +
  facet_grid(Scenario ~ domain, scales = "free_y") +
  theme_bw() +
  theme(
    text = element_text(family = "Times", size = 20, color = "black"),
    axis.text = element_text(family = "Times", size = 15),
    axis.title.x = element_blank(),
    axis.title.y = element_text(family = "Times", size = 20, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 15, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)) +
  scale_colour_gradientn(
    colours = c("#660000", "#cc0000", "#e69f00", "#9a9a9a"),
    values = scales::rescale(c(0.75, 0.5, 0.25, 0, -0.25, -0.5)),
    name = "Effect\nEstimate",
    guide = guide_colorbar(barwidth = 1.5, barheight = 15)) +
  labs(y = "Estimate")

ind.plot.v4
png("/home/hhampson/Results/Individual_Plots_Test4.png",res=300,height=4000,width=4000)
plot(ind.plot.v4)
dev.off()

#7. Sensitivity Measures ----
#GO back and make bahzing and ridge pvalues based on bcis not pvals
##A. Individual----
individual <- all_formatted_results %>% 
  filter(exposure!="Mixture")
data <- individual
sens.fxn.old <- function(data){
  BaH_ZING <- data %>% 
    filter(Model=="BaH-ZING") %>% 
    filter(Exposure!="Mixture") %>%
    # filter(grepl("species",Taxa)) %>% 
    mutate(PosNeg = case_when(indicator=="Causal" & P_Value=="*" ~ "True Positive",
                              indicator=="Causal" & P_Value!="*" ~ "False Negative",
                              indicator=="Non-Causal" & P_Value=="*" ~ "False Positive",
                              indicator=="Non-Causal" & P_Value!="*" ~ "True Negative")) %>% 
    group_by(Taxa.Level,Scenario) %>%
    # group_by(Taxa.Level) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity
    ) %>%  
    mutate(Model="BaH-ZING") %>% 
    ungroup()
  # select(-Specificity)
  
  measures <- c("Specificity","Sensitivity","PPV","NPV","FDR","TypeI","TypeII")
  
  
  ZING <- data %>% 
    # filter(Model=="ZINB" | Model=="Poisson") %>% 
    filter(Model=="ZING") %>% 
    filter(Exposure!="Mixture") %>%
    # filter(grepl("species",Taxa)) %>% 
    mutate(PosNeg = case_when(indicator=="Causal" & P_Value=="*" ~ "True Positive",
                              indicator=="Causal" & P_Value!="*" ~ "False Negative",
                              indicator=="Non-Causal" & P_Value=="*" ~ "False Positive",
                              indicator=="Non-Causal" & P_Value!="*" ~ "True Negative")) %>% 
    group_by(Taxa.Level,Scenario) %>%
    # group_by(Taxa.Level) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity
    ) %>% 
    mutate(Model="ZING") %>% 
    ungroup()
  # select(-Specificity)
  
  ZING.fdr <- data %>% 
    # filter(Model=="ZINB" | Model=="Poisson") %>% 
    filter(Model=="ZING") %>% 
    filter(Exposure!="Mixture") %>%
    # filter(grepl("species",Taxa)) %>% 
    mutate(PosNeg = case_when(indicator=="Causal" & fdr.p=="*" ~ "True Positive",
                              indicator=="Causal" & fdr.p!="*" ~ "False Negative",
                              indicator=="Non-Causal" & fdr.p=="*" ~ "False Positive",
                              indicator=="Non-Causal" & fdr.p!="*" ~ "True Negative")) %>% 
    group_by(Taxa.Level,Scenario) %>%
    # group_by(Taxa.Level) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity
    ) %>% 
    mutate(Model="Adj.ZING") %>% 
    ungroup() 
  # select(-Specificity)
  
  # SensSpec <- rbind(BaH_ZING,RBaH_ZING)
  SensSpec <- rbind(BaH_ZING,ZING)
  # SensSpec <- rbind(SensSpec,ZING)
  SensSpec <- rbind(SensSpec,ZING.fdr)
  # SensSpec <- rbind(SensSpec,ZING.fdr)
  # SensSpec.s <- SensSpec %>%
  #   filter(Taxa.Level=="Species")
  # SensSpec.g <- SensSpec %>%
  #   filter(Taxa.Level=="Genus")
  # SensSpec.f <- SensSpec %>%
  #   filter(Taxa.Level=="Family")
  # SensSpec.o <- SensSpec %>%
  #   filter(Taxa.Level=="Order")
  # SensSpec.c <- SensSpec %>%
  #   filter(Taxa.Level=="Class")
  # SensSpec.p <- SensSpec %>%
  #   filter(Taxa.Level=="Phylum")
  
  # B <- SensSpec %>%
  #   filter(Model=="BaH-ZING")
  # mean(B$Specificity,na.rm=T)#0.998
  # range(B$Specificity,na.rm=T)
  # mean(B$PPV,na.rm=T)#0.75
  # mean(B$Sensitivity,na.rm=T) #0.37
  # mean(B$TypeI,na.rm=T)*100 #0.19%
  # mean(B$TypeII,na.rm=T)*100 #63%
  # 
  # R <- SensSpec.s %>% 
  #   filter(Model=="Ridge ZING")
  # mean(R$Specificity,na.rm=T) #0.995
  # range(R$Specificity,na.rm=T)
  # mean(R$PPV,na.rm=T) #0.85
  # mean(R$Sensitivity,na.rm=T) #100
  # mean(R$TypeI,na.rm=T)*100 #0.48% 
  # mean(R$TypeII,na.rm=T)*100 #0.024% 
  # 
  # Z <- SensSpec %>% 
  #   filter(Model=="ZING")
  # mean(Z$Specificity,na.rm=T)
  # range(Z$Specificity,na.rm=T)
  # mean(Z$PPV,na.rm=T)
  # mean(Z$Sensitivity,na.rm=T)
  # mean(Z$TypeI,na.rm=T)*100 #25% 
  # mean(Z$TypeII,na.rm=T)*100 #0.006%
  # 
  # Z.a <- SensSpec %>% 
  #   filter(Model=="Adj.ZING")
  # mean(Z.a$Specificity,na.rm=T)
  # range(Z.a$Specificity,na.rm=T)
  # mean(Z.a$PPV,na.rm=T)
  # mean(Z.a$Sensitivity,na.rm=T)
  # mean(Z.a$TypeI,na.rm=T)*100 #20% 
  # mean(Z.a$TypeII,na.rm=T)*100 #0.226%
  
  # taxa_order <- c("Species","Genus","Family","Order","Class","Phylum")
  # SensSpec$Taxa.Level <- factor(SensSpec$Taxa.Level, levels = taxa_order)
  # write.csv(SensSpec,fs::path(dir.results,"SensitivitySpecificity_All.csv"))
  measures <- c("Specificity","Sensitivity","PPV","NPV","FDR","TypeI","TypeII")
  # measures <- c("PPV","FDR","TypeI")
  # SensSpec.Sum <-SensSpec %>% 
  #   tidylog::select(-c(Taxa.Level,Scenario)) %>% 
  #   tidylog::group_by(Model) %>% 
  #   tidylog::summarize(across(all_of(measures),~mean(.,na.rm=T))) 
  
  
  #   return(SensSpec)
  # }
  
  #Bias non causal and causal...and add MSE? Variance of mean?
  # Melt the data for ggplot2
  # SensSpec2 <- SensSpec %>%
  #   select(-FDR2) 
  # select(-Scenario)
  # select(FDR,Model,Scenario:Type2ErrorRate) %>% 
  # mutate(Scenario=paste0("Scenario ",Scenario))
  # SensitivitySpecificity2$FDR[1:4] <- NA
  # scenario.order <- c("Scenario 1", "Scenario 2", "Scenario 3","Scenario 4",
  #                     "Scenario 5","Scenario 6","Scenario 7","Scenario 8","Scenario 11")
  # Convert Scenario to a factor with the desired order
  # SensitivitySpecificity2$Scenario <- factor(SensitivitySpecificity2$Scenario, levels = scenario.order)
  
  # melted_data <- reshape2::melt(SensSpec2, id.vars = c('Taxa.Level', 'Model','Scenario'))
  melted_data <- reshape2::melt(SensSpec, id.vars = c('Model','Taxa.Level','Scenario'))
  # melted_data <- melted_data %>% 
  # filter(Scenario!="Scenario 11")
  # model_order <- c('Ridge ZING', 'BaH-ZING', 'ZING','Adj.ZING')
  model_order <- c('BaH-ZING', 'ZING','Adj.ZING')
  melted_data$Model <- factor(melted_data$Model, levels = model_order)
  taxa_order <- c("Species","Genus","Family","Order","Class","Phylum")
  # melted_data$Taxa.Level <- factor(melted_data$Taxa.Level, levels = taxa_order)
  
  
  #Order Taxa levels
  # taxa_order
  melted_data$Taxa.Level <- factor(melted_data$Taxa.Level, levels = taxa_order)
  return(melted_data)
}
sens.fxn <- function(data){
  BaH_ZING <- data %>% 
    filter(model=="BaHZING") %>% 
    filter(exposure!="Mixture") %>%
    mutate(PosNeg = case_when(indicator=="Associated" & sig=="*" ~ "True Positive",
                              indicator=="Associated" & sig!="*" ~ "False Negative",
                              indicator=="Not Associated" & sig=="*" ~ "False Positive",
                              indicator=="Not Associated" & sig!="*" ~ "True Negative")) %>% 
    group_by(domain, Scenario) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity,
      .groups = 'drop'
    ) %>%  
    mutate(Model="BaHZING")
  
  RBaH_ZING <- data %>% 
    filter(model=="RBaHZING") %>% 
    filter(exposure!="Mixture") %>%
    mutate(PosNeg = case_when(indicator=="Associated" & significant=="*" ~ "True Positive",
                              indicator=="Associated" & significant!="*" ~ "False Negative",
                              indicator=="Not Associated" & significant=="*" ~ "False Positive",
                              indicator=="Not Associated" & significant!="*" ~ "True Negative")) %>% 
    group_by(domain, Scenario) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity,
      .groups = 'drop'
    ) %>% 
    mutate(Model="RBaHZING")
  
  ZING <- data %>% 
    filter(model=="ZING") %>% 
    filter(exposure!="Mixture") %>%
    mutate(PosNeg = case_when(indicator=="Associated" & significant=="*" ~ "True Positive",
                              indicator=="Associated" & significant!="*" ~ "False Negative",
                              indicator=="Not Associated" & significant=="*" ~ "False Positive",
                              indicator=="Not Associated" & significant!="*" ~ "True Negative")) %>% 
    group_by(domain, Scenario) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity,
      .groups = 'drop'
    ) %>% 
    mutate(Model="ZING")
  
  ZING.fdr <- data %>% 
    filter(model=="ZING") %>% 
    filter(exposure!="Mixture") %>%
    mutate(PosNeg = case_when(indicator=="Associated" & fdr.significant=="*" ~ "True Positive",
                              indicator=="Associated" & fdr.significant!="*" ~ "False Negative",
                              indicator=="Not Associated" & fdr.significant=="*" ~ "False Positive",
                              indicator=="Not Associated" & fdr.significant!="*" ~ "True Negative")) %>% 
    group_by(domain, Scenario) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity,
      .groups = 'drop'
    ) %>% 
    mutate(Model="Adj.ZING")
  
  SensSpec <- rbind(BaH_ZING, RBaH_ZING, ZING, ZING.fdr)
  
  measures <- c("Specificity","Sensitivity","PPV","NPV","FDR","TypeI","TypeII")
  
  melted_data <- reshape2::melt(SensSpec, id.vars = c('Model','domain','Scenario'))
  
  model_order <- c('BaHZING', 'RBaHZING', 'ZING','Adj.ZING')
  melted_data$Model <- factor(melted_data$Model, levels = model_order)
  
  taxa_order <- c("Species","Genus","Family","Order","Class","Phylum")
  melted_data$domain <- factor(melted_data$domain, levels = taxa_order)
  
  return(melted_data)
}
formatted.all <- sens.fxn(individual)
rm(individual)
formatted.all <- formatted.all %>% 
  mutate(Scenario=paste0("Scenario ",Scenario))
formatted.all$Scenario <- factor(formatted.all$Scenario,levels=scenario_order)

formatted.wide <- formatted.all %>% 
  pivot_wider(names_from = Taxa.Level,
              values_from=value) %>% 
  select(-Scenario)
formatted.wide <- formatted.wide[c("Model","variable","Species","Genus",
                                   "Family","Order","Class","Phylum")]
write.csv(formatted.wide,fs::path(dir.results,"Sensitivity_Individual_Scenario15.csv"))

##B. Mixture----
mixture.sens <- formatted.data %>% 
  filter(Exposure=="Mixture")
data <- mixture.sens
sens.fxn.mix <- function(data){
  BaH_ZING <- data %>% 
    filter(Model=="BaH-ZING") %>% 
    filter(Exposure=="Mixture") %>%
    # filter(grepl("species",Taxa)) %>% 
    mutate(PosNeg = case_when(indicator=="Causal" & P_Value=="*" ~ "True Positive",
                              indicator=="Causal" & P_Value!="*" ~ "False Negative",
                              indicator=="Non-Causal" & P_Value=="*" ~ "False Positive",
                              indicator=="Non-Causal" & P_Value!="*" ~ "True Negative")) %>% 
    group_by(Taxa.Level,Scenario) %>%
    # group_by(Taxa.Level) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      # NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      # FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity
    ) %>%  
    mutate(Model="BaH-ZING") %>% 
    ungroup()
  
  measures <- c("Specificity","Sensitivity","PPV","NPV","FDR","TypeI","TypeII")
  # measures <- c("PPV","FDR")
  
  # RBaH_ZING <- data %>% 
  #   filter(Model=="Ridge ZING") %>%
  #   filter(Exposure!="Mixture") %>%
  #   # filter(grepl("species",Taxa)) %>% 
  #   mutate(PosNeg = case_when(indicator=="Causal" & P_Value=="*" ~ "True Positive",
  #                             indicator=="Causal" & P_Value!="*" ~ "False Negative",
  #                             indicator=="Non-Causal" & P_Value=="*" ~ "False Positive",
  #                             indicator=="Non-Causal" & P_Value!="*" ~ "True Negative")) %>% 
  #   group_by(Taxa.Level,Scenario) %>%
  #   # group_by(Taxa.Level) %>%
  #   summarise(
  #     Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
  #     Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
  #     PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
  #     # NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
  #     FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
  #     # FDR2 = (1 - PPV)*100,
  #     TypeI = (1-Specificity),
  #     TypeII = 1-Sensitivity
  #   ) %>% 
  #   mutate(across(all_of(measures),~ifelse(Taxa.Level!="Species",NA,.))) %>% 
  #   mutate(Model="Ridge ZING") %>% 
  #   ungroup()
  
  ZING <- data %>% 
    # filter(Model=="ZINB" | Model=="Poisson") %>% 
    filter(Model=="ZING") %>% 
    filter(Exposure=="Mixture") %>%
    # filter(grepl("species",Taxa)) %>% 
    mutate(PosNeg = case_when(indicator=="Causal" & P_Value=="*" ~ "True Positive",
                              indicator=="Causal" & P_Value!="*" ~ "False Negative",
                              indicator=="Non-Causal" & P_Value=="*" ~ "False Positive",
                              indicator=="Non-Causal" & P_Value!="*" ~ "True Negative")) %>% 
    group_by(Taxa.Level,Scenario) %>%
    # group_by(Taxa.Level) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      # NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      # FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity
    ) %>% 
    mutate(Model="ZING") %>% 
    ungroup()
  
  ZING.fdr <- data %>% 
    # filter(Model=="ZINB" | Model=="Poisson") %>% 
    filter(Model=="ZING") %>% 
    filter(Exposure=="Mixture") %>%
    # filter(grepl("species",Taxa)) %>% 
    mutate(PosNeg = case_when(indicator=="Causal" & fdr.p=="*" ~ "True Positive",
                              indicator=="Causal" & fdr.p!="*" ~ "False Negative",
                              indicator=="Non-Causal" & fdr.p=="*" ~ "False Positive",
                              indicator=="Non-Causal" & fdr.p!="*" ~ "True Negative")) %>% 
    group_by(Taxa.Level,Scenario) %>%
    # group_by(Taxa.Level) %>%
    summarise(
      Specificity = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="True Negative")) + length(which(PosNeg=="False Positive")))),
      Sensitivity = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Negative")))),
      PPV = (length(which(PosNeg=="True Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      # NPV = (length(which(PosNeg=="True Negative")) / (length(which(PosNeg=="False Negative")) + length(which(PosNeg=="True Negative")))),
      FDR =  (length(which(PosNeg=="False Positive")) / (length(which(PosNeg=="True Positive")) + length(which(PosNeg=="False Positive")))),
      # FDR2 = (1 - PPV)*100,
      TypeI = (1-Specificity),
      TypeII = 1-Sensitivity
    ) %>% 
    mutate(Model="Adj.ZING") %>% 
    ungroup()
  
  # SensSpec <- rbind(BaH_ZING,RBaH_ZING)
  SensSpec <- rbind(BaH_ZING,ZING)
  # SensSpec <- rbind(SensSpec,ZING)
  SensSpec <- rbind(SensSpec,ZING.fdr)
  
  # melted_data <- reshape2::melt(SensSpec2, id.vars = c('Taxa.Level', 'Model','Scenario'))
  melted_data <- reshape2::melt(SensSpec, id.vars = c('Model','Taxa.Level','Scenario'))
  # melted_data <- melted_data %>% 
  # filter(Scenario!="Scenario 11")
  model_order <- c( 'BaH-ZING', 'ZING','Adj.ZING')
  # model_order <- c('BaH-ZING', 'ZING','Adj.ZING')
  melted_data$Model <- factor(melted_data$Model, levels = model_order)
  taxa_order <- c("Species","Genus","Family","Order","Class","Phylum")
  # melted_data$Taxa.Level <- factor(melted_data$Taxa.Level, levels = taxa_order)
  
  melted_data$Taxa.Level <- factor(melted_data$Taxa.Level, levels = taxa_order)
  return(melted_data)
}
formatted.all.mix <- sens.fxn.mix(mixture.sens)

formatted.all.mix <- formatted.all.mix %>% 
  mutate(Scenario=paste0("Scenario ",Scenario))
formatted.all.mix$Scenario <- factor(formatted.all.mix$Scenario,levels=scenario_order)

formatted.all.mix.wide <- formatted.all.mix %>%  
  pivot_wider(names_from = Taxa.Level,
              values_from=value) %>% 
  select(-Scenario)
formatted.all.mix.wide<- formatted.all.mix.wide[c("Model","variable","Species","Genus",
                                                  "Family","Order","Class","Phylum")]
write.csv(formatted.wide,fs::path(dir.results,"Sensitivity_Mixture_Scenario15.csv"))


#8. Dendrograms ----
##A. Individual ----
#Scenarios
means <- formatted.data %>%
  tidylog::filter(Exposure!="Mixture") %>%
  tidylog::filter(Model!="Ridge ZING")  
# tidylog::filter(Model!="ZING")
means2 <- means

source("Microbiome_Cleaning.R")

key <- species.key %>%
  select(species,sim.num) %>%
  rename(Taxa=sim.num) %>%
  rename(taxa.name=species)
genus.key <- genus.key %>%
  rename(Taxa=sim.num,
         taxa.name=genus)
key <- rbind(key,genus.key)
family.key <- family.key %>%
  rename(Taxa=sim.num,
         taxa.name=family)
key <- rbind(key,family.key)
order.key <- order.key %>%
  rename(Taxa=sim.num,
         taxa.name=order)
key <- rbind(key,order.key)
class.key <- class.key %>%
  rename(Taxa=sim.num,
         taxa.name=class)
key <- rbind(key,class.key)
phylum.key <- phylum.key %>%
  rename(Taxa=sim.num,
         taxa.name=phylum)
key <- rbind(key,phylum.key)


means2 <- tidylog::left_join(means2,key,by="Taxa")
#Format for tree
# species.key$species.new
means2$taxa.name.new <- str_replace(means2$taxa.name,"_p_","--p_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_c_","--c_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_o_","--o_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_f_","--f_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_g_","--g_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_s_","--s_")

means2$Scenario <- paste0("Scenario ",means2$Scenario)
scenario14.dataA <- means2 %>%
  filter(Scenario=="Scenario 14") %>% 
  filter(Exposure=="X.1" |Exposure=="X.2" | Exposure =="X.3") %>% 
  filter(Model=="BaH-ZING") %>% 
  # filter(Exposure=="X.1") %>% 
  group_by(taxa.name.new) %>% 
  tidylog::summarize(Mean.total = mean(Mean)) %>% 
  filter(grepl("p__Bacteroidota",taxa.name.new)) %>% 
  mutate(mean.color = case_when(Mean.total<=0 ~ "<= 0",
                                Mean.total>0 & Mean.total <=0.2 ~ "0 < Mean <= 0.2",
                                Mean.total>0.2 & Mean.total <0.4 ~ "0.2 < Mean < 0.4",
                                Mean.total>=0.4 ~ "Mean >= 0.4")) %>% 
  mutate(mean.size= case_when(Mean.total<=0 ~ "<= 0",
                              Mean.total>0 & Mean.total <=0.2 ~ "0 < Mean <= 0.2",
                              Mean.total>0.2 & Mean.total <0.4 ~ "0.2 < Mean < 0.4",
                              Mean.total>=0.4 ~ "Mean >= 0.4")) %>% 
  mutate(mean.color = factor(mean.color,levels=c("<= 0","0 < Mean <= 0.2","0.2 < Mean < 0.4","Mean >= 0.4"))) %>% 
  mutate(mean.size = factor(mean.size,levels=c("<= 0","0 < Mean <= 0.2","0.2 < Mean < 0.4","Mean >= 0.4"))) %>% 
  mutate(taxa.name.new = ifelse(taxa.name.new=="d__Bacteria--p__Bacteroidota","Bacteria",taxa.name.new))

scenario15.dataA <- means2 %>%
  filter(Scenario=="Scenario 15") %>% 
  filter(Model=="BaH-ZING") %>% 
  # filter(Exposure=="X.1") %>% 
  group_by(taxa.name.new) %>% 
  tidylog::summarize(Mean.total = mean(Mean)) %>% 
  filter(grepl("p__Bacteroidota",taxa.name.new)) %>% 
  mutate(mean.color = case_when(Mean.total<=0 ~ "<= 0",
                                Mean.total>0 & Mean.total <=0.2 ~ "0 < Mean <= 0.2",
                                Mean.total>0.2 & Mean.total <0.4 ~ "0.2 < Mean < 0.4",
                                Mean.total>=0.4 ~ "Mean >= 0.4")) %>% 
  mutate(mean.size= case_when(Mean.total<=0 ~ "<= 0",
                              Mean.total>0 & Mean.total <=0.2 ~ "0 < Mean <= 0.2",
                              Mean.total>0.2 & Mean.total <0.4 ~ "0.2 < Mean < 0.4",
                              Mean.total>=0.4 ~ "Mean >= 0.4")) %>% 
  mutate(mean.color = factor(mean.color,levels=c("<= 0","0 < Mean <= 0.2","0.2 < Mean < 0.4","Mean >= 0.4"))) %>% 
  mutate(mean.size = factor(mean.size,levels=c("<= 0","0 < Mean <= 0.2","0.2 < Mean < 0.4","Mean >= 0.4"))) %>% 
  # mutate(taxa.name.new = str_replace(taxa.name.new,"d__Bacteria--p__Bacteroidota","Bacteria")) 
  mutate(taxa.name.new = ifelse(taxa.name.new=="d__Bacteria--p__Bacteroidota","Bacteria",taxa.name.new))

scen14.b <- ggtree(tree2,layout="fan",open.angle=230) %<+% scenario14.dataA +
  geom_point(aes(color=mean.color, size = mean.size), alpha=1)+
  scale_size_manual(values = c("<= 0" = 1, "0 < Mean <= 0.2" = 2, "0.2 < Mean < 0.4" = 3, "Mean >= 0.4" = 4))+
  scale_color_manual(values = c("<= 0"="#9a9a9a", "0 < Mean <= 0.2"="#e69f00", "0.2 < Mean < 0.4"="#cc0000", "Mean >= 0.4"="#660000"))+
  # scale_size_manual(values = c("small" = 1, "small.medium" = 2, "medium" = 3, "large" = 4))+
  # scale_color_manual(values = c("#9a9a9a"="#9a9a9a", "#e69f00"="#e69f00", "#cc0000"="#cc0000", "#660000"="#660000"))
  # scale_color_gradient2(low = "gray", mid = c("gray","yellow","orange","#CC0066"), 
  #                       high = "#CC0066", midpoint = 0.2, limits = c(min, max), na.value = NA) +
  theme(legend.position = "none")


scen15.b <- ggtree(tree2,layout="fan",open.angle=230) %<+% scenario15.dataA +
  geom_point(aes(color=mean.color, size = mean.size), alpha=1)+
  scale_size_manual(values = c("<= 0" = 1, "0 < Mean <= 0.2" = 2, "0.2 < Mean < 0.4" = 3, "Mean >= 0.4" = 4))+
  scale_color_manual(values = c("<= 0"="#9a9a9a", "0 < Mean <= 0.2"="#e69f00", "0.2 < Mean < 0.4"="#cc0000", "Mean >= 0.4"="#660000"))+
  
  # scen15 <- ggtree(tree2,layout="fan",open.angle=230) %<+% scenario15.dataA +
  #   geom_point(aes(color=Mean.total, size = Mean.total), alpha=1)+
  #   scale_colour_gradientn(colours = c("#660000","#cc0000","#cc0000","#e69f00","#9a9a9a"),
  #                          values = c(1.0,0.8,0.6,0.4,0.2,0))
  # scale_color_gradient2(low = "gray", mid = c("gray","yellow","orange","#CC0066"),
  #                       high = "#CC0066", midpoint = 0.2, limits = c(min, max), na.value = NA) +
  theme(legend.position = "none")

scen15.b
ggsave(plot=scen15.b,fs::path(dir.results,"Dendrogram_Individual_Scenario15.jpeg"))

##B. Mixture----
mixture2 <- formatted.data %>% 
  filter(Exposure=="Mixture") %>% 
  filter(Model=="BaH-ZING")
# mixture2 <- mixture 
# filter(Replication==1)
means2 <- tidylog::left_join(mixture2,key,by="Taxa")
#Format for tree
# species.key$species.new
means2$taxa.name.new <- str_replace(means2$taxa.name,"_p_","--p_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_c_","--c_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_o_","--o_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_f_","--f_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_g_","--g_")
means2$taxa.name.new <- str_replace(means2$taxa.name.new,"_s_","--s_")

#Scenario 1
# scenario1.data <- scenario.fxn(scenario1)
#Scenario 2
# scenario2.data <- scenario.fxn(scenario2)
means2$Scenario <- paste0("Scenario ",means2$Scenario)
scenario12.dataA <- means2 %>%
  filter(Scenario=="Scenario 12") %>% 
  filter(Model=="BaH-ZING") %>% 
  # filter(Exposure=="X.1") %>% 
  group_by(taxa.name.new) %>%
  tidylog::summarize(Mean.total = mean(Mean)) %>% 
  mutate(mean.color = case_when(Mean.total<=0 ~ "<= 0",
                                Mean.total>0 & Mean.total <=1 ~ "0 < Mean <= 1.0",
                                Mean.total>1 & Mean.total <2 ~ "1.0 < Mean < 2.0",
                                Mean.total>=2 ~ "Mean >= 2.0")) %>% 
  mutate(mean.size= case_when(Mean.total<=0 ~ "<= 0",
                              Mean.total>0 & Mean.total <=1 ~ "0 < Mean <= 1.0",
                              Mean.total>1 & Mean.total <2 ~ "1.0 < Mean < 2.0",
                              Mean.total>=2 ~ "Mean >= 2.0")) %>% 
  mutate(mean.color = factor(mean.color,levels=c("<= 0","0 < Mean <= 1.0","1.0 < Mean < 2.0","Mean >= 2.0"))) %>% 
  mutate(mean.size = factor(mean.size,levels=c("<= 0","0 < Mean <= 1.0","1.0 < Mean < 2.0","Mean >= 2.0"))) %>% 
  mutate(taxa.name.new = ifelse(taxa.name.new=="d__Bacteria--p__Bacteroidota","Bacteria",taxa.name.new))

# select(taxa.name.new,Mean)
#Scenario 13
scenario13.dataA <- means2 %>%
  filter(Scenario=="Scenario 13") %>% 
  filter(Model=="BaH-ZING") %>% 
  # filter(Exposure=="X.1") %>% 
  group_by(taxa.name.new) %>%
  tidylog::summarize(Mean.total = mean(Mean)) %>% 
  mutate(mean.color = case_when(Mean.total<=0 ~ "<= 0",
                                Mean.total>0 & Mean.total <=1 ~ "0 < Mean <= 1.0",
                                Mean.total>1 & Mean.total <2 ~ "1.0 < Mean < 2.0",
                                Mean.total>=2 ~ "Mean >= 2.0")) %>% 
  mutate(mean.size= case_when(Mean.total<=0 ~ "<= 0",
                              Mean.total>0 & Mean.total <=1 ~ "0 < Mean <= 1.0",
                              Mean.total>1 & Mean.total <2 ~ "1.0 < Mean < 2.0",
                              Mean.total>=2 ~ "Mean >= 2.0")) %>% 
  mutate(mean.color = factor(mean.color,levels=c("<= 0","0 < Mean <= 1.0","1.0 < Mean < 2.0","Mean >= 2.0"))) %>% 
  mutate(mean.size = factor(mean.size,levels=c("<= 0","0 < Mean <= 1.0","1.0 < Mean < 2.0","Mean >= 2.0"))) %>% 
  mutate(taxa.name.new = ifelse(taxa.name.new=="d__Bacteria--p__Bacteroidota","Bacteria",taxa.name.new))

# select(taxa.name.new,Mean)
#Scenario 14
scenario14.dataA <- means2 %>%
  filter(Scenario=="Scenario 14") %>% 
  filter(Model=="BaH-ZING") %>% 
  # filter(Exposure=="X.1") %>% 
  group_by(taxa.name.new) %>%
  tidylog::summarize(Mean.total = mean(Mean)) %>% 
  mutate(mean.color = case_when(Mean.total<=0 ~ "<= 0",
                                Mean.total>0 & Mean.total <=1 ~ "0 < Mean <= 1.0",
                                Mean.total>1 & Mean.total <2 ~ "1.0 < Mean < 2.0",
                                Mean.total>=2 ~ "Mean >= 2.0")) %>% 
  mutate(mean.size= case_when(Mean.total<=0 ~ "<= 0",
                              Mean.total>0 & Mean.total <=1 ~ "0 < Mean <= 1.0",
                              Mean.total>1 & Mean.total <2 ~ "1.0 < Mean < 2.0",
                              Mean.total>=2 ~ "Mean >= 2.0")) %>% 
  mutate(mean.color = factor(mean.color,levels=c("<= 0","0 < Mean <= 1.0","1.0 < Mean < 2.0","Mean >= 2.0"))) %>% 
  mutate(mean.size = factor(mean.size,levels=c("<= 0","0 < Mean <= 1.0","1.0 < Mean < 2.0","Mean >= 2.0"))) %>% 
  mutate(taxa.name.new = ifelse(taxa.name.new=="d__Bacteria--p__Bacteroidota","Bacteria",taxa.name.new))

# select(taxa.name.new,Mean)
#Scenario 15
scenario15.dataA <- means2 %>%
  filter(Scenario=="Scenario 15") %>% 
  filter(Model=="BaH-ZING") %>% 
  group_by(taxa.name.new) %>%
  tidylog::summarize(Mean.total = mean(Mean)) %>% 
  mutate(mean.color = case_when(Mean.total<=0 ~ "<= 0",
                                Mean.total>0 & Mean.total <=1 ~ "0 < Mean <= 1.0",
                                Mean.total>1 & Mean.total <2 ~ "1.0 < Mean < 2.0",
                                Mean.total>=2 ~ "Mean >= 2.0")) %>% 
  mutate(mean.size= case_when(Mean.total<=0 ~ "<= 0",
                              Mean.total>0 & Mean.total <=1 ~ "0 < Mean <= 1.0",
                              Mean.total>1 & Mean.total <2 ~ "1.0 < Mean < 2.0",
                              Mean.total>=2 ~ "Mean >= 2.0")) %>% 
  mutate(mean.color = factor(mean.color,levels=c("<= 0","0 < Mean <= 1.0","1.0 < Mean < 2.0","Mean >= 2.0"))) %>% 
  mutate(mean.size = factor(mean.size,levels=c("<= 0","0 < Mean <= 1.0","1.0 < Mean < 2.0","Mean >= 2.0"))) %>% 
  mutate(taxa.name.new = ifelse(taxa.name.new=="d__Bacteria--p__Bacteroidota","Bacteria",taxa.name.new))

#Scenario12
scen12 <- ggtree(tree2,layout="fan",open.angle=230) %<+% scenario12.dataA +
  geom_point(aes(color=mean.color, size = mean.size), alpha=1)+
  scale_size_manual(values = c("<= 0" = 1, "0 < Mean <= 1.0" = 2, "1.0 < Mean < 2.0" = 3, "Mean >= 2.0" = 4))+
  scale_color_manual(values = c("<= 0"="#9a9a9a", "0 < Mean <= 1.0"="#e69f00", "1.0 < Mean < 2.0"="#cc0000", "Mean >= 2.0"="#660000"))+
  
  # geom_point(aes(color=Mean.total, size = Mean.total), alpha=1)+
  # scale_color_gradient2(low = "gray", mid = c("gray","yellow","orange","#CC0066"),
  #                       high = "#CC0066", midpoint = 1.2, limits = c(min, max), na.value = NA) +
  theme(legend.position = "none")

#Scenario13
scen13 <- ggtree(tree2,layout="fan",open.angle=230) %<+% scenario13.dataA +
  geom_point(aes(color=mean.color, size = mean.size), alpha=1)+
  scale_size_manual(values = c("<= 0" = 1, "0 < Mean <= 1.0" = 2, "1.0 < Mean < 2.0" = 3, "Mean >= 2.0" = 4))+
  scale_color_manual(values = c("<= 0"="#9a9a9a", "0 < Mean <= 1.0"="#e69f00", "1.0 < Mean < 2.0"="#cc0000", "Mean >= 2.0"="#660000"))+
  
  # geom_point(aes(color=Mean.total, size = Mean.total), alpha=1)+
  # scale_color_gradient2(low = "gray", mid = c("gray","yellow","orange","#CC0066"), 
  #                       high = "#CC0066", midpoint = 1.2, limits = c(min, max), na.value = NA) +
  theme(legend.position = "none")

#Scenario14
scen14 <- ggtree(tree2,layout="fan",open.angle=230) %<+% scenario14.dataA +
  geom_point(aes(color=mean.color, size = mean.size), alpha=1)+
  scale_size_manual(values = c("<= 0" = 1, "0 < Mean <= 1.0" = 2, "1.0 < Mean < 2.0" = 3, "Mean >= 2.0" = 4))+
  scale_color_manual(values = c("<= 0"="#9a9a9a", "0 < Mean <= 1.0"="#e69f00", "1.0 < Mean < 2.0"="#cc0000", "Mean >= 2.0"="#660000"))+
  
  # geom_point(aes(color=Mean.total, size = Mean.total), alpha=1)+
  # scale_color_gradient2(low = "gray", mid = c("gray","yellow","orange","#CC0066"),
  #                       high = "#CC0066", midpoint = 1.2, limits = c(min, max), na.value = NA) +
  theme(legend.position = "none")

#Scenario15
scen15 <- ggtree(tree2,layout="fan",open.angle=230) %<+% scenario15.dataA +
  geom_point(aes(color=mean.color, size = mean.size), alpha=1)+
  scale_size_manual(values = c("<= 0" = 1, "0 < Mean <= 1.0" = 2, "1.0 < Mean < 2.0" = 3, "Mean >= 2.0" = 4))+
  scale_color_manual(values = c("<= 0"="#9a9a9a", "0 < Mean <= 1.0"="#e69f00", "1.0 < Mean < 2.0"="#cc0000", "Mean >= 2.0"="#660000"))+
  
  # geom_point(aes(color=Mean.total, size = Mean.total), alpha=1)+
  # scale_color_gradient2(low = "gray", mid = c("gray","yellow","orange","#CC0066"), 
  #                       high = "#CC0066", midpoint = 1.2, limits = c(min, max), na.value = NA) +
  theme(legend.position = "none")
scen15
ggsave(plot=scen15,fs::path(dir.results,"Dendrogram_Mixture_Scenario15.jpeg"))
# 
# #9. Individual Sharing Plot ----
# #New Individual Analysis ----
# means <- formatted.data %>%
#   tidylog::filter(Exposure!="Mixture")
# # tidylog::filter(Model!="Ridge ZING")
# means$Model <- factor(means$Model, levels = c("BaH-ZING", "ZING"))
# means$Scenario <- paste0("Scenario ",means$Scenario)
# # scenario_order <- c("Scenario 1","Scenario 2","Scenario 3","Scenario 4","Scenario 5",
# #                     "Scenario 6","Scenario 7","Scenario 8","Scenario 9","Scenario 10","Scenario 11",
# #                     "Scenario 12","Scenario 13","Scenario 14",
# #                     "Scenario 15")
# # means$Scenario <- factor(means$Scenario,levels=scenario_order)
# 
# taxa_order <- c("Species","Genus","Family","Order","Class","Phylum")
# means$Taxa.Level <- factor(means$Taxa.Level, levels = taxa_order)
# 
# #Filter to just a few replications for simplicity
# # means.hold <- means
# 
# means<-means %>%
#   filter(Replication==1 | Replication ==50 |Replication==150 |Replication==200|
#            Replication==30 | Replication==60|Replication==80 |
#            Replication==250|Replication==300|Replication==350|
#            Replication==400|Replication==450|
#            Replication==500|Replication==550|
#            Replication==600)
# 
# #Species & Genus Level
# names2 <- names %>% 
#   select(species,genus) %>% 
#   rename(Taxa=species,
#          group=genus)
# 
# species <- means %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Taxa.Level,Replication,indicator) %>% 
#   filter(Taxa.Level=="Species") %>% 
#   select(-Taxa.Level)
# species <- tidylog::left_join(species,names2) 
# species <- species %>% 
#   rename(Species=Taxa) %>%
#   rename(Genus=group)
# genus <- means %>% 
#   filter(Taxa.Level=="Genus") %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Replication,indicator) %>% 
#   rename(Genus=Taxa)
# species.genus <- tidylog::left_join(species,genus,by=c("Genus","Exposure","Model","Scenario","Replication"))
# species.genus <- species.genus %>% 
#   rename(Mean.species=Mean.x,
#          Mean.genus=Mean.y,
#          Indicator.species=indicator.x,
#          Indicator.genus=indicator.y)
# 
# #Add Family 
# names3 <- names %>% 
#   select(genus,family) %>% 
#   rename(Taxa=genus,
#          group=family) %>% 
#   distinct()
# 
# genus <- means %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Taxa.Level,Replication,indicator) %>% 
#   filter(Taxa.Level=="Genus") %>% 
#   select(-Taxa.Level)
# genus <- tidylog::left_join(genus,names3) 
# genus <- genus%>% 
#   rename(Genus=Taxa) %>%
#   rename(Family=group)
# family <- means %>% 
#   filter(Taxa.Level=="Family") %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Replication,indicator) %>% 
#   rename(Family=Taxa)
# genus.family <- tidylog::left_join(genus,family,by=c("Family","Exposure","Model","Scenario","Replication"))
# genus.family <- genus.family %>% 
#   rename(Mean.genus=Mean.x,
#          Mean.family=Mean.y,
#          Indicator.genus=indicator.x,
#          Indicator.family=indicator.y)
# 
# #Add order
# names4 <- names %>% 
#   select(family,order) %>% 
#   rename(Taxa=family,
#          group=order) %>% 
#   distinct()
# 
# family <- means %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Taxa.Level,Replication,indicator) %>% 
#   filter(Taxa.Level=="Family") %>% 
#   select(-Taxa.Level)
# family <- tidylog::left_join(family,names4) 
# family <- family%>% 
#   rename(Family=Taxa) %>%
#   rename(Order=group)
# order <- means %>% 
#   filter(Taxa.Level=="Order") %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Replication,indicator) %>% 
#   rename(Order=Taxa)
# family.order <- tidylog::left_join(family,order,by=c("Order","Exposure","Model","Scenario","Replication"))
# family.order <- family.order %>% 
#   rename(Mean.family=Mean.x,
#          Mean.order=Mean.y,
#          Indicator.family=indicator.x,
#          Indicator.order=indicator.y)
# 
# #Add Class
# names5 <- names %>% 
#   select(order,class) %>% 
#   rename(Taxa=order,
#          group=class) %>% 
#   distinct()
# 
# order <- means %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Taxa.Level,Replication,indicator) %>% 
#   filter(Taxa.Level=="Order") %>% 
#   select(-Taxa.Level)
# order <- tidylog::left_join(order,names5) 
# order <- order%>% 
#   rename(Order=Taxa) %>%
#   rename(Class=group)
# class <- means %>% 
#   filter(Taxa.Level=="Class") %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Replication,indicator) %>% 
#   rename(Class=Taxa)
# order.class <- tidylog::left_join(order,class,by=c("Class","Exposure","Model","Scenario","Replication"))
# order.class <- order.class %>% 
#   rename(Mean.order=Mean.x,
#          Mean.class=Mean.y,
#          Indicator.order=indicator.x,
#          Indicator.class=indicator.y)
# 
# #Add phylum
# names6 <- names %>% 
#   select(class,phylum) %>% 
#   rename(Taxa=class,
#          group=phylum) %>% 
#   distinct()
# 
# class <- means %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Taxa.Level,Replication,indicator) %>% 
#   filter(Taxa.Level=="Class") %>% 
#   select(-Taxa.Level)
# class <- tidylog::left_join(class,names6) 
# class <- class%>% 
#   rename(Class=Taxa) %>%
#   rename(Phylum=group)
# phylum <- means %>% 
#   filter(Taxa.Level=="Phylum") %>% 
#   select(Taxa,Exposure,Model,Mean,Scenario,Replication,indicator) %>% 
#   rename(Phylum=Taxa)
# class.phylum <- tidylog::left_join(class,phylum,by=c("Phylum","Exposure","Model","Scenario","Replication"))
# class.phylum <- class.phylum %>% 
#   rename(Mean.class=Mean.x,
#          Mean.phylum=Mean.y,
#          Indicator.class=indicator.x,
#          Indicator.phylum=indicator.y)
# #Combine species.genus with genus.family
# taxa.data <- tidylog::left_join(species.genus,genus.family,by=c("Genus","Exposure","Model","Scenario","Replication","Mean.genus","Indicator.genus"))
# taxa.data <- tidylog::left_join(taxa.data,family.order,by=c("Family","Exposure","Model","Scenario","Replication","Mean.family","Indicator.family"))
# taxa.data <- tidylog::left_join(taxa.data,order.class,by=c("Order","Exposure","Model","Scenario","Replication","Mean.order","Indicator.order"))
# taxa.data <- tidylog::left_join(taxa.data,class.phylum,by=c("Class","Exposure","Model","Scenario","Replication","Mean.class","Indicator.class"))
# taxa.data2 <- taxa.data 
# # filter(Scenario=="Scenario 1") %>% 
# # filter(Model=="BaH-ZING")
# NC <- taxa.data2 %>%
#   filter(Indicator.species=="Non-Causal")
# C <- taxa.data2 %>% 
#   filter(Indicator.species=="Causal")
# C <- C %>%
#   mutate(Species="Species",
#          Genus="Genus",
#          Family="Family",
#          Order="Order",
#          Class="Class",
#          Phylum="Phylum")
# NC <- NC %>%
#   mutate(Species="Species",
#          Genus="Genus",
#          Family="Family",
#          Order="Order",
#          Class="Class",
#          Phylum="Phylum")
# max <- max(C$Mean.species,C$Mean.genus,C$Mean.family,C$Mean.order,C$Mean.class,C$Mean.phylum)
# # min <- min(C$Mean.species,C$Mean.genus,C$Mean.family,C$Mean.order,C$Mean.class,C$Mean.phylum)
# min <- min(NC$Mean.species,NC$Mean.genus,NC$Mean.family,NC$Mean.order,NC$Mean.class,NC$Mean.phylum)
# 
# 
# #add error bars
# 
# individual.plot <- ggplot() +
#   geom_segment(data = NC, aes(x = Species, xend = Genus, y = Mean.species, yend = Mean.genus, color = Mean.species), linewidth = 0.25,alpha=0.5) +
#   geom_segment(data = C, aes(x = Species, xend = Genus, y = Mean.species, yend = Mean.genus, color = Mean.species), linewidth = 0.25,alpha=0.5)+
#   
#   geom_segment(data = NC, aes(x = Genus, xend = Family, y = Mean.genus, yend = Mean.family, color = Mean.genus), linewidth = 0.25,alpha=0.5) +
#   geom_segment(data = C, aes(x = Genus, xend = Family, y = Mean.genus, yend = Mean.family, color = Mean.genus), linewidth = 0.25,alpha=0.5)+
#   
#   geom_segment(data = NC, aes(x = Family, xend = Order, y = Mean.family, yend = Mean.order, color = Mean.family), linewidth = 0.25,alpha=0.5) +
#   geom_segment(data = C, aes(x = Family, xend = Order, y = Mean.family, yend = Mean.order, color = Mean.family), linewidth = 0.25,alpha=0.5)+
#   
#   geom_segment(data = NC, aes(x = Order, xend = Class, y = Mean.order, yend = Mean.class, color = Mean.order), linewidth = 0.25,alpha=0.5) +
#   geom_segment(data = C, aes(x = Order, xend = Class, y = Mean.order, yend = Mean.class, color = Mean.order), linewidth = 0.25,alpha=0.5)+
#   
#   geom_segment(data = NC, aes(x = Class, xend = Phylum, y = Mean.class, yend = Mean.phylum, color = Mean.class), linewidth = 0.25,alpha=0.5) +
#   geom_segment(data = C, aes(x = Class, xend = Phylum, y = Mean.class, yend = Mean.phylum, color = Mean.class), linewidth = 0.25,alpha=0.5)+
#   
#   
#   geom_point(data = NC, aes(x = Species, y = Mean.species, color = Mean.species), size=0.25) +
#   geom_point(data = NC, aes(x = Genus, y = Mean.genus, color = Mean.genus), size=0.25) +
#   geom_point(data = NC, aes(x = Family, y = Mean.family, color = Mean.family), size=0.25) +
#   geom_point(data = NC, aes(x = Order, y = Mean.order, color = Mean.order), size=0.25) +
#   geom_point(data = NC, aes(x = Class, y = Mean.class, color = Mean.class), size=0.25) +
#   geom_point(data = NC, aes(x = Phylum, y = Mean.phylum, color = Mean.phylum), size=0.25) +
#   
#   geom_point(data = C, aes(x = Species, y = Mean.species, color = Mean.species), size=0.25) +  
#   geom_point(data = C, aes(x = Genus, y = Mean.genus, color = Mean.genus), size=0.25) +
#   geom_point(data = C, aes(x = Family, y = Mean.family, color = Mean.family), size=0.25) +
#   geom_point(data = C, aes(x = Order, y = Mean.order, color = Mean.order), size=0.25) +
#   geom_point(data = C, aes(x = Class, y = Mean.class, color = Mean.class), size=0.25) +
#   geom_point(data = C, aes(x = Phylum, y = Mean.phylum, color = Mean.phylum), size=0.25) +
#   geom_errorbar(data=NC,
#                 aes(Species,Mean.species,group = Species),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_errorbar(data=C,
#                 aes(Species,Mean.species,group = Species),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_point(data=NC,
#              aes(Species,Mean.species,group = Species),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_point(data=C,
#              aes(Species,Mean.species,group = Species),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_errorbar(data=NC,
#                 aes(Genus,Mean.genus,group = Genus),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_errorbar(data=C,
#                 aes(Genus,Mean.genus,group = Genus),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_point(data=NC,
#              aes(Genus,Mean.genus,group = Genus),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_point(data=C,
#              aes(Genus,Mean.genus,group = Genus),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_errorbar(data=NC,
#                 aes(Family,Mean.family,group = Family),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_errorbar(data=C,
#                 aes(Family,Mean.family,group = Family),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_point(data=NC,
#              aes(Family,Mean.family,group = Family),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_point(data=C,
#              aes(Family,Mean.family,group = Family),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_errorbar(data=NC,
#                 aes(Order,Mean.order,group = Order),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_errorbar(data=C,
#                 aes(Order,Mean.order,group = Order),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_point(data=NC,
#              aes(Order,Mean.order,group = Order),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_point(data=C,
#              aes(Order,Mean.order,group = Order),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   
#   geom_errorbar(data=NC,
#                 aes(Class,Mean.class,group = Class),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_errorbar(data=C,
#                 aes(Class,Mean.class,group = Class),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_point(data=NC,
#              aes(Class,Mean.class,group = Class),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_point(data=C,
#              aes(Class,Mean.class,group = Class),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_errorbar(data=NC,
#                 aes(Phylum,Mean.phylum,group = Phylum),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_errorbar(data=C,
#                 aes(Phylum,Mean.phylum,group = Phylum),
#                 stat = "summary",
#                 fun.data = function(x) {
#                   mean_val <- mean(x)
#                   sd_val <- sd(x)
#                   data.frame(y = mean_val, ymin = mean_val - sd_val, ymax = mean_val + sd_val)
#                 },
#                 width = 0.1,
#                 linewidth = 0.35,
#                 color = "black"
#   ) +
#   geom_point(data=NC,
#              aes(Phylum,Mean.phylum,group = Phylum),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   geom_point(data=C,
#              aes(Phylum,Mean.phylum,group = Phylum),
#              stat = "summary",
#              fun = mean,   # Use the mean function to calculate the summary statistic
#              size = 1,     # Size of the points
#              shape = 16,   # Shape of the points
#              color = "black",  # Color of the points
#              fill = "black") +
#   
#   # scale_color_gradient2(low = "gray", mid = c("yellow","orange","#CC0066"), high = "#CC0066",
#   #                       midpoint = max/2, limits = c(min, max), na.value = NA)+
#   theme_bw()+
#   facet_wrap("Model",ncol=1)+
#   # scale_color_gradientn(colors = my_colors, limits = my_limits)
#   # scale_color_gradientn(colors = c("#9a9a9a","#e69f00","#cc0000","#660000"), 
#   #                       values = c(min, 0, 0.2, 0.4, max))
#   # scale_color_gradientn(colors = c("#9a9a9a","yellow","orange","#CC0066"), 
#   #                       values = c(-0.2489804, 0, 0.2, 0.4, 0.6550573))
#   # scale_colour_gradientn(colours = c("red","yellow","green","lightblue","darkblue"),
#   #                        values = c(1.0,0.8,0.6,0.4,0.2,0)) 
#   scale_colour_gradientn(colours = c("#660000","#660000","#cc0000","#cc0000","#e69f00","#e69f00","#9a9a9a","#9a9a9a","#9a9a9a","#9a9a9a"),
#                          values = c(1.0,0.8,0.6,0.4,0.2,0)) +
#   scale_x_discrete(name = "Taxonomic Level", limits = c("Species","Genus","Family","Order","Class","Phylum"))+
#   theme(text=element_text(family="Times",size=15),
#         strip.text.x = element_text(face="bold",size=15),
#         axis.text.x = element_text(angle=45,hjust=1,face="bold",size=15),
#         axis.title.y=element_text(face="bold"))+
#   ylab("Fold Change")+
#   labs(color = "Fold Change") +
#   theme(legend.position="left")+
#   facet_grid(Scenario ~ Model)
# individual.plot
# ggsave(plot=individual.plot,fs::path(dir.results,"Sharing Plot Scenario 7.jpeg"),height=4,width=8)
