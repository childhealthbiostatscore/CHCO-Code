plot_analysis = function(my_phenotype, my_id, my_outcome, sumstats, to_file = FALSE)
{
  sumstat   = sumstats %>% filter(phenotype == my_phenotype & ID == my_id)
  indata    = readRDS(sumstat$file_model)
  toplot    = indata$data
  toplot$gt = toplot[,my_id]
  toplot$y  = toplot[,my_outcome]
  my_ref    = unlist(strsplit(my_id, ":"))[[3]]
  my_alt    = unlist(strsplit(my_id, ":"))[[4]]
  my_name   = ifelse(test = sumstat$analysis [[1]] == "proteomics.standard_covariates",
                     yes  = sumstat$Target   [[1]],
                     no   = sumstat$phenotype[[1]]
                     )
  
  gt_labs   = c(paste(my_ref, my_ref, sep = "/"),
                paste(my_ref, my_alt, sep = "/"),
                paste(my_alt, my_alt, sep = "/"))

  toplot    = toplot %>%
    mutate(refalt = factor(gt, 
                           levels = 0:2, 
                           labels = gt_labs
                           ))
  
  if(to_file == TRUE)
  {
    outfile = here("output/plots", paste(my_name, my_id, my_outcome, "png", sep = "."))
    png(filename = outfile,      # file name
        width    = 11, 
        height   = 11,    
        units    = "in",
        res      = 300)                   
  }
  
  layout(rbind(1:2, 3:4))
  
  mybox = boxplot(value ~ refalt, data = toplot, 
          outline = FALSE,
          col = c("#ffdd88", "#ffcc44", "#ffa500"), 
          horizontal = FALSE,
          xlab = my_id, 
          ylab = my_name,
          main = "phenotype ~ genotype"
  )
  
  boxplot(value ~ y, data = toplot, 
          outline = FALSE,
          col = c("#dddddd", "#cc7777"), 
          horizontal = FALSE,
          xlab = my_outcome, 
          ylab = my_name,
          main = "phenotype ~ outcome"
  )
  
  mytab = as.matrix(table(toplot %>% select(gt, y)))
  mytab = mytab / rowSums(mytab)
  rownames(mytab) = gt_labs
  
  barplot(mytab[,2] * 100,
          col  = c("#ffdd88", "#ffcc44", "#ffa500"), 
          ylim = c(0,100),
          xlab = my_id, 
          ylab = my_outcome,
          main = "outcome ~ genotype"
  )
  
  plot(1,1, type = "n",
       xlim = c(0.7, 3.3),
       ylim = range(mybox$stats),
       axes = FALSE,
       xlab = my_id, 
       ylab = my_name,
       main = "phenotype ~ genotype + outcome"
  )
  
  axis(2)
  axis(1, at = 1:3, labels = gt_labs)
  
  for(my_gt in 0:2)
  {
    x0 = toplot %>% filter(gt == my_gt & y == 0)
    x1 = toplot %>% filter(gt == my_gt & y == 1)
    
    if(nrow(x0) > 0)
    {
      xpos      = my_gt + 1
      box_width = 0.3
      bx        = boxplot(x0$value, plot = FALSE)
      
      # Draw whiskers
      segments(xpos - box_width/2, bx$stats[1,1], xpos - box_width/2, bx$stats[5,1])
      
      segments(xpos - box_width/4, bx$stats[1,1], xpos - box_width*3/4, bx$stats[1,1]) # bottom cap
      segments(xpos - box_width/4, bx$stats[5,1], xpos - box_width*3/4, bx$stats[5,1]) # top cap
      
      rect(xleft = xpos - box_width, xright = xpos,
           ybottom = bx$stats[2,1], ytop = bx$stats[4,1],
           col = "#dddddd", border = "black")
      
      segments(x0 = xpos - box_width, 
               x1 = xpos , 
               y0 = bx$stats[3,1],
               y1 = bx$stats[3,1], 
               lwd = 2)
    }
    
    if(nrow(x1) > 0)
    {

      xpos      = my_gt + 1.3
      box_width = 0.3
      bx        = boxplot(x1$value, plot = FALSE)
      
      # Draw whiskers
      segments(xpos - box_width/2, bx$stats[1,1], xpos - box_width/2, bx$stats[5,1])
      
      segments(xpos - box_width/4, bx$stats[1,1], xpos - box_width*3/4, bx$stats[1,1]) # bottom cap
      segments(xpos - box_width/4, bx$stats[5,1], xpos - box_width*3/4, bx$stats[5,1]) # top cap
      
      rect(xleft = xpos - box_width, xright = xpos,
           ybottom = bx$stats[2,1], ytop = bx$stats[4,1],
           col = "#cc7777", border = "black")
      
      segments(x0 = xpos - box_width, 
               x1 = xpos , 
               y0 = bx$stats[3,1],
               y1 = bx$stats[3,1], 
               lwd = 2)
    }
    
  }
  
  if(to_file == TRUE)
  {
    dev.off()
  }
  
}
