library(readr)
library(dplyr)
library(readxl)
library(car)
library(plyr)
library(dplyr)
library(stringr)

### Opens the base file and creates a properly ordered database
fibro <- read_excel("FIBROMIALGIA.xlsx")
fibro <- fibro[c(1:6,10,7:9,11:37)]

### Labels the groups by name
fibro$grupo   <- gsub('1', 'caminhada', fibro$grupo)
fibro$grupo   <- gsub('2', 'rml', fibro$grupo)
fibro$grupo   <- gsub('3', 'controle', fibro$grupo)

### Labels menopause into two groups
fibro$menop   <- gsub('1', 'menopausa', fibro$menop)
fibro$menop   <- gsub('0', 'nao', fibro$menop)

### Transforms the continuous BMI into discrete groups and organizes columns
fibro$imccat <- cut(fibro$imc, breaks = c(-Inf, 24.9, 30, Inf),
                    labels = c("normal","sobrepeso", "obesidade"))
fibro <- subset(fibro, select = c(1:4, 38, 5:37))

# View(fibro)

### Sets the p-Value of interest that we will use to test with logical arguments
p_value = 0.05


##########


### Loops between weeks
for (n in c("0","8","16")){
  
  ### Selects the necesary columns for the analysis of each week
  positions <- c(names(fibro[1:7]), names(select(fibro, ends_with(n))))
  fibroN <- fibro %>% select(positions)
  fibroN$grupo <- ordered(fibroN$grupo, levels = c("caminhada", "rml", "controle"))
  
  ### Creates a dataframe that receives the results (names the lines)
  results <- data.frame(c("shapiro caminhada", "shapiro rml",
                          "shapiro controle", "Bartlett all","ANOVA p-value",
                          "Krustal-Wallis p-value",
                          "shapiro residuals", "Normal Residuals", "Tuckey rml-caminhada",
                          "Tuckey controle-caminhada", "Tuckey controle-rml",
                          "Chi-Squared for Menopause", "Chi-Squared for IMC",
                          "two-way ANOVA for Grupo and IMC",
                          "two-way ANOVA for Grupo and Menopausa"), row.names = TRUE)
  
  ### Creates the boxplot lattice
  namebox <- sprintf("Boxplot_das_variaveis_na_semana_%s.svg",n)
  svg(namebox, width = 16, height = 9)
  par(mfrow = c(3,5))
  
  ### Loops between variables (in columns)
  for (i in positions[-(1:2)]) {
    
    
    ### Writes results to dataframe for continuous variable
    if (i != "menop" && i != "imccat"){
      
      ### Shapiro-Wilks Test for Normality
      #   s     <- shapiro.test(fibroN[[i]])
      s1    <- shapiro.test(subset(fibroN[[i]], fibroN$grupo == "caminhada"))
      s2    <- shapiro.test(subset(fibroN[[i]], fibroN$grupo == "rml"))
      s3    <- shapiro.test(subset(fibroN[[i]], fibroN$grupo == "controle"))
      
      ### Bartlett Test for Variance
      bart  <- bartlett.test(fibroN[[i]], fibroN$grupo)
      
      ### ANOVA - One Way Test
      aovar <- aov(fibroN[[i]] ~ fibroN$grupo)
      
      ### Krustal-Wallis Non Parametric Test
      k     <- kruskal.test(fibroN[[i]] ~ fibroN$grupo)
      
      ### Shapiro-Wilks for Normality of ANOVA Residuals 
      sres <- shapiro.test(residuals(aovar))
      nres <- sres$p.value > p_value
      
      ### Tukey Test of Honest Significant Difference between groups
      tuk   <- data.frame(c(TukeyHSD(aovar)))
      tuk21 <- tuk$fibroN.grupo.p.adj[1]
      tuk31 <- tuk$fibroN.grupo.p.adj[2]
      tuk32 <- tuk$fibroN.grupo.p.adj[3]
      
      # ### Wilcoxon–Mann–Whitney Test for Pairs
      # wmw21 <- wilcox.test(fibroN[[i]] ~ fibroN$grupo, data = fibroN)
      # wmw31 <- 
      # wmw32 <- 
      
      ### ANOVA - Two Way Test
      aoimc <- aov(fibroN[[i]] ~ fibroN$grupo + fibroN$imccat)
      aomen <- aov(fibroN[[i]] ~ fibroN$grupo + fibroN$menop)
      
      results[,names(fibroN[i])] <- c(s1$p.value, s2$p.value, s3$p.value,
                                      bart$p.value, data.frame(c(summary(aovar)))$Pr..F.[1],
                                      k$p.value,
                                      sres$p.value, nres, tuk21, tuk31, tuk32, NA, NA,
                                      data.frame(c(summary(aoimc)))$Pr..F.[1],
                                      data.frame(c(summary(aomen)))$Pr..F.[1])
      
      ### Organizes Boxplots to increase readability and visualization
      if (i == "fiq0"){  
        frame()
      }
      
      if (i == "eav8" || i == "eav16"){  
        frame()
        frame()
      }
      
      ### Writes Boxplots to lattice
      boxplot(fibroN[[i]] ~ fibroN$grupo, main = names(fibroN[i]))
      
    }
    
    ### Writes results to dataframe if categorical variable
    if (i == "menop" || i == "imccat"){
      
      ### Chi-Squared Test for the discrete variables
      csmen <- chisq.test(fibro$menop, fibro$grupo)
      csimc <- chisq.test(fibro$imccat, fibro$grupo)
      
      results[,names(fibroN[i])] <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                      NA, csmen$p.value, csimc$p.value, NA, NA)
      
      
    }
    
    ### DEBUG
    #print(i)
  }
  
  ### Writes results to CSV file
  write.csv(results, sprintf("results_%s.csv", n), row.names = TRUE)
  
  dev.off()
  
  ### DEBUG
  #print(n)
}

### cleans the environment
rm(list = ls())


##########


# ### Opens reorganizes and visualizes the results for each week
# results0 <- read_csv("results_0.csv")
# transposed_res0  <- t(results_0)
# View(results_0)
# results8 <- read_csv("results_8.csv")
# transposed_res8  <- t(results_8)
# View(results_8)
# results16 <- read_csv("results_16.csv")
# transposed_res16 <- t(results_16)
# View(results_16)
# # transposed_res[7] <- as.logical(transposed_res[7])


##########


# ### Selection for Graph Creator
# n = "0" #semana desejada
# i = "estger0" #variavel desejada, verificar que variavel é da mesma semana
# g = "rml" #grupo desejado
# 
# columns <- c(names(fibro[1:2]), i)
# plotc <- fibro %>% select(columns)
# plotc$grupo <- ordered(plotc$grupo, levels = c("caminhada", "rml", "controle"))
# x <- subset(plotc[[i]], plotc$grupo == g)
# 
# ### Creates the Plot Lattice
# namehist <- sprintf("Histograma_da_variavel_%s_na_semana_%s_para_o_grupo_%s.svg", i, n, g)
# svg(namehist, width = 16, height = 9)
# par(mfrow = c(1,2))
# 
# ### Makes histograms with fitted normal curves and density lines for each variable of each week
# h <- hist(x, breaks = 13, col="white", xlab="Expressão da Variável",
#           ylab = "# of pacients", main = sprintf("variavel_%s_na_semana_%s", i, n))
# 
# ### Normal curve fit
# diff <- (max(x, na.rm = TRUE)-min(x, na.rm = TRUE))/8
# xfit <- seq((min(x, na.rm = TRUE)-diff),(max(x, na.rm = TRUE)+diff),length=100)
# yfit <- dnorm(xfit,mean=mean(x, na.rm = TRUE),sd=sd(x, na.rm = TRUE))
# yfit <- yfit*diff(h$mids[1:2])*length(x)
# lines(xfit, yfit, col="blue", lwd=2)
# 
# ### Density curve
# d <- density(x, na.rm = TRUE)
# lines(x = d$x, y = d$y * length(x) * diff(h$breaks)[1], lwd = 2, col="red")
# 
# 
# ### QQ-Plot for Normality
# qqnorm(x)
# qqline(x,col="blue")
# 
# dev.off()