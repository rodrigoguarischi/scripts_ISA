library(ggplot2)
library(dplyr)
library(MASS)

data <- read.delim("C:/Users/Camila/OneDrive/Área de Trabalho/Análises/banco707.txt")
head(data)

dta <- data.frame(data$samplefilename,
                  data$sexo,
                  data$faixa_etaria,
                  data$raca_cor,
                  data$imc,
                  data$obesidade,
                  data$idade)
head(dta)

barplot(data.frame(table(dta$data.obesidade))[,2], 
        names.arg = c(0, 1),  col="lightpink", border = "black", 
        xlab="Obesidade (sim=1 e não=0)",
        ylab="Frequência", main="Frequência de Obesidade",
        ylim = c(0,600))

###############################

new <- read.delim("C:/Users/Camila/OneDrive/Área de Trabalho/Análises/loh707.txt")
head(new)
cps <- read.delim("C:/Users/Camila/OneDrive/Área de Trabalho/Análises/cps.txt")
head(cps)

df_b = dta %>% inner_join(new,by="data.samplefilename")
df = df_b %>% inner_join(cps,by="data.samplefilename")
head(df)

ggplot(df, aes(x=LOH, y=data.obesidade)) + 
  geom_point() + 
  ggtitle("Obesidade e Porcentagem de LOH") +
  xlab("% de LOH nos Autossômicos") + ylab("Obesidade (sim=1 e não=0)")
m2=glm(data.obesidade~LOH, family = binomial(link="logit"), data = df)
summary(m2)

#######################################################

idade <- c()
mean_idade <- mean(df$data.idade)
for (i in 1:length(df$data.idade)) {
  idade[i] <- df$data.idade[i]-mean_idade
}

m3=glm(data.obesidade ~ CP1 + CP2 + data.sexo + 
         idade + data.sexo*idade + LOH,
       family = binomial(link="logit"), data = df)
summary(m3)
confint(m3)

step.model <- stepAIC(m3, direction = "both", trace = FALSE)
summary(step.model)

#########################################################