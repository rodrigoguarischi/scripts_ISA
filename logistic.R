library(ggplot2)
library(dplyr)
library(MASS)

data <- read.delim("C:/Users/Camila/OneDrive/Área de Trabalho/Ánalises/banco707.txt")

head(data)
dta <- data.frame(data$samplefilename,data$computed_gender,data$idade,
                  data$raca_cor,
                  data$imc,
                  data$obesidade,
                  data$CP1,data$CP2)
head(dta)

barplot(data.frame(table(dta$data.obesidade))[,2], 
        names.arg = c(0, 1),  col="lightpink", border = "black", 
        xlab="Obesidade (sim=1 e não=0)",
        ylab="Frequência", main="Frequência de Obesidade",
        ylim = c(0,600))

###############################

new <- read.delim("C:/Users/Camila/OneDrive/Área de Trabalho/Ánalises/loh707.txt")
head(new)

df= dta %>% inner_join(new,by="data.samplefilename")
head(df)

ggplot(df, aes(x=LOH, y=data.obesidade)) + 
  geom_point() + 
  ggtitle("Obesidade e Porcentagem de LOH") +
  xlab("% de LOH nos Autossômicos") + ylab("Obesidade (sim=1 e não=0)")
  #stat_smooth(method="glm", method.args=list(family="binomial"), se=FALSE)
m2=glm(data.obesidade~LOH, family = binomial(link="logit"), data = df)
summary(m2)$coef

#######################################################

m3=glm(data.obesidade~LOH + data.computed_gender + data.idade,
       family = binomial(link="logit"), data = df)
summary(m3)$coef
confint(m3)

step.model <- stepAIC(m3, direction = "both", trace = FALSE)
summary(step.model)

