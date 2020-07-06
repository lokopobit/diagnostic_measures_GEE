library(dplyr)
library(geepack)
library(readr)

nepal_dat <- read_csv("C:/Users/Marcos Lopez/Desktop/quinto/TFG/data/nepal.dat.txt")
View(nepal_dat)
working <- c('independence','exchangeable','ar1','unstructured','userdefined')
std.errs <- c('san.se',"fij",'jack', 'j1s')


nepal_dat <- na.omit(nepal_dat)
fit.gee.po <- geeglm(alive ~ sex+wt+ht+arm+bf+mage+lit+died, data = nepal_dat, 
                      family = poisson, std.err = std.errs[1],
                      corstr = working[2], id=nepal_dat$id)

betas <- c()
for (i in seq(2,length(fit.gee.po$coefficients))) {
  betas <- c(betas, fit.gee.po$coefficients[i])
}
#nepal_dat[nepal_dat$id == 120011,]$

id_unicos <- nepal_dat$id[!duplicated(nepal_dat$id)]
nepal_dat <- select(nepal_dat, -day,-month,-year)

##############################
###### CALCULO DE LA MATRIZ J1
##############################
for (num_ids in seq(1, length(id_unicos))) {
  num_id <- sum(nepal_dat$id == id_unicos[num_ids])
  matWi <- diag(num_id)
  matRi <- diag(num_id)
  for (i in 1:num_id) {
    for (j in 1:num_id) {
      if (i != j) {
        matRi[i,j] <- 0.1
      }
    }
  }
  
  matXi <- nepal_dat[nepal_dat$id == id_unicos[num_ids],2:9]
  
  aux <- c()
  for (i in seq(1,num_id)) {
    aux <- c(aux, betas%*%t(matXi[i,]))
  }
  
  matgT <- c()
  for (i in seq(1, length(betas))) {
    matgT <- rbind(matgT, aux)
  }
  
  matDiT <- t(matXi) / matgT 
  
  gamma <- 1
  matAi <- diag(num_id)
  aux1 <- data.frame(fitted = fit.gee.po$fitted.values, id = nepal_dat$id)
  for (i in seq(1, num_id)) {
    for (j in seq(1, num_id)) {
      if (i == j) {
        matAi[i,j] <- aux1[aux1$id == id_unicos[num_ids], ]$fitted[i]
      }
    }
  }
  
  matVi <- gamma*sqrt(matAi) %*% matRi %*% sqrt(matAi)
  
  matJ1aux <- matDiT%*%solve(matVi)%*%t(matDiT)
  if (num_ids == 1) {
    matJ1 <- matJ1aux
  } else {
    matJ1 <- matJ1 + matJ1aux
  }
  
}
matJ1 <- solve(matJ1)

#######################################
##### CALCULO DE LEVERAGES PARA CADA ID
#######################################
leverages <- c()
for (num_ids in seq(1, length(id_unicos))) {
  num_id <- sum(nepal_dat$id == id_unicos[num_ids])
  matWi <- diag(num_id)
  
  matXi <- nepal_dat[nepal_dat$id == id_unicos[num_ids],2:9]
  
  aux <- c()
  for (i in seq(1,num_id)) {
    aux <- c(aux, betas%*%t(matXi[i,]))
  }
  
  matHi <- sqrt(matWi) %*% as.matrix(matXi) %*% matJ1 %*% t(as.matrix(matXi)) %*% sqrt(matWi)
  leverages[num_ids] <- sum(diag(matHi))
}
plot(seq(1,length(id_unicos)), leverages)


########### COMPROBACION MATRIZ HAT 158, INCORRECTO
nepal_dat_test <- nepal_dat[nepal_dat$id != id_unicos[which.max(leverages)], ]
fit.gee.po_test <- geeglm(alive ~ sex+wt+ht+arm+bf+mage+lit+died, data = nepal_dat_test, 
                      family = poisson, std.err = std.errs[1],
                      corstr = working[2], id=nepal_dat_test$id)

summary(fit.gee.po_test)$coefficients$Std.err
summary(fit.gee.po)$coefficients$Std.err


########## HAT MATRIX GLM
fit.glm.po <- glm(alive ~ sex+wt+ht+arm+bf+mage+lit+died, data = nepal_dat, 
                      family = poisson)

aux <- data.frame(hat = hatvalues(fit.glm.po), id = nepal_dat$id)
leverages1 <- c()
for (i in seq(1,length(id_unicos))) {
  aux1 <- sum(aux[aux$id == id_unicos[i],]$hat)
  leverages1[i] <- aux1
}
plot(seq(1,length(id_unicos)), leverages1)

################################################
############## COOKS DISTANCE EXACTA, O FUNCIONA
################################################
cooksi <- c()
p <- length(betas)
gammaest <- 1
for (i in seq(1,length(id_unicos))) {
  nepal_dat2i <- nepal_dat[nepal_dat$id != id_unicos[i], ]
  fit.gee.po2i <- geeglm(alive ~ sex+wt+ht+arm+bf+mage+lit+died, data = nepal_dat2i, 
                        family = poisson, std.err = std.errs[1],
                        corstr = working[2], id=nepal_dat2i$id)
  
  betas2i <- c()
  for (j in seq(2,length(fit.gee.po$coefficients))) {
    betas2i <- c(betas2i, fit.gee.po2i$coefficients[j])
  
  }
  
  cooksaux <- t(betas-betas2i) %*% solve(matJ1) %*% (betas-betas2i) * (1 / (p * gammaest))
  cooksi[i] <- as.numeric(cooksaux)

}
plot(seq(1,length(id_unicos)), cooksi)

########### COMPROBACION COOKS 97, INCORRECTO
nepal_dat_test <- nepal_dat[nepal_dat$id != id_unicos[which.max(cooksi)], ]
fit.gee.po_test <- geeglm(alive ~ sex+wt+ht+arm+bf+mage+lit+died, data = nepal_dat_test, 
                          family = poisson, std.err = std.errs[1],
                          corstr = working[2], id=nepal_dat_test$id)

summary(fit.gee.po_test)$coefficients$Std.err
summary(fit.gee.po)$coefficients$Std.err

########## COOKS GLM
fit.glm.po <- glm(alive ~ sex+wt+ht+arm+bf+mage+lit+died, data = nepal_dat, 
                  family = poisson)

aux <- data.frame(cook = cooks.distance(fit.glm.po), id = nepal_dat$id)
cooksi1 <- c()
for (i in seq(1,length(id_unicos))) {
  aux1 <- sum(aux[aux$id == id_unicos[i],]$cook)
  cooksi1[i] <- aux1
}
plot(seq(1,length(id_unicos)), cooksi1)

############################
###################   DFFITS
############################
dffitsi <- c()
p <- length(betas)
gammaest <- 1
for (i in seq(1,length(id_unicos))) {
  num_id <- sum(nepal_dat$id == id_unicos[i])
  matWi <- diag(num_id)
  matXi <- nepal_dat[nepal_dat$id == id_unicos[i],2:9]
  aux <- c()
  for (j in seq(1,num_id)) {
    aux <- c(aux, betas%*%t(matXi[j,]))
  }
  matHi <- sqrt(matWi) %*% as.matrix(matXi) %*% matJ1 %*% t(as.matrix(matXi)) %*% sqrt(matWi)
  term2 <- (solve(matWi) - as.matrix(matXi) %*% matJ1 %*% t(as.matrix(matXi))) 
  
  term1 <- fit.gee.po$residuals
  matYmui <- matrix(ncol = ncol(term2))
  if (num_id == ncol(term2)) {
    matYmui <- rbind(matYmui, term1[1:num_id])
  } else {
    aux1 <- term1[1:num_id[jj]]
    while (length(aux1) != ncol(term2)) {
      aux1 <- c(aux1, 0)
    }
    matYmui <- rbind(matYmui, aux1)
  }
  term1 <- term1[num_id:length(term1)]
  matYmui <- matYmui[2:nrow(matYmui),]
  
  dffitsaux <- (matYmui) %*% term2 %*% matHi %*% (as.matrix(matYmui, nrow=5)) * (1 / (p * gammaest))
  dffitsi[i] <- as.numeric(dffitsaux)
  
}
dffitsi <- abs(dffitsi)
plot(seq(1,length(id_unicos)), abs(dffitsi))

########### COMPROBACION DFFITS
nepal_dat_test <- nepal_dat[nepal_dat$id != id_unicos[which.max(dffitsi)], ]
fit.gee.po_test <- geeglm(alive ~ sex+wt+ht+arm+bf+mage+lit+died, data = nepal_dat_test, 
                          family = poisson, std.err = std.errs[1],
                          corstr = working[2], id=nepal_dat_test$id)

summary(fit.gee.po_test)$coefficients$Std.err
summary(fit.gee.po)$coefficients$Std.err

##########  DFFITS GLM
fit.glm.po <- glm(alive ~ sex+wt+ht+arm+bf+mage+lit+died, data = nepal_dat, 
                  family = poisson)

aux <- data.frame(dffit = dffits(fit.glm.po), id = nepal_dat$id)
dffitsi1 <- c()
for (i in seq(1,length(id_unicos))) {
  aux1 <- sum(aux[aux$id == id_unicos[i],]$dffit)
  dffitsi1[i] <- aux1
}
plot(seq(1,length(id_unicos)), dffitsi1)


########################





 
# ############## COOKS DISTANCE APROXIMADA
# cooksi_aprox <- c()
# p <- length(betas)
# gammaest <- 1
# for (i in seq(1, length(id_unicos))) {
#   nepal_dat2i <- nepal_dat[nepal_dat$id != i, ]
#   fit.gee.po2i <- geeglm(alive ~ sex+wt+ht+arm+bf+mage+lit+died, data = nepal_dat2i, 
#                          family = poisson, std.err = std.errs[1],
#                          corstr = working[2], id=nepal_dat2i$id)
#   
#   
#   matXi <- nepal_dat[nepal_dat$id == id_unicos[num_ids],2:9] 
#   
#   term1 <-
#   term2 <- 
#   term3 <- 
#   
#   
#   
#   cooksaux <- t(betas-betas2i) %*% solve(matJ1) %*% (betas-betas2i) * (1 / (p * gammaest))
#   print(cooksaux)
#   cooksi[i] <- as.numeric(cooksaux)
#   
# }
# plot(seq(1,length(id_unicos)), cooksi_aprox)

dfftitsi <- c()
p <- length(betas)
gammaest <- 1
for (i in seq(1,length(id_unicos))) {
  num_id <- sum(nepal_dat$id == id_unicos[i])
  matWi <- diag(num_id)
  matXi <- nepal_dat[nepal_dat$id == id_unicos[i],2:9]
  aux <- c()
  for (j in seq(1,num_id)) {
    aux <- c(aux, betas%*%t(matXi[j,]))
  }
  matHi <- sqrt(matWi) %*% as.matrix(matXi) %*% matJ1 %*% t(as.matrix(matXi)) %*% sqrt(matWi)
  term2 <- (solve(matWi) - as.matrix(matXi) %*% matJ1 %*% t(as.matrix(matXi))) 
  
  term1 <- fit.gee.po$residuals
  matYmui <- matrix(ncol = ncol(term2))
  if (num_id == ncol(term2)) {
    matYmui <- rbind(matYmui, term1[1:num_id])
  } else {
    aux1 <- term1[1:num_id[jj]]
    while (length(aux1) != ncol(term2)) {
      aux1 <- c(aux1, 0)
    }
    matYmui <- rbind(matYmui, aux1)
  }
  term1 <- term1[num_id:length(term1)]
  matYmui <- matYmui[2:nrow(matYmui),]

  dffitsaux <- (matYmui) %*% term2 %*% matHi %*% (as.matrix(matYmui, nrow=5)) * (1 / (p * gammaest))
  dfftitsi[i] <- as.numeric(dffitsaux)
  
}
dfftisi <- abs(dfftitsi)
plot(seq(1,length(id_unicos)), dfftitsi)
