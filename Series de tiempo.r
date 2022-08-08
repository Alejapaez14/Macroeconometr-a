# Llamado de librerias
library(dplyr)
library(ggplot2)
library(lmtest)
library(gtools)
library(shiny)

# Importando archivo .csv
File = read.csv(file = ('../Ejercicio inflación/Inflación.csv'), skip = 0, header = TRUE)
head(File)
sprintf('El número de observaciones son: %s', length(File[, 1]))

# Inputs 
k = 36                                  # 'Truncado' para (T / 4), si no ingrese el valor numérico

# Parámetros y transformaciones lineales
X = File$Inflacion
T = length(X)

# Construcción de los vectores
X = cbind(X)

# Condición de k
if (k == 'Truncado'){
    k = length(File[, 1]) / 4
} else {
    k = k
}

# Cálculo de variables
Zb = sum(X) / T                                         # Media muestral
Gamma = matrix(data = NA)                               # Auto - covarianza
for (i in 0 : k){
    Gamma[i + 1] = sum((X[1 : (T - i)] - Zb) * (X[(i + 1) : T] - Zb)) / (T - 1)
}
ACF = matrix(data = NA, nrow = 36, ncol = 1)                                 # ACF
for (i in 0 : k){
    ACF[i + 1] = Gamma[i + 1] / Gamma[1]
}
PACF = matrix(data = NA)                                 # PACF
Den = matrix(data = NA, nrow = k, ncol = k)
Num = matrix(data = NA, nrow = k, ncol = k)
for(i in 1 : k){
    for(j in 1: (k - i + 1)){
        Den[j, j + i - 1] = ACF[i]
        Den[j + i - 1, j] = ACF[i]
    }
}
for (i in 2 : k + 1){
    Num = Den
    Num[1 : i, i] = ACF[2 : (i + 1)] 
    PACF[i] = det(Num[1 : i, 1 : i]) / det(Den[1 : i, 1 : i])
}


for (i in 2 : k){
    for (j in 2 : k){
        Num = Den
        Num[1 : j, j] = ACF[2 : (j + 1)] 
    }
    PACF[i] = det(Num[1 : i, 1 : i]) / det(Den[1 : i, 1 : i])
}
PACF

Num[1 : 4, 4] = ACF[2 : 5]

det(Num[1 : 4, 1 : 4])

det(Den[1 : 4, 1 : 4])


Num[1 : 2, 2] = ACF[2 : 3]

Num[1 : 2, 1 : 2]

Num[, 1 : 2]
Den[1 : 2, 1 : 2]

# Gráfica conjunta de
# Auto - covarianzas (Gamma)
plot(Gamma, type = 'h')
