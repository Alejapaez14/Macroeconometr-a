# Llamado de librerias
library(dplyr)
library(ggplot2)
library(patchwork)

# Importando archivo .csv
data = read_excel("Universidad/Semestre 6/Macroeconometria/inflacion.xlsx")

# Inputs 
K = 36                                  # 'Truncado' para (T / 4), si no ingrese el valor numérico
alpha = 0.05                            # Nivel |de significancia
Type = 'Cramer'                         # Identificación por 'Cramer' o 'YW' (Yule - Walker)
R = 12                                  #Total de observaciones por grupo (ya que es un ipc mensual y son 12 meses )


# Parámetros y transformaciones lineales
X = data$Inflacion
X = c(30.210, 30.321, 30.349, 30.430, 30.433, 30.538, 30.656, 30.690, 30.978, 31.301, 31.305, 31.541, 31.780, 31.777, 31.872, 31.914, 31.979, 32.173, 32.330, 32.481, 32.561, 32.568, 32.744, 33.021, 33.349, 33.487, 33.614, 33.786, 33.856, 34.011, 33.984, 34.294, 34.407, 34.441, 34.498, 34.660, 34.814, 34.923, 35.113, 35.335, 35.403, 35.666, 35.800, 36.037, 36.200, 36.226, 36.462, 36.586, 37.117, 37.425, 37.754, 38.352, 38.761, 39.076, 40.078, 40.722, 41.691, 42.224, 42.744, 44.405, 45.996, 47.033, 47.396, 48.041, 48.417, 48.896, 49.603, 50.128, 50.696, 51.702, 53.137, 53.552, 54.237, 54.537, 54.880, 55.344, 56.084, 57.036, 57.494, 57.992, 58.413, 58.713, 59.124, 59.606, 60.759, 61.894, 62.502, 62.939, 63.380, 63.633, 64.170, 64.787)
T = length(X)

Xm = matrix(data = X, nrow = 12)
N = length(X)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
if (is.wholenumber(T / R) == 'FALSE'){
  H = floor(T / R)
  if (T > (H * R)){
    l = T - (H * R)
    Xm[((l + 1) : R), (H + 1)] = NA
  }
} else if (is.wholenumber(T / R) == 'TRUE'){
  H = T / R
  Xm = matrix(data = X, nrow = 12)
}
Xm

#Estabilización de la varianza
#Creación del vector de lambdas 
#Creación del vector de lambdas 
lambdas = c(-2,-1,-0.5, 0, 0.5, 1,2)
nlambdas = length(lambdas)

#Creación de matrices
Zb = matrix(data = NA, ncol = H, nrow = 1)
Zvar = matrix(data = NA, ncol = H, nrow = 1) 
Xt =  matrix(data = NA, ncol = H, nrow = nlambdas) 
Ml = matrix(data = NA, nrow = nlambdas, ncol = 4)
colnames(Ml) = c("Lambda", "M(lambda)","S(lambda)","CV(lambda)")
Ml[,1] = lambdas

#Calculo de media, varianzas y coeficientes transformados por grupo
for(h in 1 : H){
  Zb[h] = mean(Xm[, h])                         #media
  Zvar[h] = sqrt((sum((Xm[, h] - Zb[h])** 2)) / (R - 1))  #Varianza
}  

for(l in 1 : nlambdas){
  for(h in 1 : H){
    Xt[h, l] = Zvar[h] / (Zb[h] ** (1 - lambdas[l]))
  }
}

for(l in 1 : nlambdas){
  Ml[l,2] = mean(Xt[,l])
  Ml[l,3] = sqrt((sum((Xt[,l]- Ml[l , 2])**2))/(H-1))
  Ml[l,4] = Ml [l,3]/ Ml[l,2]
}
Ml

#Estabilización de la media
Xm = 1/3* Xm[1:12,1:8]

for(i in 1:T){
  for(j in 0:3){
    S = (1/(T-j-1))* sum(((i+j+1)^j * Xm[,j+1] - sum(((i+j+2)^j) * Xm[,j+2]))/(T-j))^2
  }
}

S = sqrt(S)
S



#Metodo de Cramer
  # Condición de K
  if (K == 'Truncado'){
    K = length(data[, 1]) / 4
  } else {
    K = K
  }
M_Cramer = function(X, K, alpha, ACF_d = FALSE, PACF_d = FALSE){
  # Cálculo de variables
  Zb = sum(X) / T                                         # Media muestral
  Cv = qt(p = (alpha / 2), df = (T - 1), lower.tail = FALSE) / sqrt(T)        # Valor crítico del intervalo de confianza
  
  # Función de autocovarianza
  Gamma = matrix(data = NA)                              
  for (i in 0 : K){
    Gamma[i + 1] = sum((X[1 : (T - i)] - Zb) * (X[(i + 1) : T] - Zb)) / (T - 1)
  }
  
  # Construcción ACF
  ACF = data.frame(matrix(data = NA, nrow = (K + 1), ncol = 2))
  ACF[, 1] = seq(1, (K + 1), 1)
  colnames(ACF) = c('K', 'Rho')
  for (i in 0 : K){
    ACF[i + 1, 2] = Gamma[i + 1] / Gamma[1]
  }
  ACF = with(ACF, data.frame(K, Rho))
  
  # Construcción PACF
  PACF = data.frame(matrix(data = NA, nrow = (K + 1), ncol = 2))  
  PACF[, 1] = seq(1, (K + 1), 1)
  colnames(PACF) = c('K', 'Phi')                            
  for (k in 1 : K){
    Den = matrix(data = NA, nrow = k, ncol = k)
    for(i in 1 : k){
      for(j in 1 : (k - i + 1)){
        Den[j, j + i - 1] = ACF[i, 2]
        Den[j + i - 1, j] = ACF[i, 2]
      }
    }
    Num = Den
    Num[, k] = ACF[2 : (k + 1), 2]
    PACF[k, 2] = det(Num) / det(Den)
  }
  PACF = PACF[-(K + 1),]
  PACF = with(PACF, data.frame(K, Phi))
  
  # Gráfica conjunta de ACF y PACF
  G_ACF <- ggplot(data = ACF, mapping = aes(x = K, y = Rho)) +
    labs(title = 'ACF') +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = K, yend = 0)) +
    geom_hline(aes(yintercept = Cv), color="blue",
               linetype = "dashed") + 
    geom_hline(aes(yintercept = -Cv), color="blue",
               linetype = "dashed")
  
  G_PACF <- ggplot(data = PACF, mapping = aes(x = K, y = Phi)) +
    labs(title = 'PACF') +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = K, yend = 0)) + 
    geom_hline(aes(yintercept = Cv), color="blue",
               linetype = "dashed") + 
    geom_hline(aes(yintercept = -Cv), color="blue",
               linetype = "dashed")
  
  # Condición de la función
  if (PACF_d == TRUE){
    print(head(PACF))
  }
  if (ACF_d == TRUE){
    print(head(ACF))
  }
  
  # Grafiación de las pruebas
  G_ACF / G_PACF
}

# Método de identificación
if (Type == 'Cramer'){
  M_Cramer(X = X, K = K, alpha = alpha, ACF_d = FALSE, PACF_d = FALSE)     # En caso de necesitar ver los datos cambiar a TRUE
}