library(shiny)
library(shinythemes)
library(DT)
library(dplyr)
library(ggplot2)
library(patchwork)
library(forecast)
library(tseries)

ui<-fluidPage(theme = shinytheme("darkly"),
  # # Application title
  titlePanel("Serie AR3 teorica y muestral "),
  h5("Para esta shiny app se presentaran graficas y resultados por separados con el fin de que se vea mas organizado "),
  h5("Debido a que r a veces saca error con las tildes para esta shiny app no se utilizaran"),
  
  # Sidebar with a slider input for number of bins
  tabsetPanel(tabPanel("Teorica", 
                      sidebarLayout(
                        sidebarPanel(
                          h3("Parametros para la serie AR(3) teorica"),
                          numericInput("random", "Numero de Observaciones", value = 20),
                          sliderInput("K","K", min = 0, max = 100, value =2),
                          numericInput("phi1", "Phi 1", value = 0.1, min = -1, max = 1, step = 0.0001),
                          numericInput("phi2", "Phi 2", value = 0.2, min = -1, max = 1, step = 0.0001),  
                          numericInput("phi3", "Phi 3", value = .3, min = -1, max = 1, step = 0.0001),
                          numericInput("mean", "Media", value = 2, min = 0 , max = 10000)
                        ),
                        mainPanel(
                          h3("Resultados teoricos obtenidos por la serie AR (3)"),
                          verbatimTextOutput("result"),
                        )
                      )  
                    ),
              tabPanel("Teorica Grafica",
                         mainPanel(
                           h3("Presentacion de las graficas ACF, PACF y MA infinito de un AR(3) teoricas"),
                           plotOutput("distPlot")
                         )
                       ),  
                  
             tabPanel("Muestral Resultados",
                      sidebarLayout(
                        sidebarPanel(
                          h3("Parametros para la serie AR(3) muestral"),
                          numericInput("Random", "Numero de Observaciones", value = 20),
                          sliderInput("L","k", min = 0, max = 100, value =2),
                          numericInput("Phi1", "Phi 1", value = 0.1, min = -1, max = 1, step = 0.0001),
                          numericInput("Phi2", "Phi 2", value = 0.2, min = -1, max = 1, step = 0.0001),  
                          numericInput("Phi3", "Phi 3", value = 0.3, min = -1, max = 1, step = 0.0001),
                          numericInput("Mean", "Media", value = 2, min = 0 , max = 10000),
                          numericInput("P", "Datos a pronosticar",value = 2, min = 0 , max = 10000 )
                        ),
                        mainPanel(
                          h3("Resultados muestrales obtenidos por la serie AR (3)"),
                          verbatimTextOutput("Muestra"),

                          
                          
                        )
                      )
             ),  
                    
             
             tabPanel("Muestral Graficas",
                      sidebarLayout(
                        sidebarPanel(),
                        mainPanel(
                          h3("Presentacion de las graficas media no nula, ACF y PACF de un AR(3) muestral"),
                          plotOutput("Muestral"),
                          plotOutput("Fore")
                        )
                      )
             ),
             
             tabPanel("Muestral Graficas",
                      sidebarLayout(
                        sidebarPanel(),
                        mainPanel(
                          h3("Presentacion de las graficas de At, Histograma de At y el pronostico con la cantidad de datos simulados de un AR(3) muestral"),
                          plotOutput("muestra2")
                        )
                      )
             )
        )  
)

server <- function(input,output,session){
  
    observeEvent(input$random,{
      K = input$random/4
      updateSliderInput(session, "K", max = K)
    })
    
    teor = reactive({
      
        #Inputs que se usarC!n 
        Phi = c(input$phi1, input$phi2, input$phi3)               # Phi's poblacionales para un modelo AR[3]
        alpha = 0.05                                             # Nivel de significancia       
        K = input$K                                              # TamaC1o de la muestra 
        Mu = input$mean                                  # Media poblacional diferente de 0
        N = input$random
        # SimulaciC3n de ruido blanco
        set.seed(0214)                          # Semilla para la generaciC3n de datos aleatorios
        at = rnorm(N, mean = 0, sd = 1)
        sigma_a = 1
        # CC!lculo de los factores con su parte imaginaria
        Phi = matrix(data = c(-1, c(Phi)), nrow = 4, ncol = 1) * -1
        G = matrix(data = c(as.complex(polyroot(Phi)) ** -1), nrow = 3, ncol = 1)               # Polyroot devuelve las raices, ** -1 = Factores
        for (i in 1 : length(G)){
          if (sum(round(Im(G), 10)) == 0){
            G = G
          } else {
            Re_G = Re(G)
            Im_G = Im(G)
            for (i in 1 : length(G)){                                                               # Convierten los nC:meros complejos en reales
              G[i] = sqrt(Re_G[i] ** 2 + Im_G[i] ** 2)
            }
          }
        }

        # IdentificaciC3n estacionareidad
        NG = length(G)
        G = round(Re(G), 10)
        Check = matrix(data = NA, nrow = NG, ncol = 1)
        for (i in 1 : NG){
          if (abs(G[i]) < 1){
            Check[i] = 1
          } else {
            Check[i] = 0
          }
        }

        # IdentificaciC3n de factores iguales
        Check_fact = matrix(data = sort(G, decreasing = TRUE), nrow = NG, ncol = 2)
        for (i in 2 : NG){
          Check_fact[1, 2] = 1
          if (Check_fact[i, 1] == Check_fact[(i - 1), 1]){
            Check_fact[i, 2] = Check_fact[(i - 1), 2] + 1
          } else {
            Check_fact[i, 2] = 1
          }
        }

        # CC!lculo de los valores iniciales
        if (sum(Check) == 3){

          # CC!lculo de los Rho
          Phi = Phi * -1
          Rho1 = (Phi[2] + (Phi[3] * Phi[4])) / (1 - Phi[3] - (Phi[4] * Phi[2]) - Phi[4] ** 2)
          Rho2 = Rho1 * (Phi[2] + Phi[4]) + Phi[3]
          Rho = matrix(data = c(1, Rho1, Rho2), nrow = NG, ncol = 1)

          # CC!lculo de la ACF con y sin factores distintos
          if (sum(Check_fact[, 2]) == NG){
            Ginv = rbind(c(1, 1, 1), c(G), c(G ** 2))               # CC!lculo de los valores iniciales con factores diferentes
            A = as.matrix(solve(Ginv) %*% Rho)
            print('Todos los factores son diferentes.')

            # ConstrucciC3n ACF (SoluciC3n de la funciC3n generadora)
            ACF = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
            ACF[, 1] = seq(0, K, 1)
            colnames(ACF) = c('N', 'Rho')
            for (i in 0 : N){
              ACF[i + 1, 2] = as.numeric((A[1] * (G[1] ** i)) + (A[2] * (G[2] ** i)) + (A[3] * (G[3] ** i)))
            }
            ACF = with(ACF, data.frame(N, Rho))

          }
          if (sum(Check_fact[, 2]) == (NG + 1)){
            Ginv = rbind(c(1, 1, 0), c(G), c(G[1] ** 2, G[2] ** 2, 2 * G[3] ** 2))               # CC!lculo de los valores iniciales con factores iguales
            A = as.matrix(solve(Ginv) %*% Rho)
            print('Dos factores son iguales.')

            # ConstrucciC3n ACF (SoluciC3n de la funciC3n generadora)
            ACF = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
            ACF[, 1] = seq(0, N, 1)
            colnames(ACF) = c('N', 'Rho')
            FG = matrix(data = NA, ncol = NG, nrow = (N + 1))
            for (j in 1 : NG){
              for (i in 0 : N){
                if (Check_fact[j, 2] == 1){
                  FG[(i + 1), j] = A[j] * (G[j] ** i)
                } else {
                  FG[(i + 1), j] = A[j] * (G[j] ** i) * i
                }
              }
            }

            for (i in 0 : N){
              ACF[i + 1, 2] = as.numeric(sum(FG[(i + 1), ]))
            }
            ACF = with(ACF, data.frame(N, Rho))

          }
          if (sum(Check_fact[, 2]) == (NG + 3)){
            G_3 = c(G[1], (2 * G[2]), (4 * G[3]))
            Ginv = rbind(c(1, 0, 0), c(G), c(G_3 ** 2))              # CC!lculo de los valores iniciales con factores iguales
            A = as.matrix(solve(Ginv) %*% Rho)
            print('Tres factores son iguales.')

            # ConstrucciC3n ACF (SoluciC3n de la funciC3n generadora)
            ACF = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
            ACF[, 1] = seq(0, K, 1)
            colnames(ACF) = c('N', 'Rho')
            FG = matrix(data = NA, ncol = NG, nrow = (N + 1))
            for (j in 1 : NG){
              FG[(i + 1), j] = A[j] * (G[j] ** i) * i
            }
            for (i in 0 : N){
              ACF[i + 1, 2] = as.numeric(sum(FG[(i + 1), ]))
            }
            ACF = with(ACF, data.frame(N, Rho))

          }

          # ConstrucciC3n PACF (Usando Cramer)
          PACF = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
          PACF[, 1] = seq(1, (N + 1), 1)
          colnames(PACF) = c('N', 'Phi')
          for (k in 1 : N){
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
          PACF = with(PACF, data.frame(N, Phi))
          
          # Varianza del proceso (Gamma_0_1_2)
          F1 = c((1 - (Phi[2] ** 2) - (Phi[3] ** 2) - (Phi[4] * 2)), (-2 * Phi[3] * (Phi[2] + Phi[4])), (-2 * Phi[2] * Phi[4]))
          F2 = c((-Phi[2]), (1 - Phi[3]), (-Phi[4]))
          F3 = c((-Phi[3]), -(Phi[2] + Phi[4]), (1))
          C2 = matrix(data = c(sigma_a, 0, 0), ncol = 1, nrow = 3)
          Gamma = matrix(data = rbind(F1, F2, F3), ncol = 3, nrow = 3)
          Gamma_1_2_3 = matrix(data = NA, ncol = 2, nrow = 3)
          colnames(Gamma_1_2_3) = c('Gamma', 'Valor')
          Gamma_1_2_3[, 1] = seq(0, 2, 1)
          Gamma_1_2_3[, 2] = solve(Gamma) %*% C2
          
          

          # Grafica conjunta de ACF y PACF
          Cv = qt(p = (alpha / 2), df = (N - 1), lower.tail = FALSE) / sqrt(N)        # Valor cr??tico del intervalo de confianza
          G_ACF_3 <- ggplot(data = ACF, mapping = aes(x = N, y = Rho)) +
            labs(title = 'ACF AR[3] Teorica') +
            geom_hline(aes(yintercept = 0)) +
            geom_segment(mapping = aes(xend = N, yend = 0)) +
            geom_hline(aes(yintercept = Cv), color="blue",
                       linetype = "dashed") + 
            geom_hline(aes(yintercept = -Cv), color="blue",
                       linetype = "dashed")
          
          G_PACF_3 <- ggplot(data = PACF, mapping = aes(x = N, y = Phi)) +
            labs(title = 'PACF AR[3] Teorica') +
            geom_hline(aes(yintercept = 0)) +
            geom_segment(mapping = aes(xend = N, yend = 0)) + 
            geom_hline(aes(yintercept = Cv), color="blue",
                       linetype = "dashed") + 
            geom_hline(aes(yintercept = -Cv), color="blue",
                       linetype = "dashed")   
          
          
          Psi_2 = PACF[1, 2]
          Psi_3 = (PACF[1, 2] ** 2) + PACF[2, 2]
          Psi = matrix(data = c(1, Psi_2, Psi_3), nrow = 3, ncol = 1)
          
          # CC!lculo de la ACF con y sin factores distintos
          if (sum(Check_fact[, 2]) == NG){
            Ginv = rbind(c(1, 1, 1), c(G), c(G ** 2))               # CC!lculo de los valores iniciales con factores diferentes
            A = as.matrix(solve(Ginv) %*% Psi)
            
            # ConstrucciC3n ACF (SoluciC3n de la funciC3n generadora)
            ACF_MA = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
            ACF_MA[, 1] = seq(0, N, 1)
            colnames(ACF_MA) = c('N', 'Psi')
            for (i in 0 : N){
              ACF_MA[i + 1, 2] = as.numeric((A[1] * (G[1] ** i)) + (A[2] * (G[2] ** i)) + (A[3] * (G[3] ** i)))
            }
            ACF_MA = with(ACF_MA, data.frame(N, Psi))
            
          } else {
            Ginv = rbind(c(1, 1, 0), c(G), c(G ** 2))               # CC!lculo de los valores iniciales con factores iguales
            A = as.matrix(solve(Ginv) %*% Psi)
            
            # ConstrucciC3n ACF (SoluciC3n de la funciC3n generadora)
            ACF_MA = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
            ACF_MA[, 1] = seq(0, N, 1)
            colnames(ACF_MA) = c('N', 'Psi')
            FG = matrix(data = NA, ncol = NG, nrow = (N + 1))
            for (j in 1 : NG){
              for (i in 0 : N){
                if (Check_fact[j, 2] == 1){
                  FG[(i + 1), j] = A[j] * (G[j] ** i)
                } else {
                  FG[(i + 1), j] = A[j] * (G[j] ** i) * i
                }
              }
            }
            for (i in 0 : N){
              ACF_MA[i + 1, 2] = as.numeric(sum(FG[(i + 1), ]))
            }
            ACF_MA = with(ACF_MA, data.frame(N, Psi))
            
          }
          
          # GrC!fica conjunta de ACF y PACF
          G_MA <- ggplot(data = ACF_MA, mapping = aes(x = N, y = Psi)) +
            labs(title = 'MA[inf] AR[3] Te??rico') +
            geom_hline(aes(yintercept = 0)) +
            geom_segment(mapping = aes(xend = N, yend = 0)) +
            geom_hline(aes(yintercept = Cv), color="blue",
                       linetype = "dashed") + 
            geom_hline(aes(yintercept = -Cv), color="blue",
                       linetype = "dashed")
          
          
          # Salidas
          print("Funcion Generadora de un AR 3 teorico")
          print(data.frame(G))
          print(data.frame(A))
          print("")
          print("Autocovarianzas de un AR 3 teorico")
          print(data.frame(Gamma_1_2_3))
          print("")
          print("Autocorrelaciin de un AR 3 teorico")
          print(data.frame(ACF))
          print("")
          print("Autocorrelacion parcial de un AR 3 teorico")
          print(data.frame(PACF))
          plot(G_ACF_3)
          plot(G_PACF_3)
          print("")
          print("Representacion MA infinito de un AR 3 teorico")
          print(data.frame(ACF_MA))
          plot(G_MA)
          print("")
          
          return(G_ACF_3 / G_PACF_3/G_MA)
        } else {
          print(sprintf('El proceso no es estacionario al tener sus factores |%s| > 1', round(G, 4)))
        }
        
       

        })
        
    
    output$result = renderPrint(
      teor()
    )
    
    output$distPlot = renderPlot(
      teor()
    )
    
    observeEvent(input$Random,{
      L = input$Random/4
      updateSliderInput(session, "L", max = L)
    })
    
    muest =  reactive({
      
      #Inputs que se usarC!n 
      Phi = c(input$Phi1, input$Phi2, input$Phi3)               # Phi's poblacionales para un modelo AR[3]
      alpha = 0.05                                             # Nivel de significancia       
      K = input$L                                              # TamaC1o de la muestra 
      Mu = input$Mean                                  # Media poblacional diferente de 0
      N = input$Random
      
      # Datos a pronosticar
      P = input$P
      
      # SimulaciC3n de ruido blanco
      set.seed(0214)                          # Semilla para la generaciC3n de datos aleatorios
      at = rnorm(N, mean = 0, sd = 1)
      sigma_a = 1
      # CC!lculo de los factores con su parte imaginaria
      Phi = matrix(data = c(-1, c(Phi)), nrow = 4, ncol = 1) * -1
      G = matrix(data = c(as.complex(polyroot(Phi)) ** -1), nrow = 3, ncol = 1)               # Polyroot devuelve las raices, ** -1 = Factores
      for (i in 1 : length(G)){
        if (sum(round(Im(G), 10)) == 0){
          G = G
        } else {
          Re_G = Re(G)
          Im_G = Im(G)
          for (i in 1 : length(G)){                                                               # Convierten los nC:meros complejos en reales
            G[i] = sqrt(Re_G[i] ** 2 + Im_G[i] ** 2)
          }
        }
      }
      
      # IdentificaciC3n estacionareidad
      NG = length(G)
      G = round(Re(G), 10)
      Check = matrix(data = NA, nrow = NG, ncol = 1)
      for (i in 1 : NG){
        if (abs(G[i]) < 1){
          Check[i] = 1
        } else {
          Check[i] = 0
        }
      }
      
      # IdentificaciC3n de factores iguales
      Check_fact = matrix(data = sort(G, decreasing = TRUE), nrow = NG, ncol = 2)
      for (i in 2 : NG){
        Check_fact[1, 2] = 1
        if (Check_fact[i, 1] == Check_fact[(i - 1), 1]){
          Check_fact[i, 2] = Check_fact[(i - 1), 2] + 1
        } else {
          Check_fact[i, 2] = 1
        }
      }
      
      # CC!lculo de los valores iniciales
      if (sum(Check) == 3){
        
        # CC!lculo de los Rho
        Phi = Phi * -1
        Rho1 = (Phi[2] + (Phi[3] * Phi[4])) / (1 - Phi[3] - (Phi[4] * Phi[2]) - Phi[4] ** 2)
        Rho2 = Rho1 * (Phi[2] + Phi[4]) + Phi[3]
        Rho = matrix(data = c(1, Rho1, Rho2), nrow = NG, ncol = 1)
        
        # CC!lculo de la ACF con y sin factores distintos
        if (sum(Check_fact[, 2]) == NG){
          Ginv = rbind(c(1, 1, 1), c(G), c(G ** 2))               # CC!lculo de los valores iniciales con factores diferentes
          A = as.matrix(solve(Ginv) %*% Rho)
          print('Todos los factores son diferentes.')
          
          # ConstrucciC3n ACF (SoluciC3n de la funciC3n generadora)
          ACF = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
          ACF[, 1] = seq(0, K, 1)
          colnames(ACF) = c('N', 'Rho')
          for (i in 0 : N){
            ACF[i + 1, 2] = as.numeric((A[1] * (G[1] ** i)) + (A[2] * (G[2] ** i)) + (A[3] * (G[3] ** i)))
          }
          ACF = with(ACF, data.frame(N, Rho))
          
        }
        if (sum(Check_fact[, 2]) == (NG + 1)){
          Ginv = rbind(c(1, 1, 0), c(G), c(G[1] ** 2, G[2] ** 2, 2 * G[3] ** 2))               # CC!lculo de los valores iniciales con factores iguales
          A = as.matrix(solve(Ginv) %*% Rho)
          print('Dos factores son iguales.')
          
          # ConstrucciC3n ACF (SoluciC3n de la funciC3n generadora)
          ACF = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
          ACF[, 1] = seq(0, N, 1)
          colnames(ACF) = c('N', 'Rho')
          FG = matrix(data = NA, ncol = NG, nrow = (N + 1))
          for (j in 1 : NG){
            for (i in 0 : N){
              if (Check_fact[j, 2] == 1){
                FG[(i + 1), j] = A[j] * (G[j] ** i)
              } else {
                FG[(i + 1), j] = A[j] * (G[j] ** i) * i
              }
            }
          }
          
          for (i in 0 : N){
            ACF[i + 1, 2] = as.numeric(sum(FG[(i + 1), ]))
          }
          ACF = with(ACF, data.frame(N, Rho))
          
        }
        if (sum(Check_fact[, 2]) == (NG + 3)){
          G_3 = c(G[1], (2 * G[2]), (4 * G[3]))
          Ginv = rbind(c(1, 0, 0), c(G), c(G_3 ** 2))              # CC!lculo de los valores iniciales con factores iguales
          A = as.matrix(solve(Ginv) %*% Rho)
          print('Tres factores son iguales.')
          
          # ConstrucciC3n ACF (SoluciC3n de la funciC3n generadora)
          ACF = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
          ACF[, 1] = seq(0, K, 1)
          colnames(ACF) = c('N', 'Rho')
          FG = matrix(data = NA, ncol = NG, nrow = (N + 1))
          for (j in 1 : NG){
            FG[(i + 1), j] = A[j] * (G[j] ** i) * i
          }
          for (i in 0 : N){
            ACF[i + 1, 2] = as.numeric(sum(FG[(i + 1), ]))
          }
          ACF = with(ACF, data.frame(N, Rho))
          
        }
        
        # ConstrucciC3n PACF (Usando Cramer)
        PACF = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
        PACF[, 1] = seq(1, (N + 1), 1)
        colnames(PACF) = c('N', 'Phi')
        for (k in 1 : N){
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
        PACF = with(PACF, data.frame(N, Phi))
        
        # Varianza del proceso (Gamma_0_1_2)
        F1 = c((1 - (Phi[2] ** 2) - (Phi[3] ** 2) - (Phi[4] * 2)), (-2 * Phi[3] * (Phi[2] + Phi[4])), (-2 * Phi[2] * Phi[4]))
        F2 = c((-Phi[2]), (1 - Phi[3]), (-Phi[4]))
        F3 = c((-Phi[3]), -(Phi[2] + Phi[4]), (1))
        C2 = matrix(data = c(sigma_a, 0, 0), ncol = 1, nrow = 3)
        Gamma = matrix(data = rbind(F1, F2, F3), ncol = 3, nrow = 3)
        Gamma_1_2_3 = matrix(data = NA, ncol = 2, nrow = 3)
        colnames(Gamma_1_2_3) = c('Gamma', 'Valor')
        Gamma_1_2_3[, 1] = seq(0, 2, 1)
        Gamma_1_2_3[, 2] = solve(Gamma) %*% C2
        
        
        
        # GrC!fica conjunta de ACF y PACF
        Cv = qt(p = (alpha / 2), df = (N - 1), lower.tail = FALSE) / sqrt(N)        # Valor crC-tico del intervalo de confianza
        G_ACF_3 <- ggplot(data = ACF, mapping = aes(x = N, y = Rho)) +
          labs(title = 'ACF AR[3] Teorica') +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = N, yend = 0)) +
          geom_hline(aes(yintercept = Cv), color="blue",
                     linetype = "dashed") +
          geom_hline(aes(yintercept = -Cv), color="blue",
                     linetype = "dashed")
        
        G_PACF_3 <- ggplot(data = PACF, mapping = aes(x = N, y = Phi)) +
          labs(title = 'PACF AR[3] Teorica') +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = N, yend = 0)) +
          geom_hline(aes(yintercept = Cv), color="blue",
                     linetype = "dashed") +
          geom_hline(aes(yintercept = -Cv), color="blue",
                     linetype = "dashed")
        
        Psi_2 = PACF[1, 2]
        Psi_3 = (PACF[1, 2] ** 2) + PACF[2, 2]
        Psi = matrix(data = c(1, Psi_2, Psi_3), nrow = 3, ncol = 1)
        
        # CC!lculo de la ACF con y sin factores distintos
        if (sum(Check_fact[, 2]) == NG){
          Ginv = rbind(c(1, 1, 1), c(G), c(G ** 2))               # CC!lculo de los valores iniciales con factores diferentes
          A = as.matrix(solve(Ginv) %*% Psi)
          
          # ConstrucciC3n ACF (SoluciC3n de la funciC3n generadora)
          ACF_MA = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
          ACF_MA[, 1] = seq(0, N, 1)
          colnames(ACF_MA) = c('N', 'Psi')
          for (i in 0 : N){
            ACF_MA[i + 1, 2] = as.numeric((A[1] * (G[1] ** i)) + (A[2] * (G[2] ** i)) + (A[3] * (G[3] ** i)))
          }
          ACF_MA = with(ACF_MA, data.frame(N, Psi))
          
        } else {
          Ginv = rbind(c(1, 1, 0), c(G), c(G ** 2))               # CC!lculo de los valores iniciales con factores iguales
          A = as.matrix(solve(Ginv) %*% Psi)
          
          # ConstrucciC3n ACF (SoluciC3n de la funciC3n generadora)
          ACF_MA = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
          ACF_MA[, 1] = seq(0, N, 1)
          colnames(ACF_MA) = c('N', 'Psi')
          FG = matrix(data = NA, ncol = NG, nrow = (N + 1))
          for (j in 1 : NG){
            for (i in 0 : N){
              if (Check_fact[j, 2] == 1){
                FG[(i + 1), j] = A[j] * (G[j] ** i)
              } else {
                FG[(i + 1), j] = A[j] * (G[j] ** i) * i
              }
            }
          }
          for (i in 0 : N){
            ACF_MA[i + 1, 2] = as.numeric(sum(FG[(i + 1), ]))
          }
          ACF_MA = with(ACF_MA, data.frame(N, Psi))
          
        }
        
        # GrC!fica conjunta de ACF y PACF
        G_MA <- ggplot(data = ACF_MA, mapping = aes(x = N, y = Psi)) +
          labs(title = 'MA[inf] AR[3] TeC3rico') +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = N, yend = 0)) +
          geom_hline(aes(yintercept = Cv), color="blue",
                     linetype = "dashed") + 
          geom_hline(aes(yintercept = -Cv), color="blue",
                     linetype = "dashed")
        
        ## Muestral
        # Matriz de datos de Zt 
        Zt = data.frame(matrix(data = NA, ncol = 2, nrow = N))
        colnames(Zt) = c('N', 'Zt')
        Zt[, 1] = seq(1, N, 1)
        Zt[1 : 3, 2] = A
        for (i in (length(A) + 1) : N){
          Zt[i, 2] = (PACF[1, 2] * Zt[i - 1, 2]) + (PACF[2, 2] * Zt[i - 2, 2]) + (PACF[3, 2] * Zt[i - 3, 2]) + at[i]
        }
        
        # Matriz de datos de Zt con media no nula
        Ztmu = data.frame(matrix(data = NA, ncol = 2, nrow = N))
        colnames(Ztmu) = c('N', 'ZtMu')
        Ztmu[, 1] = seq(1, N, 1)
        Ztmu[, 2] = Zt[, 2] + Mu
        Media = (sum(Ztmu[, 2]) / N)
        
        # GrC!fica de la serie
        G_Zt <- ggplot(data = Ztmu, mapping = aes(x = N, y = ZtMu, group = 12)) +
          geom_line(colour = 'red', size = 1) +
          labs(title = 'Datos simulados de un AR[3] con media no nula', x = 'n', y = 'Valor de Zt') +
          geom_hline(aes(yintercept = Media))
        
        # RedefiniciC3n de los datos
        X = Ztmu[, 2]
        
      
          # ConstrucciC3n de los vectores y parC!metros
          T = length(X)
          
          # CC!lculo de variables
          Zb = sum(X) / T                                                             # Media muestral
          Cv = qt(p = (alpha / 2), df = (T - 1), lower.tail = FALSE) / sqrt(T)        # Valor crC-tico del intervalo de confianza
        
          # FunciC3n de autocovarianza
          Gamma = matrix(data = NA)                              
          for (i in 0 : N){
            Gamma[i + 1] = sum((X[1 : (N - i)] - Zb) * (X[(i + 1) : N] - Zb)) / (N - 1)
          }
          
          # ConstrucciC3n ACF
          ACF = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
          ACF[, 1] = seq(0, N, 1)
          colnames(ACF) = c('N', 'Rho')
          for (i in 0 : N){
            ACF[i + 1, 2] = Gamma[i + 1] / Gamma[1]
          }
          ACF = na.omit(ACF)
          ACF = with(ACF, data.frame(N, Rho))
          
          # ConstrucciC3n PACF
          PACF = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))  
          PACF[, 1] = seq(1, (N + 1), 1)
          colnames(PACF) = c('N', 'Phi')                            
          for (k in 1 : N){
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
          PACF = PACF[-(N + 1),]
          PACF = na.omit(PACF)
          PACF = with(PACF, data.frame(N, Phi))
          
          # GrC!fica conjunta de ACF y PACF
          G_ACF <- ggplot(data = ACF, mapping = aes(x = N, y = Rho)) +
            labs(title = 'ACF') +
            geom_hline(aes(yintercept = 0)) +
            geom_segment(mapping = aes(xend = N, yend = 0)) +
            geom_hline(aes(yintercept = Cv), color="blue",
                       linetype = "dashed", size = 0.5) + 
            geom_hline(aes(yintercept = -Cv), color="blue",
                       linetype = "dashed", size = 0.5)
          
          G_PACF <- ggplot(data = PACF, mapping = aes(x = N, y = Phi)) +
            labs(title = 'PACF') +
            geom_hline(aes(yintercept = 0)) +
            geom_segment(mapping = aes(xend = N, yend = 0)) + 
            geom_hline(aes(yintercept = Cv), color="blue",
                       linetype = "dashed", size = 0.5) + 
            geom_hline(aes(yintercept = -Cv), color="blue",
                       linetype = "dashed", size = 0.5)
          
          
          # IdentificaciC3n del modelo AR
          CheckAR = matrix(data = NA, ncol = 1, nrow = N)
          for(i in 1 : (N - 1)) {
            CheckAR[i] = between(0, (PACF[i, 2] - Cv), (PACF[i, 2] + Cv))
            if(CheckAR[i] == 1){
              next
            }
            AR_P = c(i)
          }
          
          # IdentificaciC3n del modelo MA
          CheckMA = matrix(data = NA, ncol = 1, nrow = N)
          for(i in 1 : (N - 1)){
            CheckMA[i] = between(0, (ACF[i, 2] - Cv * (sqrt(1 + 2 * sum(ACF[(i + 1), 2] ** 2)))), (ACF[i, 2] + Cv) * (sqrt(1 + 2 * sum(ACF[(i + 1), 2] ** 2))))
            if(CheckMA[i] == 1){
              next
            }
            MA_Q = c(i)
          }
          # Inputs
          AR_P = 2
          MA_Q = 1
          J_opt = 0                       # Es 0 porque la serie es estable en nivel desde la simulaciC3n
          
          
          Model = arima(X, order = c(AR_P, J_opt, MA_Q), method = c('ML'))
          At = Model$residuals
          
          # Media 0 de residuales
          t = J_opt + AR_P + 1
          T = length(At)
          Ma = sum(At[t : T]) / (T - J_opt - AR_P)
          Sigma_at = sqrt(t(At - Ma) %*% (At - Ma) / (T - J_opt - AR_P - MA_Q))
          Stat = abs((sqrt(T - J_opt - AR_P) * Ma) / (Sigma_at)) 
          Vc = qt(p = alpha, df = (T - AR_P - MA_Q - J_opt), lower.tail = FALSE)
          
          #Varianza constante de residuales
          At1 = data.frame(matrix(data = NA, nrow = T, ncol = 2))  
          colnames(At1) = c('T', 'D.At')
          At1[, 1] = seq(1, T, 1)
          At1[, 2] = At
          At1 = with(At1, data.frame(T, D.At))
          At_sigma = 1 * Sigma_at
          
          #GrC!fica
          G_At1 <- ggplot(data = At1, mapping = aes(x = T, y = D.At)) +
            geom_line(colour = 'red', size = 1) +
            labs(title = 'Inspecci??n de residuales por varianza [At]', y = 'At') +
            geom_hline(aes(yintercept = At_sigma), color="blue",
                       linetype = "dashed", size = 1) +
            geom_hline(aes(yintercept = -At_sigma), color="blue",
                       linetype = "dashed", size = 1)
          
          #Independencia de residuales
          Rho_at = matrix(data = NA, nrow = (T - N), ncol = 1)
          for (i in 1 : (T - N)){
            Rho_at[i] = sum(At[i] * At[(i + t)]) / sum(At ** 2)
          }
          
          # Prueba de Ljung - Box
          Q = (T - J_opt - AR_P) * (T - J_opt - AR_P + 2) * sum((sum(Rho_at) ** 2) / (T - J_opt - AR_P - K)) 
          Vc = qchisq(p = alpha, df = (K - AR_P - MA_Q))
          
          # Normalidad de residuales 
            #Prueba Jarque Bera (Normalidad en el error)
          S = (sum((At - mean(At)) ** 3) / T) / (sqrt(as.numeric((t(At) %*% At) / T)) ** 3)           # AsimetrC-a
          K1 = (sum((At - mean(At)) ** 4) / T) / (sqrt(as.numeric((t(At) %*% At) / T)) ** 4)          # Kurtosis
          JB = (1 / T) * ((S ** 2 / 6) + ((K1 - 3) ** 2 / 24))
          VcJB = qchisq(alpha, df = 2, lower.tail = FALSE)
          
          At1 = data.frame(matrix(data = NA, nrow = T, ncol = 2))
          colnames(At1) = c('N', 'At')
          At1[, 1] = seq(1, T, 1)
          At1[, 2] = At
          A_h = ggplot(At1, aes(x = At)) + 
            labs(title = 'Histograma de At') +
            geom_histogram(aes(y = ..density..), colour = 'darkblue', fill = 'lightblue', alpha = 0.5, binwidth = 0.1) +
            geom_density(alpha = 0.2, fill = "#FF6666") +
            geom_vline(aes(xintercept= mean(At)), color = "blue", linetype = "dashed", size = 0.7) 
          
          #No existencia de observaciones aberrantes 
          Obs.aberr <- function(X){    
            T = length(X)
            Xmean = sum(X) / T
            Xsd = sd(X)  
            LS = Xmean + (3 * Xsd)
            LI = Xmean - (3 * Xsd)  
            D.aberr = c()
            
            for (i in 1 : T){
              if (X[i] < LS & X[i] > LI){
                D.aberr[i] = NA
              } else {
                D.aberr[i] = X[i]
              }
            }
            
            D.aberr = na.omit(D.aberr)  
            
            # Salidas
            print(sprintf('Hay un total de %s observaciones aberrantes', length(D.aberr)))
          }
          
          # GrC!fica de los datos
          At1 = data.frame(matrix(data = NA, nrow = T, ncol = 2))  
          colnames(At1) = c('T', 'D.At')
          At1[, 1] = seq(1, T, 1)
          At1[, 2] = At
          At1 = with(At1, data.frame(T, D.At))
          At_sigma = 3 * Sigma_at
          
          G_At1 <- ggplot(data = At1, mapping = aes(x = T, y = D.At)) +
            geom_line(colour = '#ff0202', size = 1) +
            labs(title = 'At') +
            geom_hline(aes(yintercept = At_sigma), color="blue",
                       linetype = "dashed", size = 1) +
            geom_hline(aes(yintercept = -At_sigma), color="blue",
                       linetype = "dashed", size = 1)
          
          #Parsimonio
          Parsimonio = function(Delta, Delta.var){
            Desv = sqrt(Delta.var)
            Vec = Delta
            for (i in 1 : length(Delta)){
              LS = Delta[i] + (2 * Desv[i])
              LI = Delta[i] - (2 * Desv[i])
              if (LS > 0 & LI > 0 || LS < 0 & LI < 0){
                Vec[i] = "P. Necesario"
              } else {
                Vec[i]= "P. Innecesario"
              }
            }
            Vec = data.frame(Vec)
            colnames(Vec) = c('Test de parsimonia')
            
          }
          
          Coeff <- Model$coef
          Var.coeff <- diag(Model$var.coef)
          Parsimonio(Delta = Coeff, Delta.var = Var.coeff)
        #Pronostico
          # Datos pronosticados
          Pron = forecast(Model, h = P)
          Upper = Pron$upper
          Lower = Pron$lower
          xp = Pron$x
          
          # Media
          Mid = Pron$mean
          Mid1 = data.frame(matrix(data = NA, nrow = P, ncol = 2))
          colnames(Mid1) = c('T', 'Media')
          Mid1[, 1] = seq((length(xp) + 1), (length(xp) + P), 1)
          Mid1[, 2] = Mid
          
          # Upper 80% confidence level
          Up80 = Upper[, 1]
          Up.80 = data.frame(matrix(data = NA, nrow = P, ncol = 2))
          colnames(Up.80) = c('T', 'Up80')
          Up.80[, 1] = seq((length(xp) + 1), (length(xp) + P), 1)
          Up.80[, 2] = Up80
          
          # Lower 80% confidence level
          Low80 = Lower[, 1]
          Low.80 = data.frame(matrix(data = NA, nrow = P, ncol = 2))
          colnames(Low.80) = c('T', 'Low80')
          Low.80[, 1] = seq((length(xp) + 1), (length(xp) + P), 1)
          Low.80[, 2] = Low80
          
          # Upper 95% confidence level
          Up95 = Upper[, 2]
          Up.95 = data.frame(matrix(data = NA, nrow = P, ncol = 2))
          colnames(Up.95) = c('T', 'Up95')
          Up.95[, 1] = seq((length(xp) + 1), (length(xp) + P), 1)
          Up.95[, 2] = Up95
          
          # Lower 95% confidence level
          Low95 = Lower[, 2]
          Low.95 = data.frame(matrix(data = NA, nrow = P, ncol = 2))
          colnames(Low.95) = c('T', 'Low95')
          Low.95[, 1] = seq((length(xp) + 1), (length(xp) + P), 1)
          Low.95[, 2] = Low95
          
          # Datos muestrales
          XII = data.frame(matrix(data = NA, nrow = length(xp), ncol = 2))
          colnames(XII) = c('T', 'D.Obs')
          XII[, 1] = seq(1, length(xp), 1)
          XII[, 2] = xp
          
          G_XII <- ggplot(data = XII, mapping = aes(x = T, y = D.Obs, group = 12)) +
            geom_line(colour = 'black', size = 0.8) +
            labs(title = 'Datos simulados - Pron??stico', x = 'T', y = 'ICCV') +
            geom_line(data = Mid1, mapping = aes(x = T, y = Media, group = 12), colour = 'blue', size = 1) +
            geom_line(data = Up.80, mapping = aes(x = T, y = Up80, group = 12), colour = '#ff8c00', size = 1) +
            geom_line(data = Low.80, mapping = aes(x = T, y = Low80, group = 12), colour = '#ff8c00', size = 1) +
            geom_line(data = Up.95, mapping = aes(x = T, y = Up95, group = 12), colour = 'red', size = 1) +
            geom_line(data = Low.95, mapping = aes(x = T, y = Low95, group = 12), colour = 'red', size = 1) +
            theme(legend.position="top")
          
          
      
        # Salidas
        print(sprintf('La media del proceso es %s', Media))
        
        print("La serie con media no nula")
        print(Ztmu)
        print(data.frame(X))
        print("")
        print("Autocorrelacion de un AR 3 muestral")
        print(data.frame(ACF))
        print("")
        print("Autocorrelacion parcialde un AR 3 muestral")
        print(data.frame(PACF))
        print("")
        
        # EstimaciC3n comando ARIMA
        print("Estimacion de un ARIMA")
        print(sprintf('Los datos se comportan como un AR puro de orden %s.', AR_P))
        print(sprintf('Los datos se comportan como un MA puro de orden %s.', MA_Q))
        print("")

        #EstimaciC3n del modelo
        print("Estimacion del modelo")
        print(data.frame(At))
        print("")

        #Media 0 de Residuales
        print("Supuesto de media 0 en los residuales")
        if (Stat < Vc){
          print(sprintf('No hay suficiente evidencia para rechazar media 0 en los erorres con un alpha de %s', alpha))
        } else {
          print(sprintf('Hay suficiente evidencia para rechazar media 0 en los erorres con un alpha de %s', alpha))
        }
        print("")

        #Indepencia de residuales
        print("Supuesto de independencia en los residuales")
        if (Q < Vc){
          print(sprintf('No hay suficiente evidencia para rechazar independencia en los erorres con un alpha de %s', alpha))
        } else {
          print(sprintf('Hay suficiente evidencia para rechazar independencia en los erorres con un alpha de %s', alpha))
        }
        print("")

        print(sprintf("Q' = %s", Q))
        print(sprintf("El valor critico es :%s", Vc))
        print("")
              

        #Normalidad de residuales
        print("Supuesto de normalidad en los residuales")
        if(JB < VcJB){
          print(sprintf('No hay suficiente evidencia para rechazar normalidad en los errores a un nivel de significancia alpha de: %s', alpha))
        } else {
          print(sprintf('Hay suficiente evidencia para rechazar normalidad en los errores a un nivel de significancia alpha de: %s', alpha))
        }
        print("")
        #No existencia de observaciones aberrantes
        print("Supuesto de no existencia de observaciones aberrantes ")
        print(Obs.aberr(X = At))
        print("")

        #Parsimonio
        print("Supuesto de parsimonio")
        # print(Vec)
        print(sprintf("El coeficiente es %s",Coeff))
        print(sprintf("La varianza del coeficiente es: %s", Var.coeff))
        print(Parsimonio(Delta = Coeff, Delta.var = Var.coeff))
        print("")
        H = forecast(Model, h = 24)
        plot(H)

        # #pronC3stico
        H = forecast(Model, h = 24)
         
        return(G_Zt/ G_ACF / G_PACF)

       
      } else {
        print(sprintf('El proceso no es estacionario al tener sus factores |%s| > 1', round(G, 4)))
      }
      
      
    })
    
    output$Muestral = renderPlot(
      muest()
    )
    
    output$Muestra = renderPrint(
      muest()
    )
    muest2 =reactive({
      #Inputs que se usarC!n 
      Phi = c(input$Phi1, input$Phi2, input$Phi3)               # Phi's poblacionales para un modelo AR[3]
      alpha = 0.05                                             # Nivel de significancia       
      K = input$L                                              # TamaC1o de la muestra 
      Mu = input$Mean                                  # Media poblacional diferente de 0
      N = input$Random
      
      # Datos a pronosticar
      P = input$P
      
      # SimulaciC3n de ruido blanco
      set.seed(0214)                          # Semilla para la generaciC3n de datos aleatorios
      at = rnorm(N, mean = 0, sd = 1)
      sigma_a = 1
      # CC!lculo de los factores con su parte imaginaria
      Phi = matrix(data = c(-1, c(Phi)), nrow = 4, ncol = 1) * -1
      G = matrix(data = c(as.complex(polyroot(Phi)) ** -1), nrow = 3, ncol = 1)               # Polyroot devuelve las raices, ** -1 = Factores
      for (i in 1 : length(G)){
        if (sum(round(Im(G), 10)) == 0){
          G = G
        } else {
          Re_G = Re(G)
          Im_G = Im(G)
          for (i in 1 : length(G)){                                                               # Convierten los nC:meros complejos en reales
            G[i] = sqrt(Re_G[i] ** 2 + Im_G[i] ** 2)
          }
        }
      }
      
      # IdentificaciC3n estacionareidad
      NG = length(G)
      G = round(Re(G), 10)
      Check = matrix(data = NA, nrow = NG, ncol = 1)
      for (i in 1 : NG){
        if (abs(G[i]) < 1){
          Check[i] = 1
        } else {
          Check[i] = 0
        }
      }
      
      # IdentificaciC3n de factores iguales
      Check_fact = matrix(data = sort(G, decreasing = TRUE), nrow = NG, ncol = 2)
      for (i in 2 : NG){
        Check_fact[1, 2] = 1
        if (Check_fact[i, 1] == Check_fact[(i - 1), 1]){
          Check_fact[i, 2] = Check_fact[(i - 1), 2] + 1
        } else {
          Check_fact[i, 2] = 1
        }
      }
      
      # CC!lculo de los valores iniciales
      if (sum(Check) == 3){
        
        # CC!lculo de los Rho
        Phi = Phi * -1
        Rho1 = (Phi[2] + (Phi[3] * Phi[4])) / (1 - Phi[3] - (Phi[4] * Phi[2]) - Phi[4] ** 2)
        Rho2 = Rho1 * (Phi[2] + Phi[4]) + Phi[3]
        Rho = matrix(data = c(1, Rho1, Rho2), nrow = NG, ncol = 1)
        
        # CC!lculo de la ACF con y sin factores distintos
        if (sum(Check_fact[, 2]) == NG){
          Ginv = rbind(c(1, 1, 1), c(G), c(G ** 2))               # CC!lculo de los valores iniciales con factores diferentes
          A = as.matrix(solve(Ginv) %*% Rho)
          print('Todos los factores son diferentes.')
          
          # ConstrucciC3n ACF (SoluciC3n de la funciC3n generadora)
          ACF = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
          ACF[, 1] = seq(0, K, 1)
          colnames(ACF) = c('N', 'Rho')
          for (i in 0 : N){
            ACF[i + 1, 2] = as.numeric((A[1] * (G[1] ** i)) + (A[2] * (G[2] ** i)) + (A[3] * (G[3] ** i)))
          }
          ACF = with(ACF, data.frame(N, Rho))
          
        }
        if (sum(Check_fact[, 2]) == (NG + 1)){
          Ginv = rbind(c(1, 1, 0), c(G), c(G[1] ** 2, G[2] ** 2, 2 * G[3] ** 2))               # CC!lculo de los valores iniciales con factores iguales
          A = as.matrix(solve(Ginv) %*% Rho)
          print('Dos factores son iguales.')
          
          # ConstrucciC3n ACF (SoluciC3n de la funciC3n generadora)
          ACF = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
          ACF[, 1] = seq(0, N, 1)
          colnames(ACF) = c('N', 'Rho')
          FG = matrix(data = NA, ncol = NG, nrow = (N + 1))
          for (j in 1 : NG){
            for (i in 0 : N){
              if (Check_fact[j, 2] == 1){
                FG[(i + 1), j] = A[j] * (G[j] ** i)
              } else {
                FG[(i + 1), j] = A[j] * (G[j] ** i) * i
              }
            }
          }
          
          for (i in 0 : N){
            ACF[i + 1, 2] = as.numeric(sum(FG[(i + 1), ]))
          }
          ACF = with(ACF, data.frame(N, Rho))
          
        }
        if (sum(Check_fact[, 2]) == (NG + 3)){
          G_3 = c(G[1], (2 * G[2]), (4 * G[3]))
          Ginv = rbind(c(1, 0, 0), c(G), c(G_3 ** 2))              # CC!lculo de los valores iniciales con factores iguales
          A = as.matrix(solve(Ginv) %*% Rho)
          print('Tres factores son iguales.')
          
          # ConstrucciC3n ACF (SoluciC3n de la funciC3n generadora)
          ACF = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
          ACF[, 1] = seq(0, K, 1)
          colnames(ACF) = c('N', 'Rho')
          FG = matrix(data = NA, ncol = NG, nrow = (N + 1))
          for (j in 1 : NG){
            FG[(i + 1), j] = A[j] * (G[j] ** i) * i
          }
          for (i in 0 : N){
            ACF[i + 1, 2] = as.numeric(sum(FG[(i + 1), ]))
          }
          ACF = with(ACF, data.frame(N, Rho))
          
        }
        
        # ConstrucciC3n PACF (Usando Cramer)
        PACF = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
        PACF[, 1] = seq(1, (N + 1), 1)
        colnames(PACF) = c('N', 'Phi')
        for (k in 1 : N){
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
        PACF = with(PACF, data.frame(N, Phi))
        
        # Varianza del proceso (Gamma_0_1_2)
        F1 = c((1 - (Phi[2] ** 2) - (Phi[3] ** 2) - (Phi[4] * 2)), (-2 * Phi[3] * (Phi[2] + Phi[4])), (-2 * Phi[2] * Phi[4]))
        F2 = c((-Phi[2]), (1 - Phi[3]), (-Phi[4]))
        F3 = c((-Phi[3]), -(Phi[2] + Phi[4]), (1))
        C2 = matrix(data = c(sigma_a, 0, 0), ncol = 1, nrow = 3)
        Gamma = matrix(data = rbind(F1, F2, F3), ncol = 3, nrow = 3)
        Gamma_1_2_3 = matrix(data = NA, ncol = 2, nrow = 3)
        colnames(Gamma_1_2_3) = c('Gamma', 'Valor')
        Gamma_1_2_3[, 1] = seq(0, 2, 1)
        Gamma_1_2_3[, 2] = solve(Gamma) %*% C2
        
        
        
        # GrC!fica conjunta de ACF y PACF
        Cv = qt(p = (alpha / 2), df = (N - 1), lower.tail = FALSE) / sqrt(N)        # Valor crC-tico del intervalo de confianza
        G_ACF_3 <- ggplot(data = ACF, mapping = aes(x = N, y = Rho)) +
          labs(title = 'ACF AR[3] Teorica') +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = N, yend = 0)) +
          geom_hline(aes(yintercept = Cv), color="blue",
                     linetype = "dashed") +
          geom_hline(aes(yintercept = -Cv), color="blue",
                     linetype = "dashed")
        
        G_PACF_3 <- ggplot(data = PACF, mapping = aes(x = N, y = Phi)) +
          labs(title = 'PACF AR[3] Teorica') +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = N, yend = 0)) +
          geom_hline(aes(yintercept = Cv), color="blue",
                     linetype = "dashed") +
          geom_hline(aes(yintercept = -Cv), color="blue",
                     linetype = "dashed")
        
        Psi_2 = PACF[1, 2]
        Psi_3 = (PACF[1, 2] ** 2) + PACF[2, 2]
        Psi = matrix(data = c(1, Psi_2, Psi_3), nrow = 3, ncol = 1)
        
        # CC!lculo de la ACF con y sin factores distintos
        if (sum(Check_fact[, 2]) == NG){
          Ginv = rbind(c(1, 1, 1), c(G), c(G ** 2))               # CC!lculo de los valores iniciales con factores diferentes
          A = as.matrix(solve(Ginv) %*% Psi)
          
          # ConstrucciC3n ACF (SoluciC3n de la funciC3n generadora)
          ACF_MA = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
          ACF_MA[, 1] = seq(0, N, 1)
          colnames(ACF_MA) = c('N', 'Psi')
          for (i in 0 : N){
            ACF_MA[i + 1, 2] = as.numeric((A[1] * (G[1] ** i)) + (A[2] * (G[2] ** i)) + (A[3] * (G[3] ** i)))
          }
          ACF_MA = with(ACF_MA, data.frame(N, Psi))
          
        } else {
          Ginv = rbind(c(1, 1, 0), c(G), c(G ** 2))               # CC!lculo de los valores iniciales con factores iguales
          A = as.matrix(solve(Ginv) %*% Psi)
          
          # ConstrucciC3n ACF (SoluciC3n de la funciC3n generadora)
          ACF_MA = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
          ACF_MA[, 1] = seq(0, N, 1)
          colnames(ACF_MA) = c('N', 'Psi')
          FG = matrix(data = NA, ncol = NG, nrow = (N + 1))
          for (j in 1 : NG){
            for (i in 0 : N){
              if (Check_fact[j, 2] == 1){
                FG[(i + 1), j] = A[j] * (G[j] ** i)
              } else {
                FG[(i + 1), j] = A[j] * (G[j] ** i) * i
              }
            }
          }
          for (i in 0 : N){
            ACF_MA[i + 1, 2] = as.numeric(sum(FG[(i + 1), ]))
          }
          ACF_MA = with(ACF_MA, data.frame(N, Psi))
          
        }
        
        # GrC!fica conjunta de ACF y PACF
        G_MA <- ggplot(data = ACF_MA, mapping = aes(x = N, y = Psi)) +
          labs(title = 'MA[inf] AR[3] TeC3rico') +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = N, yend = 0)) +
          geom_hline(aes(yintercept = Cv), color="blue",
                     linetype = "dashed") + 
          geom_hline(aes(yintercept = -Cv), color="blue",
                     linetype = "dashed")
        
        ## Muestral
        # Matriz de datos de Zt 
        Zt = data.frame(matrix(data = NA, ncol = 2, nrow = N))
        colnames(Zt) = c('N', 'Zt')
        Zt[, 1] = seq(1, N, 1)
        Zt[1 : 3, 2] = A
        for (i in (length(A) + 1) : N){
          Zt[i, 2] = (PACF[1, 2] * Zt[i - 1, 2]) + (PACF[2, 2] * Zt[i - 2, 2]) + (PACF[3, 2] * Zt[i - 3, 2]) + at[i]
        }
        
        # Matriz de datos de Zt con media no nula
        Ztmu = data.frame(matrix(data = NA, ncol = 2, nrow = N))
        colnames(Ztmu) = c('N', 'ZtMu')
        Ztmu[, 1] = seq(1, N, 1)
        Ztmu[, 2] = Zt[, 2] + Mu
        Media = (sum(Ztmu[, 2]) / N)
        
        # GrC!fica de la serie
        G_Zt <- ggplot(data = Ztmu, mapping = aes(x = N, y = ZtMu, group = 12)) +
          geom_line(colour = 'red', size = 1) +
          labs(title = 'Datos simulados de un AR[3] con media no nula', x = 'n', y = 'Valor de Zt') +
          geom_hline(aes(yintercept = Media))
        
        # RedefiniciC3n de los datos
        X = Ztmu[, 2]
        
        
        # ConstrucciC3n de los vectores y parC!metros
        T = length(X)
        
        # CC!lculo de variables
        Zb = sum(X) / T                                                             # Media muestral
        Cv = qt(p = (alpha / 2), df = (T - 1), lower.tail = FALSE) / sqrt(T)        # Valor crC-tico del intervalo de confianza
        
        # FunciC3n de autocovarianza
        Gamma = matrix(data = NA)                              
        for (i in 0 : N){
          Gamma[i + 1] = sum((X[1 : (N - i)] - Zb) * (X[(i + 1) : N] - Zb)) / (N - 1)
        }
        
        # ConstrucciC3n ACF
        ACF = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))
        ACF[, 1] = seq(0, N, 1)
        colnames(ACF) = c('N', 'Rho')
        for (i in 0 : N){
          ACF[i + 1, 2] = Gamma[i + 1] / Gamma[1]
        }
        ACF = na.omit(ACF)
        ACF = with(ACF, data.frame(N, Rho))
        
        # ConstrucciC3n PACF
        PACF = data.frame(matrix(data = NA, nrow = (N + 1), ncol = 2))  
        PACF[, 1] = seq(1, (N + 1), 1)
        colnames(PACF) = c('N', 'Phi')                            
        for (k in 1 : N){
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
        PACF = PACF[-(N + 1),]
        PACF = na.omit(PACF)
        PACF = with(PACF, data.frame(N, Phi))
        
        # GrC!fica conjunta de ACF y PACF
        G_ACF <- ggplot(data = ACF, mapping = aes(x = N, y = Rho)) +
          labs(title = 'ACF') +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = N, yend = 0)) +
          geom_hline(aes(yintercept = Cv), color="blue",
                     linetype = "dashed", size = 0.5) + 
          geom_hline(aes(yintercept = -Cv), color="blue",
                     linetype = "dashed", size = 0.5)
        
        G_PACF <- ggplot(data = PACF, mapping = aes(x = N, y = Phi)) +
          labs(title = 'PACF') +
          geom_hline(aes(yintercept = 0)) +
          geom_segment(mapping = aes(xend = N, yend = 0)) + 
          geom_hline(aes(yintercept = Cv), color="blue",
                     linetype = "dashed", size = 0.5) + 
          geom_hline(aes(yintercept = -Cv), color="blue",
                     linetype = "dashed", size = 0.5)
        
        
        # IdentificaciC3n del modelo AR
        CheckAR = matrix(data = NA, ncol = 1, nrow = N)
        for(i in 1 : (N - 1)) {
          CheckAR[i] = between(0, (PACF[i, 2] - Cv), (PACF[i, 2] + Cv))
          if(CheckAR[i] == 1){
            next
          }
          AR_P = c(i)
        }
        
        # IdentificaciC3n del modelo MA
        CheckMA = matrix(data = NA, ncol = 1, nrow = N)
        for(i in 1 : (N - 1)){
          CheckMA[i] = between(0, (ACF[i, 2] - Cv * (sqrt(1 + 2 * sum(ACF[(i + 1), 2] ** 2)))), (ACF[i, 2] + Cv) * (sqrt(1 + 2 * sum(ACF[(i + 1), 2] ** 2))))
          if(CheckMA[i] == 1){
            next
          }
          MA_Q = c(i)
        }
        # Inputs
        AR_P = 2
        MA_Q = 1
        J_opt = 0                       # Es 0 porque la serie es estable en nivel desde la simulaciC3n
        
        
        Model = arima(X, order = c(AR_P, J_opt, MA_Q), method = c('ML'))
        At = Model$residuals
        
        # Media 0 de residuales
        t = J_opt + AR_P + 1
        T = length(At)
        Ma = sum(At[t : T]) / (T - J_opt - AR_P)
        Sigma_at = sqrt(t(At - Ma) %*% (At - Ma) / (T - J_opt - AR_P - MA_Q))
        Stat = abs((sqrt(T - J_opt - AR_P) * Ma) / (Sigma_at)) 
        Vc = qt(p = alpha, df = (T - AR_P - MA_Q - J_opt), lower.tail = FALSE)
        
        #Varianza constante de residuales
        At1 = data.frame(matrix(data = NA, nrow = T, ncol = 2))  
        colnames(At1) = c('T', 'D.At')
        At1[, 1] = seq(1, T, 1)
        At1[, 2] = At
        At1 = with(At1, data.frame(T, D.At))
        At_sigma = 1 * Sigma_at
        
        #GrC!fica
        G_At1 <- ggplot(data = At1, mapping = aes(x = T, y = D.At)) +
          geom_line(colour = 'red', size = 1) +
          labs(title = 'Inspecci??n de residuales por varianza [At]', y = 'At') +
          geom_hline(aes(yintercept = At_sigma), color="blue",
                     linetype = "dashed", size = 1) +
          geom_hline(aes(yintercept = -At_sigma), color="blue",
                     linetype = "dashed", size = 1)
        
        #Independencia de residuales
        Rho_at = matrix(data = NA, nrow = (T - N), ncol = 1)
        for (i in 1 : (T - N)){
          Rho_at[i] = sum(At[i] * At[(i + t)]) / sum(At ** 2)
        }
        
        # Prueba de Ljung - Box
        Q = (T - J_opt - AR_P) * (T - J_opt - AR_P + 2) * sum((sum(Rho_at) ** 2) / (T - J_opt - AR_P - K)) 
        Vc = qchisq(p = alpha, df = (K - AR_P - MA_Q))
        
        # Normalidad de residuales 
        #Prueba Jarque Bera (Normalidad en el error)
        S = (sum((At - mean(At)) ** 3) / T) / (sqrt(as.numeric((t(At) %*% At) / T)) ** 3)           # AsimetrC-a
        K1 = (sum((At - mean(At)) ** 4) / T) / (sqrt(as.numeric((t(At) %*% At) / T)) ** 4)          # Kurtosis
        JB = (1 / T) * ((S ** 2 / 6) + ((K1 - 3) ** 2 / 24))
        VcJB = qchisq(alpha, df = 2, lower.tail = FALSE)
        
        At1 = data.frame(matrix(data = NA, nrow = T, ncol = 2))
        colnames(At1) = c('N', 'At')
        At1[, 1] = seq(1, T, 1)
        At1[, 2] = At
        A_h = ggplot(At1, aes(x = At)) + 
          labs(title = 'Histograma de At') +
          geom_histogram(aes(y = ..density..), colour = 'darkblue', fill = 'lightblue', alpha = 0.5, binwidth = 0.5) +
          geom_density(alpha = 0.2, fill = "#FF6666") +
          geom_vline(aes(xintercept= mean(At)), color = "blue", linetype = "dashed", size = 0.7) 
        
        #No existencia de observaciones aberrantes 
        Obs.aberr <- function(X){    
          T = length(X)
          Xmean = sum(X) / T
          Xsd = sd(X)  
          LS = Xmean + (3 * Xsd)
          LI = Xmean - (3 * Xsd)  
          D.aberr = c()
          
          for (i in 1 : T){
            if (X[i] < LS & X[i] > LI){
              D.aberr[i] = NA
            } else {
              D.aberr[i] = X[i]
            }
          }
          
          D.aberr = na.omit(D.aberr)  
          
          # Salidas
          print(sprintf('Hay un total de %s observaciones aberrantes', length(D.aberr)))
        }
        
        # GrC!fica de los datos
        At1 = data.frame(matrix(data = NA, nrow = T, ncol = 2))  
        colnames(At1) = c('T', 'D.At')
        At1[, 1] = seq(1, T, 1)
        At1[, 2] = At
        At1 = with(At1, data.frame(T, D.At))
        At_sigma = 3 * Sigma_at
        
        G_At1 <- ggplot(data = At1, mapping = aes(x = T, y = D.At)) +
          geom_line(colour = '#ff0202', size = 1) +
          labs(title = 'At') +
          geom_hline(aes(yintercept = At_sigma), color="blue",
                     linetype = "dashed", size = 1) +
          geom_hline(aes(yintercept = -At_sigma), color="blue",
                     linetype = "dashed", size = 1)
        
        #Parsimonio
        Parsimonio = function(Delta, Delta.var){
          Desv = sqrt(Delta.var)
          Vec = Delta
          for (i in 1 : length(Delta)){
            LS = Delta[i] + (2 * Desv[i])
            LI = Delta[i] - (2 * Desv[i])
            if (LS > 0 & LI > 0 || LS < 0 & LI < 0){
              Vec[i] = "P. Necesario"
            } else {
              Vec[i]= "P. Innecesario"
            }
          }
          Vec = data.frame(Vec)
          colnames(Vec) = c('Test de parsimonia')
          
        }
        
        Coeff <- Model$coef
        Var.coeff <- diag(Model$var.coef)
        Parsimonio(Delta = Coeff, Delta.var = Var.coeff)
        #Pronostico
        # Datos pronosticados
        Pron = forecast(Model, h = P)
        Upper = Pron$upper
        Lower = Pron$lower
        xp = Pron$x
        
        # Media
        Mid = Pron$mean
        Mid1 = data.frame(matrix(data = NA, nrow = P, ncol = 2))
        colnames(Mid1) = c('T', 'Media')
        Mid1[, 1] = seq((length(xp) + 1), (length(xp) + P), 1)
        Mid1[, 2] = Mid
        
        # Upper 80% confidence level
        Up80 = Upper[, 1]
        Up.80 = data.frame(matrix(data = NA, nrow = P, ncol = 2))
        colnames(Up.80) = c('T', 'Up80')
        Up.80[, 1] = seq((length(xp) + 1), (length(xp) + P), 1)
        Up.80[, 2] = Up80
        
        # Lower 80% confidence level
        Low80 = Lower[, 1]
        Low.80 = data.frame(matrix(data = NA, nrow = P, ncol = 2))
        colnames(Low.80) = c('T', 'Low80')
        Low.80[, 1] = seq((length(xp) + 1), (length(xp) + P), 1)
        Low.80[, 2] = Low80
        
        # Upper 95% confidence level
        Up95 = Upper[, 2]
        Up.95 = data.frame(matrix(data = NA, nrow = P, ncol = 2))
        colnames(Up.95) = c('T', 'Up95')
        Up.95[, 1] = seq((length(xp) + 1), (length(xp) + P), 1)
        Up.95[, 2] = Up95
        
        # Lower 95% confidence level
        Low95 = Lower[, 2]
        Low.95 = data.frame(matrix(data = NA, nrow = P, ncol = 2))
        colnames(Low.95) = c('T', 'Low95')
        Low.95[, 1] = seq((length(xp) + 1), (length(xp) + P), 1)
        Low.95[, 2] = Low95
        
        # Datos muestrales
        XII = data.frame(matrix(data = NA, nrow = length(xp), ncol = 2))
        colnames(XII) = c('T', 'D.Obs')
        XII[, 1] = seq(1, length(xp), 1)
        XII[, 2] = xp
        
        G_XII <- ggplot(data = XII, mapping = aes(x = T, y = D.Obs, group = 12)) +
          geom_line(colour = 'black', size = 0.8) +
          labs(title = 'Datos simulados - Pron??stico', x = 'T', y = 'ICCV') +
          geom_line(data = Mid1, mapping = aes(x = T, y = Media, group = 12), colour = 'blue', size = 1) +
          geom_line(data = Up.80, mapping = aes(x = T, y = Up80, group = 12), colour = '#ff8c00', size = 1) +
          geom_line(data = Low.80, mapping = aes(x = T, y = Low80, group = 12), colour = '#ff8c00', size = 1) +
          geom_line(data = Up.95, mapping = aes(x = T, y = Up95, group = 12), colour = 'red', size = 1) +
          geom_line(data = Low.95, mapping = aes(x = T, y = Low95, group = 12), colour = 'red', size = 1) +
          theme(legend.position="top")
        
        
        
        # Salidas
        print(sprintf('La media del proceso es %s', Media))
        
        print("La serie con media no nula")
        print(Ztmu)
        print(data.frame(X))
        print("")
        print("Autocorrelacion de un AR 3 muestral")
        print(data.frame(ACF))
        print("")
        print("Autocorrelacion parcialde un AR 3 muestral")
        print(data.frame(PACF))
        print("")
        
        # EstimaciC3n comando ARIMA
        print("Estimacion de un ARIMA")
        print(sprintf('Los datos se comportan como un AR puro de orden %s.', AR_P))
        print(sprintf('Los datos se comportan como un MA puro de orden %s.', MA_Q))
        print("")
        
        #EstimaciC3n del modelo
        print("Estimacion del modelo")
        print(data.frame(At))
        print("")
        
        #Media 0 de Residuales
        print("Supuesto de media 0 en los residuales")
        if (Stat < Vc){
          print(sprintf('No hay suficiente evidencia para rechazar media 0 en los erorres con un alpha de %s', alpha))
        } else {
          print(sprintf('Hay suficiente evidencia para rechazar media 0 en los erorres con un alpha de %s', alpha))
        }
        print("")
        
        #Indepencia de residuales
        print("Supuesto de independencia en los residuales")
        if (Q < Vc){
          print(sprintf('No hay suficiente evidencia para rechazar independencia en los erorres con un alpha de %s', alpha))
        } else {
          print(sprintf('Hay suficiente evidencia para rechazar independencia en los erorres con un alpha de %s', alpha))
        }
        print("")
        
        print(sprintf("Q' = %s", Q))
        print(sprintf("El valor critico es :%s", Vc))
        print("")
        
        
        #Normalidad de residuales
        print("Supuesto de normalidad en los residuales")
        if(JB < VcJB){
          print(sprintf('No hay suficiente evidencia para rechazar normalidad en los errores a un nivel de significancia alpha de: %s', alpha))
        } else {
          print(sprintf('Hay suficiente evidencia para rechazar normalidad en los errores a un nivel de significancia alpha de: %s', alpha))
        }
        print("")
        #No existencia de observaciones aberrantes
        print("Supuesto de no existencia de observaciones aberrantes ")
        print(Obs.aberr(X = At))
        print("")
        
        #Parsimonio
        print("Supuesto de parsimonio")
        # print(Vec)
        print(sprintf("El coeficiente es %s",Coeff))
        print(sprintf("La varianza del coeficiente es: %s", Var.coeff))
        print(Parsimonio(Delta = Coeff, Delta.var = Var.coeff))
        print("")
        H = forecast(Model, h = 24)
        plot(H)
        
        # #pronC3stico
        H = forecast(Model, h = 24)
        
        return(G_At1/A_h/G_XII)
        
        
      } else {
        print(sprintf('El proceso no es estacionario al tener sus factores |%s| > 1', round(G, 4)))
      }
      
      
    })
    
    output$muestra2 = renderPlot(
      muest2()
    )
      
  }
  
  
  


  


shinyApp(ui=ui, server=server)

