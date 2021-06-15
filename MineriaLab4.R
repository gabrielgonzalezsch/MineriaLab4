require(bnlearn)
library(bnlearn)
library(bnviewer)


########################### 1. Lectura de Datos ###########################
data <- data.frame(alarm)

########################### 2. Preprocesamiento ###########################


cantidad.nulos <-  sapply(data, function(x) sum(is.na(data))) 
#No hay NA en el dataset

summary(data)
str(data)
#2000 observaciones
#37 variables
#Datos categoricos (no hay números)
#Datos son factores de 2,3 y 4 niveles según la variable

bn_df <- data

res <- hc(bn_df) #Algoritmo Hill-Climbing

# Visualizacion de las redes bayesianas 

viewer(res,
       bayesianNetwork.width = "100%",
       bayesianNetwork.height = "80vh",
       bayesianNetwork.layout = "layout_with_sugiyama",
       bayesianNetwork.title="Discrete Bayesian Network - Alarm",
       bayesianNetwork.subtitle = "Monitoring of emergency care patients",
       bayesianNetwork.footer = "Fig. 1 - Layout with Sugiyama"
)

viewer(res,
       bayesianNetwork.width = "100%",
       bayesianNetwork.height = "80vh",
       bayesianNetwork.layout = "layout_on_grid",
       bayesianNetwork.title="Discrete Bayesian Network - Alarm",
       bayesianNetwork.subtitle = "Monitoring of emergency care patients",
       bayesianNetwork.footer = "Fig. 2 - Layout on grid",
       
       node.colors = list(background = "#f4bafd",
                          border = "#2b7ce9",
                          highlight = list(background = "#97c2fc",
                                           border = "#2b7ce9"))
)



bl<-data.frame("M..Work","Family") #Lista negra de relaciones (Origen, Destino)
res <- hc(bn_df,blacklist = bl)
plot(res)

print(res)

sc<-score(res,bn_df) # BIC por default
print(sc)


fittedbn <- bn.fit(res, data = bn_df) # Se obtiene la tabla de probabilidades condicionales mediante EM. (Máxima Expectación, propagación de la evidencia)
print(fittedbn$CVP) #se obtiene la información respecto del nodo Proteins


cpquery(fittedbn, event = (Proteins=="<3"), evidence = ( Smoking=="no") ) 

cpquery(fittedbn, event = (Pressure==">140"), evidence = ( Proteins=="<3" ) )

require(bnlearn)
bn_df <- data.frame(alarm)


res <- hc(bn_df)
plot(res)
sc<-score(res,bn_df)
print(sc)


fittedbn <- bn.fit(res, data = bn_df)
print(fittedbn$E)


cpquery(fittedbn, event = (B=="yes"), evidence = ( S=="no") ) 

cpquery(fittedbn, event = (B=="yes"), evidence = ( S=="yes") )

bn_df <- data.frame(alarm)
res <- mmhc(bn_df)
plot(res)
sc<-score(res,bn_df)
print(sc)


fittedbn <- bn.fit(res, data = bn_df) #idéntico al anterior
print(fittedbn$E)

cpquery(fittedbn, event = (B=="yes"), evidence = ( S=="no") )

cpquery(fittedbn, event = (B=="yes"), evidence = ( S=="yes") )

bn_df <- data.frame(alarm)
res <- mmpc(bn_df)
plot(res)

fittedbn <- bn.fit(res, data = bn_df) #no hay direccionalidad en los arcos



#bl <- matrix(c("Age","GoodStudent","Age","SeniorTrain","SocioEcon","OtherCar","ThisCarDam","Theft","MakeModel","VehicleYear","SeniorTrain","Age","Mileage","Mileage"),ncol=2,byrow = TRUE)