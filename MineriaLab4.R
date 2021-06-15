require(bnlearn)
library(bnlearn)
library(bnviewer)



########################### Funciones ########################### 
graficar1 <- function(res){
        # Visualizacion de las redes bayesianas 
        viewer(res,
               bayesianNetwork.width = "100%",
               bayesianNetwork.height = "80vh",
               bayesianNetwork.layout = "layout_with_sugiyama",
               bayesianNetwork.title="Discrete Bayesian Network - Alarm",
               bayesianNetwork.subtitle = "Monitoring of emergency care patients",
               bayesianNetwork.footer = "Fig. 1 - Layout with Sugiyama"
        )   
}

graficar2 <- function(res){
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
        
        
}




########################### 1. Lectura de Datos ###########################
data <- data.frame(alarm)

modelstring = paste0("[HIST|LVF][CVP|LVV][PCWP|LVV][HYP][LVV|HYP:LVF][LVF]",
                     "[STKV|HYP:LVF][ERLO][HRBP|ERLO:HR][HREK|ERCA:HR][ERCA][HRSA|ERCA:HR][ANES]",
                     "[APL][TPR|APL][ECO2|ACO2:VLNG][KINK][MINV|INT:VLNG][FIO2][PVS|FIO2:VALV]",
                     "[SAO2|PVS:SHNT][PAP|PMB][PMB][SHNT|INT:PMB][INT][PRSS|INT:KINK:VTUB][DISC]",
                     "[MVS][VMCH|MVS][VTUB|DISC:VMCH][VLNG|INT:KINK:VTUB][VALV|INT:VLNG]",
                     "[ACO2|VALV][CCHL|ACO2:ANES:SAO2:TPR][HR|CCHL][CO|HR:STKV][BP|CO:TPR]")

dag = model2network(modelstring)

########################### 2. Preprocesamiento ###########################


cantidad.nulos <-  sapply(data, function(x) sum(is.na(data))) 
#No hay NA en el dataset

summary(data)
str(data)
#2000 observaciones
#37 variables
#Datos categoricos (no hay números)
#Datos son factores de 2,3 y 4 niveles según la variable



#Lista negra


bl <- matrix(c("HREK","HRBP","PCWP","LVV","PRSS","VTUB","ERLO","HR","HRBP","HR","KINK","DISC","INT","KINK","INT","DISC",
               "PRSS","KINK","MINV","VALV","MINV","INT","MINV","VLNG","MINV","VTUB"),ncol=2,byrow=TRUE)

wl <- matrix(c("CCHL","HR","CVP","CO","PAP","PMB"),ncol=2,byrow=TRUE)


########################### 3. Aplicar Método ###########################

bn_df <- data

# Algoritmo Hill-Climbing
resHC <- hc(bn_df,blacklist = bl, whitelist = wl)      
print(resHC)
#Graficar
graficar1(resHC)
graficar2(resHC)
#Obtener BIC
scHC <-score(resHC,bn_df )  
print(scHC)
#Propagación de la evidencia: EM (Ejecuación de Máxima Exprectación)
fittedbn.HC <- bn.fit(resHC, data = bn_df)

# Algoritmo MaxMin Hill-Climbing
resMMHC <- mmhc(bn_df,blacklist = bl, whitelist = wl)  
print(resMMHC)
#Graficar
graficar1(resMMHC)
graficar2(resMMHC)
#Obtener BIC
scMMHC <- score(resMMHC,bn_df)
print(scMMHC)
#Propagación de la evidencia: EM (Ejecuación de Máxima Exprectación)
fittedbn.MMHC <- bn.fit(resMMHC, data = bn_df)

# Algoritmo MaxMin Parents and Children
resMMPC <- mmpc(bn_df,blacklist = bl, whitelist = wl)  
print(resMMPC)
#Graficar
graficar1(resMMHC)
graficar2(resMMPC)
#Obtener BIC
scMMPC <- score(resMMPC,bn_df)
print(scMMPC)
#Propagación de la evidencia: EM (Ejecuación de Máxima Exprectación)
fittedbn.MMPC <- bn.fit(resMMPC, data = bn_df)



print(fittedbn.HC$E)
print(fittedbn.MMHC$E)
print(fittedbn.MMPC$E)

par(mfrow = c(1, 2))
graphviz.compare(dag, resHC, layout = "fdp" ,shape = "ellipse", main = c("DAG original", "DAG propio"))
graphviz.compare(fittedbn.HC, resHC, layout = "fdp" ,shape = "ellipse", main = c("DAG fitted", "DAG propio"))

#Obtención de tabla de probabilidades condicionales mediante EM


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