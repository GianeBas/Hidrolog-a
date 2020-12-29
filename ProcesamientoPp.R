library(raster)
library(sf)
library(dplyr)
library(jsonlite)
library(sf)
library(readxl)
library(tidyverse)
library(climatol)
library(ncdf4)

# Lectura de hojas de calculo ---------------------------------------------
setwd("D:/VI CICLO/hidro/Proyecto/datos/Precipitacion")
tabla <- read.csv("Procesamiento/Crudo/Estaciones1988_2013.csv")
tabla2 <- read.csv("Procesamiento/Crudo/Metadata.csv", stringsAsFactors = F)
# Funcion limpiar -99.9 ----------------------------------------------------------------

Limp <- function (data) {
  for (j in 1:length(data)) {
    vec = c()
    for (i in 1:nrow(data)) {
      if (is.na(data[i,j])) {
        data[i,j] = NA
      } else if (data[i,j] == -99.9) {
        data[i,j] = NA
      } else {
        data[i,j] = data[i,j]
      }
    }
  }
  return(data)
}

data <- Limp(tabla)

write(as.matrix(data), "Procesamiento/InsumoClimatol/Pp_1988-2013.dat")
write.table(tabla2, "Procesamiento/InsumoClimatol/Pp_1988-2013.est", row.names = F, col.names = F)


# Climatol ----------------------------------------------------------------

setwd("D:/VI CICLO/hidro/Proyecto/datos/Precipitacion/Procesamiento/InsumoClimatol")
homogen("Pp", 1988, 2013, expl = T)
dahstat("Pp", 1988, 2013, stat = 'series')



# RegRAIN -----------------------------------------------------------------

RegRAIN <- function(datos,dem,ini,fin,crossv){
  
  if (class(datos$FECHA) !="Date") {#Identificar si la columna fecha en el dataframe con los datos se encuentra en formato "Date"
    
    print(paste("Error: The column FECHA in your dataframe is not of the class Date and in the format %y-%m-%d", sep="")) #Imprime mensaje de error
    
  } else { #Si se encuentra en formato "Date" sigue el proceso
    
    options(digits=20) #Convertir las columnas numericas con 20 digitos
    datos$LONGITUD <- as.numeric(as.character(datos$LONGITUD)) #Convertir la columna de "LONGITUD" al tipo "numeric"
    datos$LATITUD <- as.numeric(as.character(datos$LATITUD)) #Convertir la columna de "LATITUD" al tipo "numeric"
    fechas = unique(datos$FECHA) #Indicar que cada fecha en la columna "FECHA" es un identificador unico
    fechas = sort(fechas) #Ordenar las fechas
    
    #Obtener dataframe con estaciones totales ingresadas
    DATA1 <- datos[!duplicated(datos$CODIGO), ] #Obtener el numero de estaciones a procesar
    coordinates(DATA1) = ~LONGITUD + LATITUD #Convetir a objeto de la clase "SpatialPointsDataFrame"
    
    #DERIVAR ASPECTO DEL DEM
    ASPECTO <- terrain(dem, opt=c('aspect'), unit='degrees', neighbors=8)
    #DERIVAR PENDIENTE DEL DEM
    PENDIENTE <- terrain(dem, opt=c('slope'), unit='degrees', neighbors=8)
    #EXTRAER INFORMACION DE ASPECTO DEL RASTER DE ASPECTO GENERADO
    aspecto <- data.frame(coordinates(DATA1), DATA1$CODIGO, extract(ASPECTO, DATA1)) #Extraer los datos en un data frame
    names(aspecto)[1:4] <- c("X", "Y","CODIGO", "ASP") #Nombrar las columnas
    aspecto <- data.frame(aspecto$CODIGO, aspecto$ASP) #Seleccionar columnas de interes
    names(aspecto)[1:2] <- c("CODIGO", "ASP") #Nombrar las columnas
    aspecto<-aspecto[complete.cases(aspecto),] #Suprimir filas con valores NA del archivo datos
    datos = merge(datos,aspecto,by="CODIGO") #Unir el data frame de extraccion con el archivo de datos totales
    
    #EXTRAER INFORMACION DE ELEVACION DEL RASTER DE ASPECTO GENERADO
    dem1 <- data.frame(coordinates(DATA1), DATA1$CODIGO, extract(dem, DATA1)) #Extraer los datos en un data frame
    names(dem1)[1:4] <- c("X", "Y","CODIGO", "Z") #Nombrar las columnas
    dem1 <- data.frame(dem1$CODIGO, dem1$Z) #Seleccionar columnas de interes
    names(dem1)[1:2] <- c("CODIGO", "DEM") #Nombrar las columnas
    dem1<-dem1[complete.cases(dem1),] #Suprimir filas con valores NA del archivo datos
    datos = merge(datos,dem1,by="CODIGO") #Unir el data frame de extraccion con el archivo de datos totales
    
    #EXTRAER INFORMACION DE PENDIENTES DEL RASTER DE PENDIENTE GENERADO
    pendiente <- data.frame(coordinates(DATA1), DATA1$CODIGO, extract(PENDIENTE, DATA1)) #Extraer los datos en un data frame
    names(pendiente)[1:4] <- c("X", "Y","CODIGO", "PEN") #Nombrar las columnas
    pendiente <- data.frame(pendiente$CODIGO, pendiente$PEN) #Seleccionar columnas de interes
    names(pendiente)[1:2] <- c("CODIGO", "PEN") #Nombrar las columnas
    pendiente<-pendiente[complete.cases(pendiente),] #Suprimir filas con valores NA del archivo datos
    datos = merge(datos,pendiente,by="CODIGO") #Unir el data frame de extraccion con el archivo de datos totales
    
    DATA1<- datos #Crear un data frame de respaldo al de datos totales con las columnas extraidas
    
    p <- raster(extent(bbox(dem))) #Crea un objeto de tipo raster con extent igual al DEM de entrada
    res(p) <- res(dem) #Asigna una resolucion de pixel acorde con el DEM de entrada
    p <- resample(p, dem) #Iguala la estructura y caracteristicas a las del DEM de entrada
    res <- data.frame(res(p)) #Extrae en un data frame la resolucion del DEM de entrada
    res1<-as.numeric(as.character(res[1:1,])) #Extrae como un vector la resolucion en X del DEM
    res2<-as.numeric(as.character(res[2:2,])) #Extrae como un vector la resolucion en Y del DEM
    xy <- data.frame(datos$LATITUD,datos$LONGITUD) #Crea un dataframe con las coordenadas de los datos de entrada
    xy <- unique(xy) #Asigna valores unicos a las coordenadas
    suppressWarnings(spline1 <- Tps(xy,xy$datos.LONGITUD)) #Suprime advertencias y mensajes del proceso de interpolacion
    spline.lat <- interpolate(p,spline1) #Asigna los valores interpolados al raster "p" creado inicialmente en un nuevo raster
    
    suppressWarnings(spline2 <- Tps(xy,xy$datos.LATITUD)) #Suprime advertencias y mensajes del proceso de interpolacion
    spline.long <- interpolate(p,spline2) #Asigna los valores interpolados al raster "p" creado inicialmente en un nuevo raster
    
    #Crear Grid para la interpolacion IDW
    project <- proj4string(dem) #define un objeto de proyeccion geografica a partir del DEM de entrada
    xmin <- xmin(dem) #Extrae xmin del DEM de entrada
    xmax <- xmax(dem) #Extrae xmax del DEM de entrada
    ymin <- ymin(dem) #Exrae ymin del DEM de entrada
    ymax <- ymax(dem) #Extrae ymax del DEM de entrada
    x.range <- as.numeric(c(xmin,xmax)) #Genera un rango de X a partir de xmin y xmax
    y.range <- as.numeric(c(ymin,ymax)) #Genera un rango de Y a partir de ymin y ymax
    grd <- expand.grid(x = seq(from = x.range[1], to = x.range[2], by = res1), y = seq(from = y.range[1], to = y.range[2], by = res2)) #Expande la grilla a partir de los valores minimos y maximos de X - Y
    coordinates(grd) <- ~x + y #Definir coordenadas espaciales para crear un objeto georreferenciado
    gridded(grd) <- TRUE #Asignar tipo grilla
    crs(grd) <- project #Proyectar la grilla
    
    for (i in ini:fin)  { #Bucle de proceso para la interpolacion RegRAIN
      PPT = subset(datos,datos$FECHA==fechas[i]) #Extraccion de la base de datos total para una fecha especifica
      print(paste("Starting RegRAIN Interpolation",fechas[i], sep=" ")) #Imprime mensaje en la consola
      PPT[PPT == 99999] <- NA #Reemplaza todos los valores 99999 por NAs
      PPT = as.data.frame(na.omit(PPT)) #Elimina NAs del analisis
      PPT$FECHA <- as.Date(PPT$FECHA) #Convierte a Date la columna FECHA
      RLM = lm(PPT ~ ASP + DEM + PEN + LONGITUD + LATITUD, data=PPT) #Aplica la Regresion Lineal Multiple
      
      cons <- summary(RLM)$coefficients[1,1] #Extrae la constante de la regresion como vector
      asp  <- summary(RLM)$coefficients[2,1] #Extrae el coeficiente de aspecto como vector
      demc  <- summary(RLM)$coefficients[3,1] #Extrae el coeficiente de elevacion como vector
      pen  <- summary(RLM)$coefficients[4,1] #Extrae el coeficiente de pendiente como vector
      long <- summary(RLM)$coefficients[5,1] #Extrae el coeficiente de longitud como vector
      lat <- summary(RLM)$coefficients[6,1] #Extrae el coeficiente de latitud como vector
      
      #Calcula el residuo de la regresion lineal multiple
      residuo = PPT$PPT - (asp*PPT$ASP) - (demc*PPT$DEM) - (pen*PPT$PEN) - (long*PPT$LONGITUD) - (lat*PPT$LATITUD)
      residuo2 <- as.data.frame(residuo) #convierte el residuo en data frame
      residuo3 <-as.data.frame(na.omit(residuo2)) #Omite NAs
      d <- as.data.frame(row.names(residuo3)) #Prepara el data frame para interpolacion
      residuo4 <- data.frame(residuo3, d) #Prepara el data.frame para interpolacion
      names(residuo4)[1:2] <- c("residuo","ID") #Asigna nombres a columnas
      f <- data.frame(PPT$LONGITUD, PPT$LATITUD, residuo2) #Prepara el data frame para interpolacion
      names(f)[1:3] <- c("X","Y", "residuo") #Asigna nombres a columnas
      coordinates(f) <- ~X + Y #Definir coordenadas espaciales para crear un objeto georreferenciado
      crs(f) <- project #Proyecta el data frame de residuos
      z= capture.output(idw.res <- idw(formula = residuo ~ 1, locations = f, newdata = grd)) # Aplica el interpolador IDW al residuo eliminando mensaje IDW
      raster.idw <- raster(idw.res) #Convierte la superficie IDW a raster
      Residuo <- resample(raster.idw, dem) #Asigna caracteristicas del DEM al raster de IDW generado
      
      RegRAIN_PPT <- (ASPECTO*asp) + (dem*demc) + (PENDIENTE*pen) + (spline.lat*lat) + (spline.long*long) + Residuo #Realiza la interpolacion RegRAIN
      RegRAIN_PPT[RegRAIN_PPT<0] <- 0 #Asigna valores menores a 0 = 0
      
      if (crossv == FALSE) { #Condicional para Cross Validation = FALSE
        
        print(paste("Exporting RegRAIN Raster to working directory", sep="")) #Imprime mensaje en la consola
        writeRaster(RegRAIN_PPT, paste("RegRAIN_PPT_",fechas[i],".tif",sep=""), drivername = "GTiff", overwrite=TRUE) #Exporta raster al directorio de trabajo
        
      } else {
        
        print(paste("Exporting RegRAIN Raster to working directory",sep="")) #Imprime mensaje en la consola
        writeRaster(RegRAIN_PPT, paste("RegRAIN_PPT_",fechas[i],".tif",sep=""), drivername = "GTiff", overwrite=TRUE) #Exporta raster al directorio de trabajo
        
        #Crear data frame para datos de validacion cruzada
        PPT_cv <- PPT
        PPT_cv$PPT_cross <- NA
        PPT_cv[] <- NA
        PPT_cv$FECHA <- as.Date(PPT_cv$FECHA)
        #Crear data frame para datos de validacion cruzada
        codigo = unique(PPT$CODIGO) #Indicar que cada campo en la columna CODIGO es un identificador unico
        codigo = sort(codigo) #Ordenar los codigos
        
        print(paste("Starting Cross Validation Process ",fechas[i],sep="")) #Imprime mensaje en la consola
        
        for (i in 1:nrow(PPT))  {
          CODIGO=NULL #Asigna un identificador a la variable CODIGO
          PPT_cross <- subset(PPT, !(CODIGO %in% codigo[i])) #Extraccion de la base de datos de la fecha interpolada para un codigo especifico
          PPT_cross1 <- subset(PPT, (CODIGO %in% codigo[i])) #Extraccion de la base de datos de la informacion asociada al codigo especifico
          coordinates(PPT_cross1) = ~LONGITUD + LATITUD #Convetir a objeto de la clase "SpatialPointsDataFrame"
          RLM = lm(PPT ~ ASP + DEM + PEN + LONGITUD + LATITUD, data=PPT_cross) #Aplica la Regresion Lineal Multiple
          
          cons <- summary(RLM)$coefficients[1,1] #Extrae la constante de la regresion como vector
          asp  <- summary(RLM)$coefficients[2,1] #Extrae el coeficiente de aspecto como vector
          demc  <- summary(RLM)$coefficients[3,1] #Extrae el coeficiente de elevacion como vector
          pen  <- summary(RLM)$coefficients[4,1] #Extrae el coeficiente de pendiente como vector
          long <- summary(RLM)$coefficients[5,1] #Extrae el coeficiente de longitud como vector
          lat <- summary(RLM)$coefficients[6,1] #Extrae el coeficiente de latitud como vector
          
          #Calcula el residuo de la regresion lineal multiple
          residuo = PPT_cross$PPT - (asp*PPT_cross$ASP) - (demc*PPT_cross$DEM) - (pen*PPT_cross$PEN) - (long*PPT_cross$LONGITUD) - (lat*PPT_cross$LATITUD)
          residuo2 <- as.data.frame(residuo) #Convierte el residuo en data.frame
          residuo3 <-as.data.frame(na.omit(residuo2)) #Omite NAs
          d <- as.data.frame(row.names(residuo3)) #Prepara el data frame para interpolacion
          residuo4 <- data.frame(residuo3, d) #Prepara el data frame para interpolacion
          names(residuo4)[1:2] <- c("residuo","ID") #Asigna nombres a columnas
          f <- data.frame(PPT_cross$LONGITUD, PPT_cross$LATITUD, residuo2) #Prepara el data frame para interpolacion
          names(f)[1:3] <- c("X","Y", "residuo") #Asigna nombres a columnas
          coordinates(f) <- ~X + Y #Definir coordenadas espaciales para crear un objeto georreferenciado
          crs(f) <- project #Proyecta el data frame de residuos
          j= capture.output(idw.res <- idw(formula = residuo ~ 1, locations = f, newdata = grd)) # Aplica el interpolador IDW al residuo eliminando mensaje IDW
          raster.idw <- raster(idw.res) #Convierte la superficie IDW a raster
          Residuo <- resample(raster.idw, dem) #Asigna caracteristicas del DEM al raster de IDW generado
          
          RegRAIN_PPT_cross <- (ASPECTO*asp) + (dem*demc) + (PENDIENTE*pen) + (spline.lat*lat) + (spline.long*long) + Residuo #Realiza la interpolacion RegRAIN
          RegRAIN_PPT_cross[RegRAIN_PPT_cross<0] <- 0 #Asigna valores menores a 0 = 0
          
          PPT_test <- data.frame(coordinates(PPT_cross1), PPT_cross1$CODIGO, extract(RegRAIN_PPT_cross, PPT_cross1)) #Extraer los datos en un data frame
          names(PPT_test)[1:4] <- c("x", "Y","CODIGO", "PPT_cross") #Nombrar las columnas
          PPT_test <- data.frame(PPT_test$CODIGO, PPT_test$PPT_cross) #Seleccionar columnas de interes
          names(PPT_test)[1:2] <- c("CODIGO", "PPT_cross") #Nombrar las columnas
          PPT_test<-PPT_test[complete.cases(PPT_test),]##Suprimir filas con valores NA del archivo datos
          PPT_test1 = merge(PPT,PPT_test,by="CODIGO") #Unir el data frame de extraccion con el archivo de datos totales
          PPT_cv[i,] <- PPT_test1 #Selecciona la fila i del datanfrmae PPT_cv
          fecha <- PPT_cv$FECHA[1] #Extrae la primera fecha del data frame
          avance <- (i/nrow(PPT))*100 #Calcula el avance de la validacion cruzada
          print(paste("Running Cross Validation ",round(avance, digits = 1),"%", sep=" ")) #Imprime mensaje en la consola
          
        }
        fecha <- PPT_cross$FECHA[1] #Extrae la primera fecha del data.frame
        if (sum(PPT_cv$PPT) > 0) { #Condicional que evalua que la suma de la columna PPT sea mayor a "0"
          png(filename= paste("RegRAIN_Cross_Validation_Plot_",fecha,".png",sep=""), #Abre el comando de generacion del png
              units="cm", #Define unidades
              width=10, #Define el ancho del grafico
              height=10, #Define el alto del grafico
              pointsize=6, #Define unidades de tamano de puntos en el grafico
              res=250) #Define resolucion del grafico
          diagram <- plot(PPT_cv$PPT_cross, PPT_cv$PPT, xlab = "Predicted Precipitation (mm)",
                          ylab = "Observed Precipitation(mm)", main= paste("RegRAIN Cross Validation Plot", fecha, sep=" ")) #Crea Plot de Validacion Cruzada
          mod <- nls(PPT_cv$PPT_cross ~ exp(a + b * PPT_cv$PPT), data=PPT_cv, start = list(a = 0, b = 0) ) #Aplica modelo exponencial al grafico
          lines(PPT_cv$PPT, predict(mod, list(x = PPT_cv$PPT),col="red")) # Adiciona curva ajustada al modelo exponencial
          RSS.p<-sum(residuals(mod)^2) #Suma residual de cuadrados
          TSS<-sum((PPT_cv$PPT-mean(PPT_cv$PPT_cross))^2) #Suma total de cuadrados
          r.squared<-1-(RSS.p/TSS) #Calcula coeficiente de determinacion R2
          #Genera leyenda al grafico
          legend("topright", bty="n", legend=c(paste("R2 is", format(r.squared, digits=4)),paste("RMSE is", format(rmse(PPT_cv$PPT_cross, PPT_cv$PPT), digits=4)),paste("MAE is", format(mae(PPT_cv$PPT_cross, PPT_cv$PPT, na.rm=TRUE), digits=4))))
          print(paste("Exporting Cross Validation Plot ",fecha,sep=" ")) #Imprime mensaje en la consola
          dev.off() #Cierra el comando de generacion del png
          print(paste("End of Cross Validation Process ",fecha,sep=" ")) #Imprime mensaje en la consola
        } else {
          png(filename= paste("RegRAIN_Cross_Validation_Plot_",fecha,".png",sep=""), #Abre el comando de generacion del png
              units="cm", #Define unidades
              width=10, #Define el ancho del grafico
              height=10, #Define el alto del grafico
              pointsize=6, #Define unidades de tamano de puntos en el grafico
              res=250) #Define resolucion del grafico
          diagram <- plot(PPT_cv$PPT_cross, PPT_cv$PPT, xlab = "Predicted Precipitation (mm)",
                          ylab = "Observed Precipitation(mm)", main= paste("RegRAIN Cross Validation Plot", fecha, sep=" ")) #Crea Plot de Validacion Cruzada
          #Genera leyenda al grafico
          legend("topright", bty="n", legend=c(paste("RMSE is", format(rmse(PPT_cv$PPT_cross, PPT_cv$PPT), digits=4)),paste("MAE is", format(mae(PPT_cv$PPT_cross, PPT_cv$PPT, na.rm=TRUE), digits=4))))
          dev.off() #Cierra la generacion del png
          print(paste("Exporting Cross Validation Plot ",fecha,sep=" ")) #Imprime mensaje en la consola
          print(paste("End of Cross Validation Process ",fecha,sep=" ")) #Imprime mensaje en la consola
        }
        
      }
      
    }
    
  }
  
}
setwd("D:/VI CICLO/hidro/Proyecto/datos/Precipitacion/Procesamiento/RegRAIN/Resultados")
RegRAIN(datos_prec, dem, 4384, 9496, crossv = F) #Correr la función.

# Aplicación por centroides -----------------------------------------------

subs <- st_read("D:/VI CICLO/hidro/Proyecto/datos/Precipitacion/Procesamiento/RegRAIN/Automatic/Subs.shp")
lista <- list.files()
centros <- st_centroid(subs)[1]
# st_write(centros,"centroidslvl8.shp")
setwd("D:/VI CICLO/hidro/Proyecto/datos/Precipitacion/Procesamiento/RegRAIN/Automatic/Raster")
a <- data.frame(matrix(ncol = 5))
names(a) <- c("C49897685","C49897686","C49897687","C49897688","C49897689")
for (i in 1:length(lista)) {
  extract <- raster::extract(raster(lista[i]), centros, df = T)
  for (j in 1:5){
    a[i,j] <- extract[j,2]
  } 
  print(i)
}

a["Fecha"] <- seq(as.Date("2000-01-01"), as.Date("2009-12-31"), 1)
write.csv(a, "D:/VI CICLO/hidro/Proyecto/datos/Precipitacion/Procesamiento/RegRAIN/ResultadosTbl/PpDiaria.csv")



# Estadisticas de zona ----------------------------------------------------

setwd("D:/VI CICLO/hidro/Proyecto/datos/Precipitacion/Procesamiento/RegRAIN/Automatic/Raster")
subs <- st_read("D:/VI CICLO/hidro/Proyecto/datos/Precipitacion/Procesamiento/RegRAIN/Automatic/Subs.shp")
lista <- list.files()
a <- data.frame(matrix(ncol = 5))
names(a) <- c("C49897685","C49897686","C49897687","C49897688","C49897689")
for (i in 1:length(lista)) {
  extract <- raster::extract(raster(lista[i]), subs, df = T)
  names(extract)[2] <- "Pp"
  Pprom <- group_by(extract, ID) %>%
    summarise(Pp = mean(Pp))
  for (j in 1:5){
    a[i,j] <- Pprom[j,2]
  } 
}
a["Fecha"] <- seq(as.Date("2000-01-01"), as.Date("2009-12-31"), 1)

