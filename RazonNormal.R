install.packages("jsonlite")
install.packages("openxlsx", dependencies = TRUE)
library(jsonlite)
library(sf)
library(readxl)


# Lectura de Json ---------------------------------------------------------
EH <- fromJSON("D:/CModelamientoSwatR/Data/roy/PE_senamhi_historic.json")
class(EH)
df  <-  st_as_sf(EH, coords = c('lon', 'lat'), crs = 4326)
st_write(df, "D:/CModelamientoSwatR/Data/roy/EH.shp")


# Lectura de hojas de calculo ---------------------------------------------
setwd("D:/CModelamientoSwatR/Data/EstacionesData")
sanm <- read.csv("sanmarcos.csv")
shug <- read.csv("shugar.csv")
chug <- read.csv("chugur.csv")
webe <- read.csv("weberbauer.csv")
sanmarcos <- sanm %>% gather(DATE, PP, "ï..1965":"X2014")
shugar <- shug %>% gather(DATE, PP, "ï..1965":"X2014")
chugur <- chug %>% gather(DATE, PP, "ï..1965":"X2014")
weberbauer <- webe %>% gather(DATE, PP, "ï..1965":"X2014")
Date <- seq(as.Date("1965/1/1"), as.Date("2014/12/1"), "months")
Estaciones <- data.frame(Date,Sanmarcos = sanmarcos$PP, Shugar = shugar$PP, Chugur = chugur$PP, Weberbauer = weberbauer$PP)

attach(Estaciones)

vect = c()
COMP <- function(a,b,c,an,bn,cn,n) {
  vect = c()
  for (i in 1:600) {
    if(sum(c(a[i],b[i],c[i]) %in% NA) == 1) {
      p1 = (n/an)*a[i]
      p2 = (n/bn)*b[i]
      p3 = (n/cn)*c[i]
      vect[i] = 1/2*sum(c(p1,p2,p3), na.rm = T)
    } else if(sum(c(a[i],b[i],c[i]) %in% NA) == 2) {
      p1 = (n/an)*a[i]
      p2 = (n/bn)*b[i]
      p3 = (n/cn)*c[i]
      vect[i] = sum(c(p1,p2,p3), na.rm = T)
    } else {
      p1 = (n/an)*a[i]
      p2 = (n/bn)*b[i]
      p3 = (n/cn)*c[i]
      vect[i] = 1/3*sum(c(p1,p2,p3))
    }
    # if(sum(c(a[i],b[i],c[i]) %in% NA) == 2) {
    #   p1 = (47533.32/n)*a[i]
    #   p2 = (6301.96/n)*b[i]
    #   p3 = (2702.78/n)*c[i]
    #   vect[i] = 1*sum(c(p1,p2,p3), na.rm = T)
    # }
  }
  return(vect)
}
# Sa = 3098.67
# Sh = 4009.96
# Ch = 6301.96
# We = 2702.78
Estaciones["Sa"] <- COMP(Shugar,Chugur,Weberbauer,4009.96,6301.96,2702.78,3098.67)
Estaciones["Sh"] <- COMP(Sanmarcos, Chugur, Weberbauer,3098.67,6301.96,2702.78,4009.96)
Estaciones["Ch"]<- COMP(Sanmarcos, Shugar, Weberbauer,3098.67,4009.96,2702.78,6301.96)
Estaciones["We"] <- COMP(Sanmarcos, Shugar, Chugur,3098.67,4009.96,6301.96,2702.78)

write.csv(Estaciones, "EstacionesComp.csv")


library(ggplot2)


fecha<-as.data.frame(as.POSIXct(datos.orig$Fecha))

# # Build a new and shorter data frame for temperature
# datos.graf<-cbind.data.frame(fecha,datos.orig$Temp_Max,datos.orig$Temp_Min)
# colnames(datos.graf)<-c("fecha","TMax","Tmin")

# Plot
ggplot() + 
  geom_line(data=datos.graf,aes(x=fecha, y=TMax),colour="red") +
  geom_line(data=datos.graf,aes(x=fecha, y=Tmin),colour="blue") +
  ylab("Temperature axis title") + 
  xlab(" ") +
  scale_x_datetime(
    expand=c(0,0),                           # avoid blank space around plot
    date_breaks = "1 month",                 # main axis breaks
    date_labels="%d/%m/%Y") +                # axis labels date format
  scale_y_continuous(expand=c(0,0), limits=c(-10,40)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle("Temperature main title") + 
  theme(plot.title = element_text(size=10, vjust=0.5, hjust=0.5))
