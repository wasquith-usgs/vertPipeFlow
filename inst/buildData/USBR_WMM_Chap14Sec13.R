D2  <- read.table("2IN-ID.dat",  header=TRUE);
D3  <- read.table("3IN-ID.dat",  header=TRUE);
D4  <- read.table("4IN-ID.dat",  header=TRUE);
D5  <- read.table("5IN-ID.dat",  header=TRUE);
D6  <- read.table("6IN-ID.dat",  header=TRUE);
D7  <- read.table("7IN-ID.dat",  header=TRUE);
D8  <- read.table("8IN-ID.dat",  header=TRUE);
D10 <- read.table("10IN-ID.dat", header=TRUE);
D12 <- read.table("12IN-ID.dat", header=TRUE);

lgH <- seq(log10(1.5), log10(60), by=0.01);
lgQ <- seq(log10(30), log10(3000), by=0.01);

H2  <- approx(log10(D2$FLOWgpm),  log10(D2$Hinches),  lgQ)$y;
H3  <- approx(log10(D3$FLOWgpm),  log10(D3$Hinches),  lgQ)$y;
H4  <- approx(log10(D4$FLOWgpm),  log10(D4$Hinches),  lgQ)$y;
H5  <- approx(log10(D5$FLOWgpm),  log10(D5$Hinches),  lgQ)$y;
H6  <- approx(log10(D6$FLOWgpm),  log10(D6$Hinches),  lgQ)$y;
H7  <- approx(log10(D7$FLOWgpm),  log10(D7$Hinches),  lgQ)$y;
H8  <- approx(log10(D8$FLOWgpm),  log10(D8$Hinches),  lgQ)$y;
H10 <- approx(log10(D10$FLOWgpm), log10(D10$Hinches), lgQ)$y;
H12 <- approx(log10(D12$FLOWgpm), log10(D12$Hinches), lgQ)$y;

.USBRfig14_12GIF <-
        data.frame(lgQ=lgQ,
                   lgH2=H2, lgH3=H3, lgH4=H4, lgH5=H5,
                   lgH6=H6, lgH7=H7, lgH8=H8, lgH10=H10, 
                   lgH12=H12);
.USBRfig14_12ID <- c(2,3,4,5,6,7,8,10,12);

save(.USBRfig14_12GIF, .USBRfig14_12ID, 
     file="sysdata.rda")


plot(.USBRfig14_12GIF$lgQ, .USBRfig14_12GIF$lgH2,
     xlim=c(log10(30), log10(3000)),
     ylim=c(log10(1.5), log10(60)), type="l")
lines(.USBRfig14_12GIF$lgQ, .USBRfig14_12GIF$lgH3)
lines(.USBRfig14_12GIF$lgQ, .USBRfig14_12GIF$lgH4)
lines(.USBRfig14_12GIF$lgQ, .USBRfig14_12GIF$lgH5)
lines(.USBRfig14_12GIF$lgQ, .USBRfig14_12GIF$lgH6)
lines(.USBRfig14_12GIF$lgQ, .USBRfig14_12GIF$lgH7)
lines(.USBRfig14_12GIF$lgQ, .USBRfig14_12GIF$lgH8)
lines(.USBRfig14_12GIF$lgQ, .USBRfig14_12GIF$lgH10)
lines(.USBRfig14_12GIF$lgQ, .USBRfig14_12GIF$lgH12)


#www.usbr.gov/pmts/hydraulics_lab/pubs/wmm/chap14_13.html
