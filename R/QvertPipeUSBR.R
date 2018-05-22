"QvertPipeUSBR" <- 
function(h=NULL, d=NULL, digits=3,
         hunits=c("inch", "feet"),
         dunits=c("inch", "feet"),
         Qunits=c("gpm",  "cfs"),
         forceTable=FALSE) {

  hunits <- match.arg(hunits);
  dunits <- match.arg(dunits);
  Qunits <- match.arg(Qunits);
 
  the.h <- h; the.d <- d;
  if(hunits == "feet") h <- h*12;
  if(dunits == "feet") d <- d*12;
  
  if(length(h) != 1) {
    warning("h can only be a scalar");
    return(NA);
  }
  if(length(d) != 1) {
    warning("d can only be a scalar");
    return(NA);
  }

  z <- list(h=the.h, d=the.d, Q=NA,
            message=NA,
            equation=NA,
            error=NA,
            hunits=hunits,
            dunits=dunits,
            Qunits=Qunits);

  if(forceTable == FALSE) {          
    if(h < 0.37*d) {
      Q <- 6.17 * d^1.25 * h^1.35;
      if(Qunits == "cfs") Q <- Q*0.13368/60;
      z$message <- "circular-weir flow";
      z$equation <- "Q=6.17*d^1.25*h^1.35";
      z$Q <- round(Q, digits=digits);
      return(z);	
    } else if(h > 1.4*d) {
      Q <- 5.01 * d^1.99 * h^0.53;
      if(Qunits == "cfs") Q <- Q*0.13368/60;
      z$equation <- "Q=5.01*d^1.99*h^0.53";
      z$message <- "jet flow";
      z$Q <- round(Q, digits=digits);
      return(z);
    }
    z$message <- "quasi-jet flow";
  } else {
    z$message <- "forcing quasi-jet flow";
  }
  
  z$error <- "(h [1.5,60] or d[2,12] outside of table)";
  if(d < 2) {
    warning("d is < 2 inches");
    return(z);
  }
  if(d > 12) {
    warning("d is > 12 inches");
    return(z);
  }
  if(h < 1.5) {
    warning("h is < 1.5 inches");
    return(z);
  }
  if(h > 60) {
    warning("h is > 60 inches");
    return(z);
  }
  
  z$error <- NA;
  z$equation <- "table interp.";
  
  nID   <- seq(1:length(.USBRfig14_12ID));
  
  IDlow <- max(.USBRfig14_12ID[.USBRfig14_12ID <= d]);
  IDhi  <- min(.USBRfig14_12ID[.USBRfig14_12ID >= d]);

  if(is.na(IDlow)) {
    txt <- paste(c("variable IDlow is NA for h=", h, " and d=", d), collapse="");
    warning(txt);
    z$error <- txt;
    return(z);
  }
  if(is.na(IDhi)) {
    txt <- paste(c("variable IDhi is NA for h=", h, " and d=", d), collapse="");
    warning(txt);
    z$error <- txt;
    return(z);
  }

  idxIDlow <- nID[.USBRfig14_12ID == IDlow];
  idxIDhi  <- nID[.USBRfig14_12ID == IDhi ];

  lgQ    <- .USBRfig14_12GIF$lgQ;
  lgHlow <- .USBRfig14_12GIF[, (idxIDlow+1)];
  lgHhi  <- .USBRfig14_12GIF[, (idxIDhi+1) ];
  lgQlow <- approx(lgHlow, lgQ, log10(h), ties="ordered")$y;
  lgQhi  <- approx(lgHhi,  lgQ, log10(h), ties="ordered")$y;
  
  if(is.na(lgQlow)) {
    txt <- paste(c("variable lqQlow is NA for h=", h," and d=", d), collapse="");
    warning(txt);
    z$error <- txt;
    return(z);
  }
  if(is.na(lgQhi)) {
    txt <- paste(c("variable lqQhi is NA for h=", h, " and d=", d), collapse="");
    warning(txt);
    z$error <- txt;
    return(z);
  }
  
  xs <- log10(c(IDlow,IDhi));
  ys <- c(lgQlow,lgQhi);
  my.lgQ <- approx(xs, ys, log10(d), ties="ordered")$y;
  Q <- 10^my.lgQ;
  if(Qunits == "cfs") Q <- Q*0.13368/60;
  z$Q <- round(Q, digits=digits);
  return(z);
}

