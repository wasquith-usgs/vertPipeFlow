\name{QvertPipeUSBR}
\alias{QvertPipeUSBR}
\title{ Discharge from a Vertical Pipe by U.S. Bureau of Reclamation }
\description{
This function solve equations for circular-weir and jet flow issuing from a vertical pipe. For conditions of quasi-jet flow an interpolation into a table derived from a graphical image of discharge curves is used. The height of water above the pipe is defined as \eqn{h}, and the internal diameter of the pipe is defined as \eqn{d}. If \eqn{h < 0.37d} then the pipe is flowing as a circular weir and \eqn{Q = 6.17d^{1.25}h^{1.35}}. If \eqn{h > 1.4d} then the pipe is flowing as a jet and \eqn{Q = 5.01d^{1.99}h^{0.53}}. Otherwise the pipe is flowing as a quasi-jet and explicit equations are not available, but are derived from figure 14-12 in USBR (2001) that is attributed to Stock (1955).  Bos (1989) reports that the jet equation has an accuracy of \eqn{\pm 10} and the weir equation has an accuracy of \eqn{\pm 15}.
}
\usage{
QvertPipeUSBR(h=NULL, d=NULL, digits=3,
              hunits=c("inch", "feet"),
              dunits=c("inch", "feet"),
              Qunits=c("gpm",  "cfs"),
              forceTable=FALSE)}
\arguments{
  \item{h}{A scalar value for the height of water above the pipe,}
  \item{d}{A scalar value for the internal diameter of the pipe,}
  \item{digits}{The number of digits to round the discharge result,}
  \item{hunits}{The units on the \code{h},}
  \item{dunits}{The units on the \code{d},}
  \item{Qunits}{The units on the discharge, and}
  \item{forceTable}{A logical controlling whether the quasi-jet table is used instead of equations.}
}
\value{
  A list of results including the discharge.
}
\references{
Bos, M.G., 1989, Discharge Measurement Structures, ILRI Publication 20, Third Revised Edition, Wageningen, The Netherlands. 

Stock, E.M., 1955, Measuring of Irrigation Water, Bulletin No. 5, Utah State Engineering Experiment Station and Utah Cooperative Extension Service.

United States Bureau of Reclamation, 2001, Water Measurement Manual---A water resources technical publication: United States Department of the Interior Bureau of Reclamation, \url{www.usbr.gov/pmts/hydraulics_lab/pubs/wmm/chap14_13.html}
}
\author{ W.H. Asquith}
\examples{
str(QvertPipeUSBR(h=.15, d=8))
str(QvertPipeUSBR(h=3,   d=8))
str(QvertPipeUSBR(h=15,  d=8))

replicate(10,
            QvertPipeUSBR(h=runif(1, min=1.5, max=60),
                          d=runif(1, min=2,   max=12)));

H <- seq(.01,100,by=0.1)
Q <- sapply(H, function(x) {
  	        return(QvertPipeUSBR(h=x,d=8)$Q) } )
plot(H,Q, type="l")

}
\keyword{misc}
