meshgrid <- function(a,b) {
  ##=========================================
  # Meshgrid function which has the same    #
  # functionality as the one used in Matlab #
  ##=========================================
  list(
    x=outer(b*0,a,FUN="+"),
    y=outer(b,a*0,FUN="+")
  )
}