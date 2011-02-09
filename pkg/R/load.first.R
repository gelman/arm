.onAttach <- function(...) {
  mylib <- dirname(system.file(package = "arm"))
  ver <- packageDescription("arm", lib.loc = mylib)$Version
  builddate <- packageDescription("arm", lib.loc = mylib)$Date
  cat(paste("\narm (Version ", ver, ", built: ", builddate, ")\n", sep = ""))
  cat ("Working directory is", getwd(), "\n")
#  options(digits = 2, papersize="letter")
#  cat ("options( digits = ", getOption("digits"), ")\n")
  if(!any(search()=="package:car"))
    require(car)
  if(!any(search()=="package:foreign"))
    require(foreign) 
  if(!any(search()=="package:MASS"))
    require(MASS) 
  if(!any(search()=="package:Matrix"))
    require(Matrix)
  if(!any(search()=="package:lme4"))
    require(lme4) 
  if(!any(search()=="package:R2WinBUGS"))
    require(R2WinBUGS)
}
