\name{array.append}
\alias{array.append}

\title{Array Merging}
\description{
 Add arrays to an array
}

\usage{
array.append(array1, array2, d)
}

\arguments{
  \item{array1}{The array to be appended}.
  \item{array2}{The array to be included in the appended array}.
  \item{d}{append after which dimension}
}

value{
array
}

\author{Masanao Yajima \email{yajima@stat.columbia.edu},
        Yu-Sung Su \email{ys463@columbia.edu}
}


\examples{
a <- array(1:27, c(3,3,3))
b <- array(28:54, c(3,3,3))
array.append(array1=a, array2=b, d=3)
}
\keyword{manip}
