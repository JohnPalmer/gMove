\name{n.locs}

\alias{n.locs}
\alias{n.locs,SpatialPointsDataFrame-method}
\alias{n.locs,.MoveTrackStack-method}

\docType{methods}

\title{Extract the number of locations of a Move or MoveStack object}

\description{The n.locs method returns the number of locations of a track from a Move or MovesStack object.}

\usage{
  \S4method{n.locs}{SpatialPointsDataFrame}(obj)
  \S4method{n.locs}{.MoveTrackStack}(obj)
}

\arguments{
  \item{obj}{Move or MoveStack object}
}

\author{Marco Smolla}

\examples{
load(system.file("extdata", "move.RData", package="move"), .GlobalEnv)
stack <- moveStack(list(leroy,leroy))

n.locs(ricky) #number of locations of Move object
n.locs(stack) #number of locations of MoveStack object
}
