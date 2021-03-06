\name{angle}

\alias{angle}
\alias{angleSummary}
\alias{angle,.MoveTrackSingle-method}
\alias{angle,.MoveTrackStack-method}
\alias{angleSummary,.MoveTrackSingle-method}
\alias{angleSummary,.MoveTrackStack-method}

\docType{methods}

\title{angle information from a track or track stack}

\description{This function returns a summary about angle related measurements of a track or track stack. These are: average azimuth, variance of azimuth, standard error of azimuth. }

\usage{
  \S4method{angle}{.MoveTrackSingle}(x)
  \S4method{angle}{.MoveTrackStack}(x)
  \S4method{angleSummary}{.MoveTrackSingle}(x)
  \S4method{angleSummary}{.MoveTrackStack}(x)
}

\arguments{
  \item{x}{Move or MoveStack object}
}

\author{Marco Smolla}

\value{
Angles in degrees
}

\examples{
  load(system.file("extdata", "move.RData", package="move"), .GlobalEnv)
\dontrun{  
  stack <- moveStack(list(leroy,leroy))
  angle(leroy) #angles from a Move object 
  angle(stack) #angles from a MoveStack object
  angleSummary(leroy) # summary of angle measures of a Move object
  angleSummary(stack) # summary of angle measures of a MoveStack object
}# failed on cran build
}
