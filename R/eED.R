#' Calculate expected Evolutionary Distinctiveness (expected ED) for every tip of a phylogenetic tree
#'
#' \code{eED} calculates the expected Evolutionary Distinctiveness (expected ED) for
#' each tip of a phylogenetic tree given probabilities of extinction for each
#' tip, speciation and extinction rates, and a time into the future. It also
#' calculates several other values for each edge of the tree like probability of
#' extinction and number of subtending daughter tips.
#'
#'
#' @param tree An object of class phylo.
#'
#' @param tip.extinction.probabilities A named numeric vector of extinction (not survival) probabilities for each tip of the tree. Names must match the tip labels of the tree. If no probabilites are supplied, all tips are given a 0 probability of extinction and ED is calculated.
#'
#' @param lambda Numeric. A single instantaneous speciation rate in lineages per million species years.
#'
#' @param mu Numeric. A single instantaneous extinction rate in lineages per million species years.
#'
#' @param tMa Numeric. How many millions of years into the future will expected ED be calculated at? For normal function use where one is just looking at expected ED without any future evolution, tMa should be 0.
#'
#' @param auto.save Logical. Automatically write the output of the function to the working directory?
#'
#' @param source.of.data Character. Optional data tag to include in the function output.
#'
#' @section Background: Evolutionary Distinctiveness (ED) (Redding \emph{et
#'   al.}, 2014) fairly apportions Phylogenetic Diversity (PD) among tips
#'    of a phylogenetic tree and can be calculated in the
#'   function \code{\link[picante]{evol.distinct}}. From (Redding \emph{et al.}, 2014),
#'   ED "...for species i is the sum of edge lengths along the path from i to
#'   the root, each edge divided by the number of species ultimately subtending
#'   it". Thus, ED is the amount of unique evolutionary history that can be
#'   attributed to each tip. ED assumes that all tips have a probability of 1 of
#'   being measured, though.
#'
#'   Expected Evolutionary Distinctiveness (expected ED) is the probabilistic
#'   implementation of ED so that expected ED is the expected amount of unique
#'   evolutionary history that can be attributed to each tip. The sum of all
#'   tip's ED values in a tree will equal PD but the sum of all tip's expected ED values
#'   will equal expected PD (Faith, 2008). Each tip is given a probability of extinction from
#'   0 to 1 that could reflect the taxon's actual extinction probability (e.g. IUCN) or
#'   the probability of \emph{not} sampling this taxon in a certain community. Note that if each tip is given a 0 probability of extinction, expected ED simplifies to ED and this function will return the same results as \code{\link[picante]{evol.distinct}}, \code{type="fair.proportion"}.
#'
#' @section Normal Use:
#' Typical usage is
#'
#' \code{eED(tree, tip.extinction.probabilities)}
#'
#' This will calculate expected ED on a tree without any projected future evolution. Most users will only need to designate a tree and a vector of extinction (\strong{not survial!}) probabilties named with labels that match the tip labels in the tree. The time \code{tMa=0} is set by default.
#'
#' @section Future Evolution:
#' This function can also calculate expected ED given future evolution using a birth-death framework developed by (Mooers \emph{et al.}, 2012). The user must also enter an extinction rate (mu) and specation rate (lambda) in lineages per million species years and a timespan (tMa) in millions of years. The function calculates average expected new branch lengths (evolution in the future) for each tip and probabilites that lineages will go extinct within the timespan tMa. These values are incorportated into the calculation of expected ED. When considering future evolution, the initial extinction probabilities that are loaded into the function are the probabilities that the tips are extinct at 0 million years in the future (i.e. the present), not at some time in the distant future which is set by \code{tMa=}.
#'
#' @section Warning:
#' This function has been tested only on ultrametric, fully resolved phylogenetic trees. Technically, expected ED could be measured on non-ultrametric trees where branch lengths are scaled to something besides time (e.g. number of nucleotide substitutions) but results will be meaningless if you include future evolution. Use non-resovled and non-ultrametric trees at your own peril. Ultrametricity is checked by a call to \code{\link[ape]{is.ultrametric}} but the default tolerance has been set to 0.000001 because a phylogeny where tip-to-root distances vary by no more than 1 millionth of the age of the tree seems ultrametric enough.
#'
#' @section Acknowledgements:
#' This function uses code and internal functions from the \pkg{picante} (Kembel \emph{et al.}, 2010) and \pkg{ape} (Paradis \emph{et al.}, 2004) packages.
#'
#'
#' @return A list with components:
#' \itemize{
#' \item \strong{tip.values} A data.frame containing values for each tip
#' \describe{
#'   \item{Labels}{Names for each tip that match tree labels.}
#'   \item{Prob.Tip.Extinct.0}{The probability that a tip is extinct at the start of the analysis. These are the tip extinction probabilities that were entered by the user during the function call.}
#'   \item{Prob.Tip.Extinct.t}{The probability that a tip (or rather the lineage descended from that tip) is extinct by \code{tMa} into the future. If \code{tMa} is not designated, these probabilities will be the same as Prob.Tip.Extinct.0.}
#'   \item{Prob.Tip.Survive.t}{The probability that a tip (or rather the lineage descended from that tip) survives \code{tMa} into the future.}
#'   \item{Expected.Evol.Distinct.Ma.no.new}{The Expected Evolutionary Distinctiveness (expected ED) for a tip (in units of millions of years of unique evolution) without any new evolution into the future. This is the value most users will want to use for expected ED.}
#'   \item{Expected.Evol.Distinct.Ma}{Expected ED for each tip (or rather the lineage descended from that tip) that includes new branch lengths generated through evolution until \code{tMa} into the future. Measured in units millions of years of unique evolution. If \code{tMa=0}, this will be the same as Expected.Evol.Distinct.Ma.no.new as new evolution was not considered.}
#'   \item{Expected.Evol.Distinct.perc}{Expected ED for each tip (or rather the lineage descended from that tip) as a percent of total tree Phylogenetic Diversity (PD). In other words, what percentage of the total tree can be attributed to each tip. Note that this will include new evolution if \code{tMa} does not equal 0!}
#'   \item{Expected.Evol.Distinct.Ma.loss}{How much of the expected loss of PD can be attributed to each tip (in units of millions of years of unique evolution). Expected.Evol.Distinct.Ma is how much ED we expect a tip to retain, and Expected.Evol.Distinct.Ma.loss is how much ED we expect a tip to lose. Note that this will not include new evolution per se because new lineages that don't survive until \code{tMa} are not counted as losses. A tip can't lose more ED than it started out with at \code{tMa=0}. If \code{tMa=0}, the sum of each tip's Expected.Evol.Distinct.Ma and Expected.Evol.Distinct.Ma.loss will equal PD.}
#'   \item{Expected.Evol.Distinct.perc.loss}{The expected loss in ED for each tip as a percentage of total tree PD. Note that this will not include new evolution per se because new lineages that don't survive until \code{tMa} are not counted as losses. A tip can't lose more ED than it started out with at \code{tMa=0}.}
#'   \item{Pendant.Edge.Ma}{The length of a tip's pendant edge in units of millions of years of unique evolution. If tips are species, this is the age of the species. Note that this \strong{does not} include future evolution regardless of what \code{tMa} equals.}
#'   \item{Expected.Pendant.Edge.Ma.loss}{How much of a tip's pendant edge is expected to be lost (in units of millions of years of unique evolution). For example, if a species has a .74 probability of extinction by \code{tMa}, we would expect 74 \% of the pendant edge to be lost. Note that this \strong{does not} include branch lengths from future evolution regardless of what \code{tMa} equals but it will use the Prob.Tip.Extinct.t to calculate the expected loss \strong{not} Prob.Tip.Extinct.0 if \code{tMa} does not equal 0. This means that Expected.Pendant.Edge.Ma.loss is the amount of the tip's \emph{original} pendant edge that we would expect to lose given any future evolution.}
#'
#' }
#'
#'
#' \item \strong{edge.values} A data.frame containing values for each edge
#' \describe{
#'   \item{Node}{The name of the tipward node connected to the edge. This is where an edge ends.}
#'   \item{Daughters}{The number of daughter tips subtending a node.}
#'   \item{Edge.Length.Ma}{The length of the edge in units of millions of years.}
#'   \item{Node.Age.Ma}{The age of the node in millions of years in the past. This is the age where the edge starts, the rootward node of an edge.}
#'   \item{Prob.Edge.Extinct.t}{The probability that an edge will go extinct by \code{tMa} into the future. An edge only goes extinct if all tips subtending that edge go extinct.}
#'   \item{Prob.Edge.Survive.t}{The probability that an edge will survive until \code{tMa} into the future.}
#'   \item{Expected.Edge.Length.Ma}{The expected length of an edge in units of millions of years given the probability the edge will survive until \code{tMa}. For example, if an edge only has a 25 \% chance of surviving until \code{tMa}, we would only expect 25 \% of its length to remain.}
#'   \item{Expected.Edge.Length.Ma.loss}{The expected loss of length of an edge in units of millions of years given the probability the edge will survive until \code{tMa}. For example, if an edge only has a 25 \% chance of surviving until \code{tMa}, we would expect 75 \% of its length to be lost.}
#'
#' }
#'
#' \item \strong{tree} The original phylo object input
#'
#' \item \strong{lambda} The original speciation rate input
#'
#' \item \strong{mu} The original extinction rate input
#'
#' \item \strong{tMa} The original timespan input in millions of years
#'
#' \item \strong{source.of.data} Optional data tag from input
#'
#' }
#'
#'
#' @examples
#' data(bear_tree)
#' data(bear_probs)
#'
#' # Normal usage (without future evolution)
#' eED(tree=bear_tree, tip.extinction.probabilities=bear_probs)
#'
#' # Usage with future evolution
#' eED(tree=bear_tree, tip.extinction.probabilities=bear_probs, lambda=0.276, mu=0.272, tMa=2)

#'
#'
#' @author
#'  Matt Davis
#'
#' @references
#'  Faith, D. P. (2008). Threatened species and the potential loss of phylogenetic diversity: conservation scenarios based on estimated extinction probabilities and phylogenetic risk analysis. Conservation Biology, 22(6), 1461–1470.
#'
#'  Kembel, S. W., Cowan, P. D., Helmus, M. R., Cornwell, W. K., Morlon, H., Ackerly, D. D., et al. (2010). Picante: R tools for integrating phylogenies and ecology. Bioinformatics, 26(11), 1463–1464.
#'
#'  Mooers, A., Gascuel, O., Stadler, T., Li, H., & Steel, M. (2012). Branch lengths on birth–death trees and the expected loss of phylogenetic diversity. Systematic Biology, 61(2), 195–203.
#'
#'  Paradis, E., Claude, J. and Strimmer, K. (2004) APE: analyses of phylogenetics and evolution in R language. Bioinformatics, 20, 289–290.
#'
#'  Redding, D. W., Mazel, F., & Mooers, A. Ø. (2014). Measuring evolutionary isolation for conservation. PLoS ONE, 9(12), e113490.
#'
#' @export
#'
#'



#############################################################

#Expected Evolutionary Distinctiveness

#Function written by Matt Davis matt.davis@bios.au.dk


#Version 2.1 has been renamed eED to avoid namespace problems
#Version 2.0 has improved error checking and documentation. Auto plotting has been removed.
#Version 1.2 has optional plotting and auto save commands
#Version 1.1 has improvements to improve speed

###############################################################




eED <- function(tree=NA, tip.extinction.probabilities=NULL, lambda=NULL, mu=NULL, tMa=0, auto.save=F, source.of.data=NA){

  #What is the function version number
  version.number <- 2.1



  # Check inputs ##############################################################

  # Check the tree
  # Is there a tree?
  if(missing(tree)){stop("You need a phylogenetic tree to calculate expected ED")}

  # Is the tree in the right format?
  if(!ape::is.rooted.phylo(tree) | !ape::is.binary(tree)){stop("Tree must be a rooted, fully resolved phylogeney of class phylo")}

  # Is the tree ultrametric
  if(!ape::is.ultrametric(tree, tol=.000001)){warning("This function is not tested on non-ultrametric trees. Proceed at your own peril")}


  #####


  # Check the extinction probabilities
  # Are there extinction probabilities? If not, give every species an extinction probability of 0
  if(is.null(tip.extinction.probabilities)){
    warning("No extinction probabilities entered. All tips given 0 probability of going extinct")

    tip.extinction.probabilities <- rep(0, times=length(tree$tip.label))
    names(tip.extinction.probabilities) <- as.character(tree$tip.label)
  }

  # Are the extinction probabilities actually probabilities?
  if(any(tip.extinction.probabilities > 1) |
     any(tip.extinction.probabilities < 0) |
     !is.numeric(tip.extinction.probabilities)){
    stop("Check to make sure your tip extinction probabilities are from 0 to 1 and numeric")
  }

  # Do the names on the extinction probabilities match the names on the tree
  if(!all(names(tip.extinction.probabilities) %in% tree$tip.label) |
     !length(tip.extinction.probabilities)==length(tree$tip.label)){
    stop("Check that your tree tip names match the extinction probability names")
  }

  if(any(duplicated(names(tip.extinction.probabilities))) |
    any(duplicated(tree$tip.label))){
    stop("You have duplicate taxa names")
  }


  #####


  # Check to make sure the extinction and speciation rates are included
  if(tMa!=0 & (is.null(lambda) | is.null(mu) )){
    stop("You must provide lambda and mu values if tMa does not equal 0")
  }

  #If t = 0, make dummy lambda and mu values. These will be cancelled out by the t of 0 in equations so they will not affect results
  if(tMa==0 & (is.null(lambda) | is.null(mu) )){
    lambda <- 1
    mu <- .5
  }




  # Prep inputs ###############################################################

  thistree <- tree

  # Make sure the names of the extinction probabilities are in the same order as the tree tip labels
  tip.extinction.probabilities2 <- tip.extinction.probabilities[thistree$tip.label]

  tipprobs <- data.frame(Labels=names(tip.extinction.probabilities2),
                         Prob.Tip.Extinct.0=tip.extinction.probabilities2,
                         Prob.Tip.Extinct.t = NA,
                         Prob.Tip.Survive.t = NA,
                         Expected.Evol.Distinct.Ma.no.new = NA,
                         Expected.Evol.Distinct.Ma = NA,
                         Expected.Evol.Distinct.perc = NA,
                         Expected.Evol.Distinct.Ma.loss = NA,
                         Expected.Evol.Distinct.perc.loss = NA,
                         Pendant.Edge.Ma = NA,
                         Expected.Pendant.Edge.Ma.loss = NA)
  rownames(tipprobs) <- tipprobs$Labels
  #head(tipprobs)





  # Define internal functions ################################################

  # Make a function that gives the probability of a lineage surviving during the time period

  # From the proof of Theorem 4 in Mooers, A., Gascuel, O., Stadler, T., Li, H., & Steel, M. (2012). Branch lengths on birth–death trees and the expected loss of phylogenetic diversity. Systematic Biology, 61(2), 195–203.
  probsurvive <- function(lambda=NA, mu=NA, tMa=NA){

    # First some book keeping
    r <- lambda-mu
    e <- exp(1)

    # Now the theorem
    prob1 <- r/(lambda-(mu*(e^(-r*tMa))))
    return(prob1)


  }# End probsurvive function


  ###############################


  # From Theorem 4.i Mooers, A., Gascuel, O., Stadler, T., Li, H., & Steel, M. (2012). Branch lengths on birth–death trees and the expected loss of phylogenetic diversity. Systematic Biology, 61(2), 195–203.
  #This is half of the value put forth by Mooers et al. because we are starting from a single lineage instead of two lineages with length 0

  # This calculates the total expected length of a simple birth death tree
  theorem4 <-  function(lambda=NA, mu=NA, tMa=NA){

    # First some book keeping
    rho <- lambda/mu
    r <- lambda-mu
    e <- exp(1)

    # Now the frho function
    frho <- function(s){
      ((rho*(e^s))-1) / ((rho-1)*(e^s))
    }# End frho function

    # Now the theorem
    treelength <- (((e^(r*tMa))/mu) * log(frho(r*tMa)))

    return(treelength)




  }# End theorem4


  ##################################




  # Calculate tree values ##############################################

  # What is the expected average growth per lineage till time t
  #If t is 0, this length should be 0
  new.lineage.growth <- theorem4(lambda = lambda, mu=mu, tMa=tMa)

  # What is the probability of a tip that is already extant at time 0 surviving until time t
  # If t is 0, this probability should be 1
  prob.lineage.future <- probsurvive(lambda = lambda, mu=mu, tMa=tMa)

  # Calculate the probability of a tip surviving until time t
  tipprobs$Prob.Tip.Survive.t <- (1-tipprobs$Prob.Tip.Extinct.0)*prob.lineage.future

  # Calculate the probability of a tip going extinct before time t
  tipprobs$Prob.Tip.Extinct.t <- 1-tipprobs$Prob.Tip.Survive.t


  # Pull out those values as a vector to speed up downstream calculations. They will be in the same order as the tip labels
  tipprobs.extinct.t <- tipprobs[ , "Prob.Tip.Extinct.t"]

  #Get the number ot tips
  ntaxa <- length(thistree$tip.label)


  # Pull out the nodes
  # Pull out the second column as that is tipward end of an edge. The first column is the node where an edge starts (the rootward side)
  nodes <- thistree$edge[,2]
  numnodes <- length(nodes)
  #head(nodes)
  #View(thistree$edge)



  # Make a matrix for node memberships and extinction probabilities to go into
  Prob.Tip.Extinct.t <- matrix(data=NA, nrow=length(nodes), ncol=ntaxa, dimnames=list(nodes, thistree$tip.label))
  #head(Prob.Tip.Extinct.t)
  #View(Prob.Tip.Extinct.t)



  # Calculate node values ##############################################



  # Make a matrix to store information for each edge (ending at the listed node)
  variables <- c('Node', "Daughters", "Edge.Length.Ma", "Node.Age.Ma", "Prob.Edge.Extinct.t", "Prob.Edge.Survive.t",  "Expected.Edge.Length.Ma", "Expected.Edge.Length.Ma.loss", "sum.Prob.Tip.Extinct.t", "sum.Prob.Tip.Survive.t", "Expected.Edge.Length.rel", "Expected.Edge.Length.rel.loss")

  thesenodelengths <- matrix(data=NA, nrow=length(nodes), ncol=length(variables), dimnames=list(NULL, variables))
  #head(thesenodelengths)
  #View(thesenodelengths)

  # For progress report
  start.time <- proc.time()
  # Will take about an 30 seconds for a tree with 6,000 tips


  for (j in 1:length(nodes)) {

    # Now for every node
    # What is this node?
    thesenodelengths[j,'Node'] <- nodes[j]

    # Which nodes/tips subtend this node?
    sons <- picante::.node.desc(thistree, nodes[j])

    # Which tips subtend this node?
    # Fill in those slots in the matrix with the right probabilities
    Prob.Tip.Extinct.t[j, sons$tips] <- tipprobs.extinct.t[sons$tips]


    # How many tips subtend this node?
    thesenodelengths[j, "Daughters"] <- length(sons$tips)

    # Now what is the branchlength connecting to that node
    thesenodelengths[j, "Edge.Length.Ma"] <- thistree$edge.length[j]


    # Update progress report every 1000th node
    if(j %% 1000==0){
      loop.elapsed.time <- proc.time()-start.time
      time.per.node <- loop.elapsed.time/j

      time.left <- round((time.per.node*(numnodes-j))[3]/60, digits=1)

      message(paste0("##### Just finished node ", j, " of ", numnodes, ".  Possibly ", time.left, " minutes left"))
    }



  }# End j loop
  #View(thesenodelengths)

  # Get the probabilities of survival
  Prob.Tip.Survive.t <- 1-Prob.Tip.Extinct.t


  rownames(thesenodelengths) <- thesenodelengths[, "Node"]


  # Get the age of each node
  node.age <- ape::branching.times(thistree)

  # These nodes are labeled with the number of the ancester node = the node on the rootward side of an edge
  # The first column of the tree edge matrix corresponds to the rootward side of an edge = the node where the edge starts. So each node will have two edges subtending it
  # So we need to double the node ages and match them up with the order of the edge labels

  # This creates the order we need pull the nodes out and match them up with edges
  node.order <- match(x=thistree$edge[,1], table=names(node.age))
  #length(node.order)# Same as the number of edges we have
  thesenodelengths[, "Node.Age.Ma"] <- node.age[node.order]


  # Now multiply tip probabilites to get the probability of extinction of each edge
  # An edge will only go extinct if every tip subtending it goes extinct. So the probability of an edge going extinct is the product of all the extinction probabilities of the tips subtending it
  thesenodelengths[, "Prob.Edge.Extinct.t"] <-  matrixStats::rowProds(Prob.Tip.Extinct.t, na.rm=T)
  # Slightly slower version below
  #thesenodelengths[, "Prob.Edge.Extinct.t"] <-  apply(Prob.Tip.Extinct.t, 1, function(x) prod(x, na.rm=T))

  # The probability that the edge survives to time t
  thesenodelengths[, "Prob.Edge.Survive.t"] <-  1 - thesenodelengths[, "Prob.Edge.Extinct.t"]

  # Now we can calculate expected edge length by multiplying the original edge length by the probability that that edge is actually there
  thesenodelengths[, "Expected.Edge.Length.Ma"] <- thesenodelengths[, "Edge.Length.Ma"] * thesenodelengths[, "Prob.Edge.Survive.t"]

  # We can also calculate the expected loss of length for each edge by multiplying the length of the edge by the probability that the edge goes extinct
  thesenodelengths[, "Expected.Edge.Length.Ma.loss"] <- thesenodelengths[, "Edge.Length.Ma"] * thesenodelengths[, "Prob.Edge.Extinct.t"]

  # The sum of the tip extinction probabilities for dividing the edge length.  Why we need these is discussed below
  thesenodelengths[, "sum.Prob.Tip.Extinct.t"] <-  rowSums(Prob.Tip.Extinct.t, na.rm=T)

  # The sum of the tip survival probabilities for dividing the edge length. Why we need these is discussed below
  thesenodelengths[, "sum.Prob.Tip.Survive.t"] <-  rowSums(Prob.Tip.Survive.t, na.rm=T)

  # Now divide the expected length by the sum of tip survival probabilities to get one unit of relative length. This is so the length can be fairly split among tip probabilities later on
  # Note that some of these will come out NaN because we are dividing a zero probability of surival by the zero expected length. This taxa correctly won't have any contribution to expected length so it won't have evolutionary distinctiveness
  thesenodelengths[, "Expected.Edge.Length.rel"] <- thesenodelengths[, "Expected.Edge.Length.Ma"] / thesenodelengths[, "sum.Prob.Tip.Survive.t"]

  # Do the same for loss
  thesenodelengths[, "Expected.Edge.Length.rel.loss"] <- thesenodelengths[, "Expected.Edge.Length.Ma.loss"] / thesenodelengths[, "sum.Prob.Tip.Extinct.t"]

  # Quick sanity check, as the number of daughters increases, the probably of going extinct should drop down
  #plot(log(thesenodelengths[, "Daughters"]), log(thesenodelengths[ , "Prob.Edge.Extinct.t"]))



  # Calculate tip values ##############################################
  # Now calculate expected evolutionary distinctiveness loss for each taxa


  # So for each taxa

  # For progress report
  start.time2 <- proc.time()
  # Will take a few seconds

  for(q in 1:ntaxa){

    thistaxa <- thistree$tip.label[q]

    # Get all the edges a taxa subtends
    this.taxas.nodes <- nodes[!is.na( Prob.Tip.Extinct.t[, thistaxa])]

    # Now multiply the relative loss of length by extinction probability of this taxon to see how much of the lost length it is responsible for
    tipprobs[tipprobs$Labels==thistaxa, "Expected.Evol.Distinct.Ma.loss"] <- sum(tipprobs[thistaxa, "Prob.Tip.Extinct.t"] * thesenodelengths[nodes %in% this.taxas.nodes, "Expected.Edge.Length.rel.loss"], na.rm = T)

    # The sum of that is the expected evolutionary distinctiveness of the expected loss of PD for this taxa


    # Do the same for expected evolutionary distinctiveness (excluding new evolution of branch lengths)
    tipprobs[tipprobs$Labels==thistaxa, "Expected.Evol.Distinct.Ma.no.new" ] <- sum(tipprobs[thistaxa, "Prob.Tip.Survive.t"] * thesenodelengths[nodes %in% this.taxas.nodes, "Expected.Edge.Length.rel"], na.rm = T)

    # Get the expected growth of tree length per lineage to time t and multiply it by the probability that a tip even survived to time 0. Then add that to the expected evolutionary distinctiveness excluding new evolution
    tipprobs[tipprobs$Labels==thistaxa, "Expected.Evol.Distinct.Ma" ] <- tipprobs[tipprobs$Labels==thistaxa, "Expected.Evol.Distinct.Ma.no.new" ] + (new.lineage.growth * (1-tipprobs[tipprobs$Labels==thistaxa, "Prob.Tip.Extinct.0" ]))


    # Get the length of the pendant edge of this species = the species age = unique evolution
    # Find which edge has that species and then what that edgelength is
    tipprobs[tipprobs$Labels==thistaxa, "Pendant.Edge.Ma" ] <- thistree$edge.length[which(nodes==q)]



    # Update progress report every 500th taxa
    if(q%%500==0){
      loop.elapsed.time2 <- proc.time()-start.time2
      time.per.taxa <- loop.elapsed.time2/q

      time.left2 <- round((time.per.taxa*(ntaxa-q))[3]/60, digits=1)

      message(paste0("##### Just finished taxa ", q, " of ", ntaxa, ".  Possibly ", time.left2, " minutes left"))
    }

  }#End q loop



  # Check
  # PD should equal expected PD + the sum of expected evolutionary distinctiveness loss
  # Expected PD equals the sum of expected evolutionary distinctiveness without new evolution
  # So these two quantities should be the same
  #picante::pd(tree=thistree, samp=matrix(data=rep(1, ntaxa), nrow=1, dimnames=list(NULL, tipprobs$Labels) ))
  #sum(tipprobs$Expected.Evol.Distinct.Ma.no.new, tipprobs$Expected.Evol.Distinct.Ma.loss)


  # Calculate percentages
  tipprobs$Expected.Evol.Distinct.perc <- 100*tipprobs$Expected.Evol.Distinct.Ma/sum(tipprobs$Expected.Evol.Distinct.Ma)

  tipprobs$Expected.Evol.Distinct.perc.loss <- 100*tipprobs$Expected.Evol.Distinct.Ma.loss/sum(tipprobs$Expected.Evol.Distinct.Ma.loss)

  # Get the expected loss of length from the pendant edge
  tipprobs$Expected.Pendant.Edge.Ma.loss <- tipprobs$Pendant.Edge.Ma * tipprobs$Prob.Tip.Extinct.t




  # Save results #########################################################################

  # Remove some of the values used for internal calculations that will just be confusing to most users
  thesenodelengths2 <- as.data.frame(thesenodelengths)
  thesenodelengths2$sum.Prob.Tip.Extinct.t <- NULL
  thesenodelengths2$sum.Prob.Tip.Survive.t <- NULL
  thesenodelengths2$Expected.Edge.Length.rel <- NULL
  thesenodelengths2$Expected.Edge.Length.rel.loss <- NULL
  #head(thesenodelengths2)

  #If t = 0, don't include speciation rates and extinction rates as that will just be confusing.
  if(tMa==0){

    results <- list(tip.values=tipprobs, edge.values=thesenodelengths2, tree=tree, lambda=NA, mu=NA, tMa=tMa, source.of.data=source.of.data)

  }else{

    results <- list(tip.values=tipprobs, edge.values=thesenodelengths2, tree=tree, lambda=lambda, mu=mu, tMa=tMa, source.of.data=source.of.data)

    }


  #Save results if autosave is on
  if(auto.save==T){

    #Save results
    save(results, file=paste0("Results of eED v.", version.number, " with Lambda = ", lambda, ", Mu = ", mu, ", and t = ", tMa, " million years from ", source.of.data ))

  }#End save if statement




  return(results)

}#End eED function
