#' Calculate expected Phylogenetic Diversity (expected PD) for one or more communities
#'
#' \code{ePD} calculates the expected Phylogenetic Diversity (expected PD) for one or more communities of taxa and can include expected future evolution.
#'
#'
#' @param tree An object of class phylo.
#'
#' @param probabilities.tips.present.matrix A taxon by community matrix of presence probabilities for each tip of the tree. Each row is a taxon. Each column is a community (or temporal bin). Taxa names must match the tip labels of the tree. If no probabilities are supplied, all tips are assumed present and PD is calculated.
#'
#' @param lambda Numeric. A single instantaneous speciation rate in lineages per million species years.
#'
#' @param mu Numeric. A single instantaneous extinction rate in lineages per million species years.
#'
#' @param tMa Numeric. How many millions of years into the future will expected PD be calculated at? For normal function use where one is just looking at expected PD without any future evolution, tMa should be 0.
#'
#' @param auto.save Logical. Automatically write the output of the function to the working directory?
#'
#' @param save.tree Logical. Automatically save the tree used with the model output?
#'
#' @param save.tip.probabilites Logical. Automatically save the probabilites that tips are present in the model output?
#'
#' @param source.of.data Character. Optional data tag to include in the function output.
#'
#' @section Background: Phylogenetic Diversity (PD) (Faith, 1992) is a biodiversity metric that measures the total amount of evolutionary history in a community. It is defined as the sum of all branch lengths neccesary to connect a given set of taxa to the root of a phylogenetic tree and can be calculated in the function \code{\link[picante]{pd}}. PD assumes that all tips have a probability of 1 of being measured, though.
#'
#'   Expected Phylogenetic Diversity (expected PD) is the probabilistic implementation of PD (Faith, 2008). Expected PD is the expected amount of evolutionary history contained in a community. Each tip is given a probability of being present (from 0=absent to 1=present) that could reflect the taxon's actual survival probability (e.g. IUCN Red List Rank), output of a species distribution model, or the probability of sampling this taxon in a certain community. Note that if each tip is given a probability of 1 of being present, expected PD simplifies to PD and this function will return the same results as \code{\link[picante]{pd}}.
#'
#' @section Normal Use:
#' Typical usage is
#'
#' \code{ePD(tree, probabilities.tips.present.matrix)}
#'
#' This will calculate expected PD on a tree without any projected future evolution. Most users will only need to designate a tree and a matrix of probabilties that species are present named with labels that match the tip labels in the tree. The time \code{tMa=0} is set by default.
#'
#' Note that expected PD is simply the sum of each tip's expected Evolutionary Distinctiveness (expected ED) in a community. \code{ePD} is a convenience function to save you time if you don't care about individual taxon values. However, if you have already calculated expected Evolutionary Distinctiveness in \code{\link[mallorn]{eED}}, you can get the same results as \code{ePD} by using \code{sum(tip.values$Expected.Evol.Distinct.Ma)}. This may save you time if you are working with very large trees. If all your tip extinction probabilities are binary (0 or 1), you should just use \code{\link[picante]{pd}}, which will likely be much faster than \code{ePD}.
#'
#' \code{ePD} should not be confused with the \pkg{picante} function \code{\link[picante]{expected.pd}}. \code{\link[picante]{expected.pd}} calculates the amount of PD expected given taxonomic richness. It is not a implementation of Faith's (2008) probabilistic PD like \code{ePD}.
#'
#' @section Future Evolution:
#' This function can also calculate expected PD given future evolution using a birth-death framework developed by (Mooers \emph{et al.}, 2012). The user must also enter an extinction rate (mu) and specation rate (lambda) in lineages per million species years and a timespan (tMa) in millions of years. The function calculates average expected new branch lengths (evolution in the future) for each tip and probabilites that lineages will go extinct within the timespan tMa. These values are incorportated into the calculation of expected PD. When considering future evolution, the initial presence probabilities that are loaded into the function are the probabilities that the tips are present at 0 million years in the future (i.e. the present), not at some time in the distant future which is determined by the function iteself once \code{tMa} is set. Note that considering future evolution really only makes sense on large global phylogenies. This is not a feature that a typical user will need.
#'
#' @section Warning:
#' This function has been tested only on ultrametric, fully resolved phylogenetic trees. Technically, expected PD could be measured on non-ultrametric trees where branch lengths are scaled to something besides time (e.g. number of nucleotide substitutions) but results will be meaningless if you include future evolution. Use non-resovled and non-ultrametric trees at your own peril. Ultrametricity is checked by a call to \code{\link[ape]{is.ultrametric}} but the default tolerance has been set to 0.000001 because a phylogeny where tip-to-root distances vary by no more than 1 millionth of the age of the tree seems ultrametric enough.
#'
#' @section Acknowledgements:
#' This function uses code and internal functions from the \pkg{picante} (Kembel \emph{et al.}, 2010) and \pkg{ape} (Paradis \emph{et al.}, 2004) packages.
#'
#'
#' @return A list with components:
#' \itemize{
#' \item \strong{community.values} A matrix with expected Phylogenetic Diversity (expected PD) and expected species richness for each community
#' \describe{
#'   \item{Expected.Phylo.Diversity.Ma}{The expected PD in millions of years of evolution. This will include new evolution if \code{tMa} does not equal 0.}
#'   \item{Expected.Species.Diversity}{The expected taxonomic richness. This is really the expected number of tips so depending on the taxonomic resolution of your tree, this could be species, subspecies, or whole orders. This will include newly evolved tips if \code{tMa} does not equal 0.}
#'   }
#'
#'
#' \item \strong{lambda} The original speciation rate input
#'
#' \item \strong{mu} The original extinction rate input
#'
#' \item \strong{tMa} The original timespan input in millions of years
#'
#' \item \strong{source.of.data} Optional data tag from input
#'
#' \item \strong{probabilities.tips.present.matrix} The original taxon by community matrix of presence probabilities. Only present in ouput if \code{save.tip.probabilites=T}
#'
#' \item \strong{tree} The original phylo object input. Only present in output if \code{save.tree=T}
#'
#' }
#'
#'
#' @examples
#' data(bear_tree)
#' data(bear_matrix)
#'
#' # Normal usage (without future evolution)
#' ePD(tree=bear_tree, probabilities.tips.present.matrix=bear_matrix)
#'
#' # Usage with future evolution
#' # Note that it would not make sense to consider future evolution on a tree this small
#' ePD(tree=bear_tree, probabilities.tips.present.matrix=bear_matrix, lambda=0.276, mu=0.272, tMa=2)
#'
#'
#' @author
#'  Matt Davis
#'
#' @references
#'  Faith, D. P. (1992). Conservation evaluation and phylogenetic diversity. Biological Conservation, 61(1), 1–10.
#'
#'  Faith, D. P. (2008). Threatened species and the potential loss of phylogenetic diversity: conservation scenarios based on estimated extinction probabilities and phylogenetic risk analysis. Conservation Biology, 22(6), 1461–1470.
#'
#'  Kembel, S. W., Cowan, P. D., Helmus, M. R., Cornwell, W. K., Morlon, H., Ackerly, D. D., et al. (2010). Picante: R tools for integrating phylogenies and ecology. Bioinformatics, 26(11), 1463–1464.
#'
#'  Mooers, A., Gascuel, O., Stadler, T., Li, H., & Steel, M. (2012). Branch lengths on birth–death trees and the expected loss of phylogenetic diversity. Systematic Biology, 61(2), 195–203.
#'
#'  Paradis, E., Claude, J. and Strimmer, K. (2004) APE: analyses of phylogenetics and evolution in R language. Bioinformatics, 20, 289–290.
#'
#'
#' @export
#'
#'



#############################################################

#Expected Phylogenetic Diversity

#Function written by Matt Davis matt.davis@bios.au.dk

#Version 2.0 has inputs in probability of presence, not absence
#Version 1.0

###############################################################




ePD <- function(tree=NA, probabilities.tips.present.matrix=NULL, lambda=NULL, mu=NULL, tMa=0, auto.save=F, save.tree=F, save.tip.probabilites=F, source.of.data=NA){

  #What is the function version number
  version.number <- "2.0"


  # Check inputs ##############################################################

  # Check the tree
  # Is there a tree?
  if(missing(tree)){stop("You need a phylogenetic tree to calculate expected PD")}

  # Is the tree in the right format?
  if(!ape::is.rooted.phylo(tree) | !ape::is.binary(tree)){stop("Tree must be a rooted, fully resolved phylogeney of class phylo")}

  # Is the tree ultrametric
  if(!ape::is.ultrametric(tree, tol=.000001)){warning("This function is not tested on non-ultrametric trees. Proceed at your own peril")}


  #####


  # Check the presence probabilities
  # Are there presence probabilities? If not, stop
  if(is.null(probabilities.tips.present.matrix)){
    stop("You need to enter a taxon by commmunity matrix giving the probabilities that tips are present")
  }

  # Are the presence probabilities actually probabilities?
  if(any(probabilities.tips.present.matrix > 1) |
     any(probabilities.tips.present.matrix < 0) |
     !is.numeric(probabilities.tips.present.matrix)){
    stop("Check to make sure your tip probabilities are from 0 to 1 and numeric")
  }

  # Do the names on the probabilities match the names on the tree?
  if(!all(rownames(probabilities.tips.present.matrix) %in% tree$tip.label) |
     !dim(probabilities.tips.present.matrix)[1]==length(tree$tip.label)){
    stop("Check that your tree tip names match the row names for your tip probability matrix")
  }

  if(any(duplicated(rownames(probabilities.tips.present.matrix))) |
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


  # Change presence matrix to an extinction probability matrix
  tip.extinction.probabilities.matrix <- 1-probabilities.tips.present.matrix


  # How many communities are there?
  communities <- colnames(tip.extinction.probabilities.matrix)
  numcommunities <- length(communities)


  # Make sure the names of the extinction probabilities are in the same order as the tree tip labels
  # Add drop=F to maintain matrix structure
  tip.extinction.probabilities.matrix2 <- tip.extinction.probabilities.matrix[thistree$tip.label, , drop=F]


  tipprobs0 <- data.frame(Label=rownames(tip.extinction.probabilities.matrix2),
                         Prob.Tip.Extinct.0=NA,
                         Prob.Tip.Extinct.t = NA,
                         Prob.Tip.Survive.t = NA)
  rownames(tipprobs0) <- tipprobs0$Label
  #head(tipprobs0)



  # Make the matrix to store final results in
  community.values <- matrix(data=NA, nrow=numcommunities, ncol=2, dimnames=list(colnames(tip.extinction.probabilities.matrix2), c("Expected.Phylo.Diversity.Ma", "Expected.Species.Diversity")))



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

  #From Mooers, A., Gascuel, O., Stadler, T., Li, H., & Steel, M. (2012). Branch lengths on birth–death trees and the expected loss of phylogenetic diversity. Systematic Biology, 61(2), 195–203. (And many other authors who have used this equation)

  #This function calculates the average number of tips (species) in the reconstructed birth-death tree. Trees that don't make it to time t, i.e. number of tips = 0 are included in this average
  #This equation is half the value of that used by Mooers et al. because they are starting with two lineages

  expectedtips <-  function(lambda=NA, mu=NA, tMa=NA){

    numberoftips <- exp((lambda-mu)*tMa)

    return(numberoftips)

  }#End expectedtips



  # Calculate global tree values ##############################################



  # What is the expected average growth per lineage till time t
  #If t is 0, this length should be 0
  new.lineage.growth <- theorem4(lambda = lambda, mu=mu, tMa=tMa)

  # What is the probability of a tip that is already extant at time 0 surviving until time t
  # If t is 0, this probability should be 1
  prob.lineage.future <- probsurvive(lambda = lambda, mu=mu, tMa=tMa)

  # How many tips (species) will be produced for each existing lineage
  # If t is 0, this should be 1
  num.tips.per.lineage <- expectedtips(lambda = lambda, mu=mu, tMa=tMa)





  # Split into different communities ##############################################

  #Values for each community will be calculated in a big for loop

  for(i in 1:numcommunities){

    # Pull out this community
    thiscom <- communities[i]

    # Pull out the extinction probabilities
    tipprobs <- tipprobs0

    tipprobs$Prob.Tip.Extinct.0 <- tip.extinction.probabilities.matrix2[ ,i]






    # Calculate tree values ##############################################


    # Calculate the probability of a tip surviving until time t
    tipprobs$Prob.Tip.Survive.t <- (1-tipprobs$Prob.Tip.Extinct.0)*prob.lineage.future

    # Calculate the probability of a tip going extinct before time t
    tipprobs$Prob.Tip.Extinct.t <- 1-tipprobs$Prob.Tip.Survive.t


    # Pull out those values as a vector to speed up downstream calculations. They will be in the same order as the tip labels
    tipprobs.extinct.t <- tipprobs[ , "Prob.Tip.Extinct.t"]



    # To get the expected new evolution, multiply the probability that a tip even survives to the present (t=0) times the expected average new growth per lineage
    #If t is 0, these lengths should all be 0
    #The sum of these lengths is the total new evolution on the tree
    total.new.lineage.growth <- sum((1-tipprobs$Prob.Tip.Extinct.0) * new.lineage.growth)


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
    variables <- c('Node', "Edge.Length.Ma", "Prob.Edge.Extinct.t", "Prob.Edge.Survive.t",  "Expected.Edge.Length.Ma")

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


      # Now what is the branchlength connecting to that node
      thesenodelengths[j, "Edge.Length.Ma"] <- thistree$edge.length[j]


      # Update progress report every 1000th node
      if(j %% 1000==0){
        loop.elapsed.time <- proc.time()-start.time
        time.per.node <- loop.elapsed.time/j

        time.left <- round((time.per.node*(numnodes-j))[3]/60, digits=1)

        message(paste0("##### Just finished node ", j, " of ", numnodes, " for ", thiscom, ".  Possibly ", time.left, " minutes left for ", thiscom))
      }



    }# End j loop
    #View(thesenodelengths)


    rownames(thesenodelengths) <- thesenodelengths[, "Node"]


    # Now multiply tip probabilites to get the probability of extinction of each edge
    # An edge will only go extinct if every tip subtending it goes extinct. So the probability of an edge going extinct is the product of all the extinction probabilities of the tips subtending it
    thesenodelengths[, "Prob.Edge.Extinct.t"] <-  matrixStats::rowProds(Prob.Tip.Extinct.t, na.rm=T)


    # The probability that the edge survives to time t
    thesenodelengths[, "Prob.Edge.Survive.t"] <-  1 - thesenodelengths[, "Prob.Edge.Extinct.t"]

    # Now we can calculate expected edge length by multiplying the original edge length by the probability that that edge is actually there
    thesenodelengths[, "Expected.Edge.Length.Ma"] <- thesenodelengths[, "Edge.Length.Ma"] * thesenodelengths[, "Prob.Edge.Survive.t"]



    # The sum of the expected edge lengths in the expected phylogenetic diversity for the community without new evolution
    Expected.Phylo.Diversity.Ma.no.new <- sum(thesenodelengths[, "Expected.Edge.Length.Ma"])

    # So add in new evolution to get the total expected phylogenetic diversity
    community.values[i, "Expected.Phylo.Diversity.Ma"] <- Expected.Phylo.Diversity.Ma.no.new + total.new.lineage.growth

    # The sum of tip survival probabilities at t=0 is the number of expected tips. Multiply that by the number of tips generated per exisiting lineage to get the total tip (species) richness
    community.values[i, "Expected.Species.Diversity"] <- sum(1-tipprobs$Prob.Tip.Extinct.0) * num.tips.per.lineage


  }#End i community loop



  # Save results #########################################################################


  # If t = 0, don't include speciation rates and extinction rates as that will just be confusing.
  if(tMa==0){
    lambda <- NA
    mu <- NA
  }


  results <- list(community.values=community.values, lambda=lambda, mu=mu, tMa=tMa, source.of.data=source.of.data)

  if(save.tip.probabilites==T){results$probabilities.tips.present.matrix <- probabilities.tips.present.matrix}

  if(save.tree==T){results$tree <- tree}



  #Save results if autosave is on
  if(auto.save==T){

    #Save results
    save(results, file=paste0("Results of ePD v.", version.number, " with Lambda = ", lambda, ", Mu = ", mu, ", and t = ", tMa, " million years from ", source.of.data ))

  }#End save if statement




  return(results)

}#End ePD function


