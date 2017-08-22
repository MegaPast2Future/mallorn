# mallorn
`mallorn` is an R package for calculating expected Phylogenetic Diversity and Evolutionary Distinctiveness

This package provides functions for calculating probabilistic phylogenetic diversity metrics. These functions are designed to work with extinction probabilities but they would also work in cases where species have variable probabilities of being found in a plot or grid cell.

Warning! This package is still in developement. Use it with caution.

To install it from its
[GitHub repository](https://github.com/ursus-americanus/mallorn). You first need to
install the [devtools](https://github.com/hadley/devtools) package.

```r
install.packages("devtools")
```

Then install `mallorn` using the `install_github` function in the
[devtools](https://github.com/hadley/devtools) package.

```r
library(devtools)
install_github("ursus-americanus/mallorn")
```

#### Example use

To calculate expected Phylogenetic Diversity, use the `ePD` function. Expected PD is the expected amount of evolutionary history represeneted in a community.

```r
library(ape)
library(mallorn)

data(bear_tree)
data(bear_matrix)

ePD(tree=bear_tree, tip.extinction.probabilities.matrix=bear_matrix)
```


To calculate expected Evolutionary Distinctiveness, use the `eED` function. Expected ED is the expected amount of evolutionary history represented by each individual taxon.

```r
library(ape)
library(mallorn)

data(bear_tree)
data(bear_probs)

eED(tree=bear_tree, tip.extinction.probabilities=bear_probs)
```


`mallorn` can also predict expected PD and ED given future evolution.

```r
library(ape)
library(mallorn)

data(bear_tree)
data(bear_probs)
data(bear_matrix)

ePD(tree=bear_tree, tip.extinction.probabilities.matrix=bear_matrix, lambda=0.276, mu=0.272, tMa=2)

eED(tree=bear_tree, tip.extinction.probabilities=bear_probs, lambda=0.276, mu=0.272, tMa=2)
```


The output of `eED` can be used to plot phylogenies where each edge is colored by its probability of extinction.

```r
library(ape)
library(mallorn)

# To install the ggtree package from Bioconductor
source("https://bioconductor.org/biocLite.R")
biocLite("ggtree")
library(ggtree)

data(bear_tree)
data(bear_probs)

# Calculate expected ED
res <- eED(tree=bear_tree, tip.extinction.probabilities=bear_probs)

# Extract internal node and edge values
edge.values <- res$edge.values

# Have to change the name of the label column so that data matches up with ggtree
colnames(edge.values)[1] <- "node"


# Create ggtree object
p <- ggtree(bear_tree, layout="rectangular", aes(color=Prob.Edge.Extinct.t), size=1)

#  Put in the edge values
p2 <- p %<+% as.data.frame(edge.values)

# Construct plot
p3 <- p2+
  geom_tiplab()+
  # Add a color scale to show extinction probability
  scale_color_gradient2(name="Probability\nof extinction", low="dodgerblue2", mid="bisque", high="red1", midpoint=.50, guide="colorbar")+
  # Must explicitly place the legend
  theme(legend.position="right")+
  # Leave some room for tip labels
  xlim_tree(max(branching.times(bear_tree))*1.7)

print(p3)
```
