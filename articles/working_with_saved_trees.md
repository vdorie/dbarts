# Working with dbarts Saved Trees

When a `dbarts` model is fit with the `keepTrees` option set to `TRUE`,
the sampler will return those trees in a “flat” format that corresponds
to a depth-first, left-hand traversal. Recursion can be used to map
functions the nodes in this structure, aggregating statistics of
interest.

## Simulation

To illustrate, we generate a fake dataset according to the “Friedman 1”
model (Friedman 1991).

``` r

f <- function(x)
    10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
        10 * x[,4] + 5 * x[,5]

set.seed(99)
sigma <- 1.0
n     <- 100

x <- matrix(runif(n * 10), n, 10)
y <- rnorm(n, f(x), sigma)

data <- data.frame(x, y)
```

## Model Fitting

In order to interrogate the trees, they must be saved when the model is
fit. This is accomplished by setting:

- For `bart`: `keeptrees = TRUE`
- For `bart2`: `keepTrees = TRUE`
- For a custom `dbartsSampler`,
  `control = dbartsControl(keepTrees = TRUE)`

In the context of our fake data, with small sample numbers for
illustrative purposes, a model can be fit thusly:

``` r

library(dbarts, quietly = TRUE)

bartFit <- bart(
    y ~ ., data,
    ndpost = 4,   # number of posterior samples
    nskip = 1000, # number of "warmup" samples to discard
    nchain = 2,   # number of independent, parallel chains
    nthread = 1,  # units of parallel execution
    ntree = 3,    # number of trees per chain
    seed = 56,    # chosen to generate a deep tree
    keeptrees = TRUE,
    verbose = FALSE)
```

## Extracting Trees

The extract function accepts as a `type` the value `"trees"`. If
present, the arguments `chainNums`, `sampleNums`, and/or `treeNum`s can
be used to extract only a subset of trees.

``` r

trees <- extract(bartFit, "trees")
```

## Flattened Trees

The `trees` data frame corresponds to a depth-first, left-hand-side tree
traversal.

``` r

print(head(trees, n = 10))
```

    ##    chain sample tree   n var        value
    ## 1      1      1    1 100   5  0.172099022
    ## 2      1      1    1  19  -1 -0.152703748
    ## 3      1      1    1  81   2  0.553265862
    ## 4      1      1    1  34   2  0.218104606
    ## 5      1      1    1   9   3  0.706622887
    ## 6      1      1    1   7  -1 -0.222561794
    ## 7      1      1    1   2  -1 -0.004843434
    ## 8      1      1    1  25   5  0.897998292
    ## 9      1      1    1  23  -1  0.011393382
    ## 10     1      1    1   2  -1  0.008832009

The columns refer to:

- `chain`, `sample`, `tree` - index variables
- `n` - number of observations in node (from the training data, or from
  `newdata` when it is supplied; see “Counts for New Data” below)
- `var` - either the index of the variable used for splitting or -1 if
  the node is a leaf
- `value` - either the value such that observations less than or equal
  to it are sent down the left path of the tree or the predicted value
  for a leaf node

The mapping between the values of `var` and the variable names can be
looked up in the internal copy of the data that the sampler stores. This
can be found in a fitted model as the element `fit$data@x`, as seen
below.

## Tree Traversal

A useful technique for processing trees is
[recursion](https://en.wikipedia.org/wiki/Recursion_%28computer_science%29).
By having the function return the number of nodes it has processed, it
is possible to advance from the left-hand-side to the right by skipping
ahead the appropriate number of rows. For example:

``` r

# Turns a flatted tree data frame into a list of lists, or a "natural" tree
# structure.
rebuildTree <- function(tree, object) {
    # Define a worker function that will be recursively called on every node.
    rebuildTreeRecurse <- function(tree) {
        node <- list(
            value = tree$value[1],
            n     = tree$n[1]
        )
        # Check node if is a leaf, and if so return early.
        if (tree$var[1] == -1) {
            node$n_nodes <- 1
            return(node)
        }
        
        node$var <- variableNames[tree$var[1]]
        
        # By removing the current row, we can recurse down the left branch.
        headOfLeftBranch <- tree[-1,]
        left <- rebuildTreeRecurse(headOfLeftBranch)
        n_nodes.left <- left$n_nodes
        left$n_nodes <- NULL
        node$left <- left
        
        # The right branch is obtained by advancing past the left nodes.
        headOfRightBranch <- tree[seq.int(2 + n_nodes.left, nrow(tree)),]
        right <- rebuildTreeRecurse(headOfRightBranch)
        n_nodes.right <- right$n_nodes
        right$n_nodes <- NULL
        node$right <- right
        
        node$n_nodes <- 1L + n_nodes.left + n_nodes.right
        
        return(node)
    }
    variableNames <- colnames(object$fit$data@x)
    
    result <- rebuildTreeRecurse(tree)
    result$n_nodes <- NULL
    return(result)
}

treeOfInterest <- subset(trees, chain == 1 & sample == 3 & tree == 1)
print(rebuildTree(treeOfInterest, bartFit))
```

    ## $value
    ## [1] 0.172099
    ## 
    ## $n
    ## [1] 100
    ## 
    ## $var
    ## [1] "X5"
    ## 
    ## $left
    ## $left$value
    ## [1] -0.125139
    ## 
    ## $left$n
    ## [1] 19
    ## 
    ## 
    ## $right
    ## $right$value
    ## [1] 0.5532659
    ## 
    ## $right$n
    ## [1] 81
    ## 
    ## $right$var
    ## [1] "X2"
    ## 
    ## $right$left
    ## $right$left$value
    ## [1] 0.2181046
    ## 
    ## $right$left$n
    ## [1] 34
    ## 
    ## $right$left$var
    ## [1] "X2"
    ## 
    ## $right$left$left
    ## $right$left$left$value
    ## [1] 0.7066229
    ## 
    ## $right$left$left$n
    ## [1] 9
    ## 
    ## $right$left$left$var
    ## [1] "X3"
    ## 
    ## $right$left$left$left
    ## $right$left$left$left$value
    ## [1] 0.2277311
    ## 
    ## $right$left$left$left$n
    ## [1] 7
    ## 
    ## $right$left$left$left$var
    ## [1] "X9"
    ## 
    ## $right$left$left$left$left
    ## $right$left$left$left$left$value
    ## [1] -0.1986161
    ## 
    ## $right$left$left$left$left$n
    ## [1] 2
    ## 
    ## 
    ## $right$left$left$left$right
    ## $right$left$left$left$right$value
    ## [1] -0.2271963
    ## 
    ## $right$left$left$left$right$n
    ## [1] 5
    ## 
    ## 
    ## 
    ## $right$left$left$right
    ## $right$left$left$right$value
    ## [1] -0.003246283
    ## 
    ## $right$left$left$right$n
    ## [1] 2
    ## 
    ## 
    ## 
    ## $right$left$right
    ## $right$left$right$value
    ## [1] 0.8979983
    ## 
    ## $right$left$right$n
    ## [1] 25
    ## 
    ## $right$left$right$var
    ## [1] "X5"
    ## 
    ## $right$left$right$left
    ## $right$left$right$left$value
    ## [1] 0.02943045
    ## 
    ## $right$left$right$left$n
    ## [1] 23
    ## 
    ## 
    ## $right$left$right$right
    ## $right$left$right$right$value
    ## [1] 0.04386729
    ## 
    ## $right$left$right$right$n
    ## [1] 2
    ## 
    ## 
    ## 
    ## 
    ## $right$right
    ## $right$right$value
    ## [1] 0.1121237
    ## 
    ## $right$right$n
    ## [1] 47

Using a `by` statement, it is possible to “rebuild” all trees at once:

``` r

allTrees <- by(
    data    = trees,
    INDICES = trees[,c("chain", "sample", "tree")],
    FUN     = rebuildTree, 
    object  = bartFit)

# One way to index the result of this:
#    allTrees[chain = "1", sample = "2", tree = "3"]
```

## Plotting Trees

`dbartsSampler` objects have a `plotTree` method that can be used to
visualize single trees:

``` r

bartFit$fit$plotTree(chainNum = 1, sampleNum = 3, treeNum = 1)
```

![](working_with_saved_trees_files/figure-html/plotTree-1.png)

## Getting Tree Predictions

The following function traverses a flattened tree, splits observations
while going the branches, and populates a vector giving the predicted
value of that tree on input data. It requires a data in the same format
as the fitted bart model so that it can evaluate the splits.

``` r

getPredictionsForTree <- function(tree, x) {
    predictions <- rep(NA_real_, nrow(x))
    
    getPredictionsForTreeRecursive <- function(tree, indices) {
        if (tree$var[1] == -1) {
            # Assigns in the calling environment by using <<-
            predictions[indices] <<- tree$value[1]
            return(1)
        }

        goesLeft <- x[indices, tree$var[1]] <= tree$value[1]
        headOfLeftBranch <- tree[-1,]
        n_nodes.left <- getPredictionsForTreeRecursive(
            headOfLeftBranch, indices[goesLeft])
        
        headOfRightBranch <- tree[seq.int(2 + n_nodes.left, nrow(tree)),]
        n_nodes.right <- getPredictionsForTreeRecursive(
            headOfRightBranch, indices[!goesLeft])
        
        return(1 + n_nodes.left + n_nodes.right)
    }

    getPredictionsForTreeRecursive(tree, seq_len(nrow(x)))

    return(predictions)
}

getPredictionsForTree(treeOfInterest, bartFit$fit$data@x[1:5,])
```

    ## [1] -0.2271963  0.1121237  0.1121237 -0.1251390 -0.1251390

A `by` statement can be used to obtain all predictions for all trees.

## Counts for New Data

The recursion above routes new data through the trees by hand. To simply
count how a new dataset distributes across the (frozen) partition, pass
`newdata` to `extract`: the tree structure, split rules, and leaf
predictions are unchanged, and only the `n` column is recomputed by
routing that data through each tree. The predictors are supplied in the
same form as for `predict`.

``` r

newData <- x[1:5,]
newTrees <- extract(bartFit, "trees", newdata = newData, sampleNums = 3, treeNums = 1)
print(subset(newTrees, chain == 1))
```

    ##    chain sample tree n var        value
    ## 1      1      3    1 5   5  0.172099022
    ## 2      1      3    1 2  -1 -0.125139008
    ## 3      1      3    1 3   2  0.553265862
    ## 4      1      3    1 1   2  0.218104606
    ## 5      1      3    1 1   3  0.706622887
    ## 6      1      3    1 1   9  0.227731139
    ## 7      1      3    1 0  -1 -0.198616087
    ## 8      1      3    1 1  -1 -0.227196258
    ## 9      1      3    1 0  -1 -0.003246283
    ## 10     1      3    1 0   5  0.897998292
    ## 11     1      3    1 0  -1  0.029430446
    ## 12     1      3    1 0  -1  0.043867291
    ## 13     1      3    1 2  -1  0.112123719

## Advanced Traversal

The following function can be used to map an arbitrary function over a
tree.

``` r

mapOverNodes <- function(tree, f, ...) {
    mapOverNodesRecurse <- function(tree, depth, f, ...) {
        node <- list(
            value = tree$value[1],
            n = tree$n[1],
            depth = depth
        )
        if (tree$var[1] == -1) {
            node$n_nodes <- 1
            node$f.x <- f(node, ...)
            return(node)
        }
        node$var <- tree$var[1]
        node$f.x <- f(node, ...)
        
        headOfLeftBranch <- tree[-1,]
        left <- mapOverNodesRecurse(headOfLeftBranch, depth + 1, f, ...)
        n_nodes.left <- left$n_nodes
        left$n_nodes <- NULL
        node$left <- left

        
        headOfRightBranch <- tree[seq.int(2 + n_nodes.left, nrow(tree)),]
        right <- mapOverNodesRecurse(headOfRightBranch, depth + 1, f, ...)
        n_nodes.right <- right$n_nodes
        right$n_nodes <- NULL
        node$right <- right
        
        node$n_nodes <- 1 + n_nodes.left + n_nodes.right
        return(node)
    }
    result <- mapOverNodesRecurse(tree, 1, f, ...)
    result$n_nodes <- NULL
    return(result)
}
```

As an example of its usage, the following function aggregates all
ancestor/descendant relationships in a tree. In a data object - here an
R environment - it keeps track of the current state of traversal. This
includes the variables that have been used for splits above the current
node and also includes the current node’s depth, which is used to detect
backtracking.

``` r

observeInteractions <- function(node, ...) {
    if (is.null(node$var)) return(NULL)

    interactionData <- list(...)$interactionData
    # Make the current node visibile inside the environment.
    interactionData$node <- node
    with(interactionData, {
        if (node$depth <= currentDepth) {
            # If true, we have backtracked to go down the right branch, so we
            # remove the variables from the left branch.
            currentVariables <- currentVariables[seq_len(node$depth - 1)]
        }
        if (length(interactionData$currentVariables) > 0) {
            # This is a brute-force way of updating the following indices,
            # relying on the column-major storage order that R uses:
            #     hasInteraction[currentVariables,,drop = FALSE][,node$var]
            updateIndices <- currentVariables +
                (node$var - 1) * nrow(hasInteraction)
            hasInteraction[updateIndices] <- TRUE
        }
        currentVariables <- c(currentVariables, node$var)
        currentDepth <- node$depth
    })
    rm("node", envir = interactionData)
    
    # Since the function is used for its side effects, there isn't a return
    # value.
    return(NULL)
}

numVariables  <- ncol(bartFit$fit$data@x)
variableNames <- colnames(bartFit$fit$data@x)

# Define this as an environment as they are mutable
interactionData <- list2env(list(
    currentDepth = 0,
    currentVariables = integer(),
    hasInteraction = matrix(
        data = FALSE,
        ncol = numVariables, nrow = numVariables,
        dimnames = list(ancestor = variableNames, descendant = variableNames)
    )
))

invisible(mapOverNodes(
    treeOfInterest,
    observeInteractions,
    interactionData = interactionData
))
```

After execution, the boolean matrix in the interaction data environment
will contain all ancestor/descendant relationships in this tree.

``` r

print(interactionData$hasInteraction)
```

    ##         descendant
    ## ancestor    X1    X2    X3    X4    X5    X6    X7    X8    X9   X10
    ##      X1  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ##      X2  FALSE  TRUE  TRUE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE
    ##      X3  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE
    ##      X4  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ##      X5  FALSE  TRUE  TRUE FALSE  TRUE FALSE FALSE FALSE  TRUE FALSE
    ##      X6  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ##      X7  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ##      X8  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ##      X9  FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ##      X10 FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE

Finally, the entire process can be wrapped in a `by` statement and the
results aggregated across all trees in order to count the number of
times variables have specific relationships.

## References

Friedman, Jerome H. 1991. “Multivariate Adaptive Regression Splines.”
*The Annals of Statistics* 19 (1): 1–67.
<https://doi.org/10.1214/aos/1176347963>.
