hc <- function(data, node.sizes, scoring.func = 0, cpc, cont.nodes = c(), ess = 1, tabu.tenure = 100,
               max.parents = length(node.sizes) - 1,
               init.net = NULL, wm.max = 15, layering = NULL, layer.struct = NULL,
               mandatory.edges = NULL) {
  n.nodes <- ncol(data)
  n.cases <- nrow(data)

  # just to be sure
  storage.mode(node.sizes) <- "integer"

  # quantize data of continuous nodes
  levels <- rep(0, n.nodes)
  levels[cont.nodes] <- node.sizes[cont.nodes]

  # data <- quantize.with.na.matrix( data, levels )
  # data <- quantize.matrix( data, levels )
  out.data <- quantize.matrix(data, levels)
  data <- out.data$quant
  # quantiles(bn) <- out.data$quantiles
  # quantiles(dataset) <- out.data$quantiles

  # apply mandatory edges (if present)
  m.edges <- matrix(0L, n.nodes, n.nodes)
  if (!is.null(mandatory.edges)) {
    m.edges <- m.edges | mandatory.edges
  }
  storage.mode(m.edges) <- "integer"

  # start with init.net if not NULL, otherwise with empty matrix
  if (!is.null(init.net)) {
    curr.g <- init.net
    # just to be sure
    storage.mode(curr.g) <- "integer"
    # add init edges out of cpc
    cpc <- cpc | init.net | t(init.net)
  } else {
    curr.g <- matrix(0L, n.nodes, n.nodes)
  } # integers!
  curr.g <- curr.g | m.edges

  # apply layering
  if (!is.null(layer.struct) && is.null(layering)) {
    stop("layer.struct provided without layering.\n")
  }
  if (!is.null(layering) && length(unique(layering)) > 1) {
    n.layers <- length(unique(layering))

    if (is.null(layer.struct)) {
      layer.struct <- matrix(1L, n.layers, n.layers)
      layer.struct[lower.tri(layer.struct)] <- 0
      layer.struct[1, 1] <- 0
    }

    layers <- matrix(1L, n.nodes, n.nodes)
    for (i in 1:n.layers) {
      for (j in 1:n.layers) {
        layers[which(layering == i), which(layering == j)] <- layer.struct[i, j]
      }
    }
    diag(layers) <- 0

    # keep only edges allowed by both the CPC / initial network and the layering
    cpc <- cpc & layers

    if (sum(m.edges | t(m.edges)) > 0L && sum(curr.g & layers) == 0L) {
      stop("the mandatory edges provided are inconsistent with the given layering.\n")
    }
  }
  # end apply layering

  curr.score.nodes <- array(0, n.nodes)
  for (i in 1L:n.nodes) {
    curr.score.nodes[i] <- .Call("bnstruct_score_node", data, node.sizes, i - 1L, which(curr.g[, i] != 0) - 1L,
      scoring.func, ess,
      PACKAGE = "bnstruct"
    )
  }

  # global best solution
  global.best.g <- curr.g
  global.best.score <- sum(curr.score.nodes)

  # tabu list
  tabu <- array(0L, c(n.nodes, n.nodes, tabu.tenure))
  tabu.pt <- 1

  # worsening moves
  wm.count <- 0
  # wm.max <- 15
  while (wm.count < wm.max) {
    next.score.diff <- rep(-Inf, n.nodes)
    next.pert <- rep(-1L, n.nodes)
    # print(tabu[,,1:10])
    # try all possible perturbations
    for (node in 1L:n.nodes)
    {
      for (par in 1L:n.nodes)
      {
        if (cpc[par, node] == 1L && (m.edges[par, node] == 0 && m.edges[node, par] == 0)) {
          next.g <- curr.g
          s.diff <- -Inf

          if (curr.g[par, node] == 1L) # edge removal
            {
              next.g[par, node] <- 0L
              if (.Call("bnstruct_is_acyclic", next.g, PACKAGE = "bnstruct") & !.Call("bnstruct_in_tabu", next.g, tabu, PACKAGE = "bnstruct")) {
                #               cat(node.sizes,'\t',node-1L,'\t',which(next.g[,node]!=0)-1L,'\n')
                s.diff <- .Call("bnstruct_score_node", data, node.sizes, node - 1L, which(next.g[, node] != 0) - 1L,
                  scoring.func, ess,
                  PACKAGE = "bnstruct"
                ) - curr.score.nodes[node]
              }
            }
          # edge addition
          # check also if there is room for one more parent
          else if (curr.g[node, par] == 0L && sum(curr.g[, node]) < max.parents) {
            next.g[par, node] <- 1L
            # print(c(node,par,.Call("is_acyclic", next.g, PACKAGE = "bnstruct"),!.Call("in_tabu", next.g, tabu, PACKAGE = "bnstruct")))
            if (.Call("bnstruct_is_acyclic", next.g, PACKAGE = "bnstruct") & !.Call("bnstruct_in_tabu", next.g, tabu, PACKAGE = "bnstruct")) {
              # print("here\n");

              #               cat(node.sizes,'\t',node-1L,'\t',which(next.g[,node]!=0)-1L,'\n')
              s.diff <- .Call("bnstruct_score_node", data, node.sizes, node - 1L, which(next.g[, node] != 0) - 1L,
                scoring.func, ess,
                PACKAGE = "bnstruct"
              ) - curr.score.nodes[node]
            }
          }
          # edge reversal
          # check also if there is room for one more parent
          else if (sum(curr.g[, node]) < max.parents) {
            next.g[par, node] <- 1L
            next.g[node, par] <- 0L
            if (.Call("bnstruct_is_acyclic", next.g, PACKAGE = "bnstruct") & !.Call("bnstruct_in_tabu", next.g, tabu, PACKAGE = "bnstruct")) {
              #               cat(node.sizes,'\t',node-1L,'\t',which(next.g[,node]!=0)-1L,'\t',par-1L,'\t',which(next.g[,par]!=0)-1L,'\n')
              s.diff <- .Call("bnstruct_score_node", data, node.sizes, node - 1L, which(next.g[, node] != 0) - 1L,
                scoring.func, ess,
                PACKAGE = "bnstruct"
              ) +
                .Call("bnstruct_score_node", data, node.sizes, par - 1L, which(next.g[, par] != 0) - 1L,
                  scoring.func, ess,
                  PACKAGE = "bnstruct"
                ) -
                (curr.score.nodes[node] + curr.score.nodes[par])
            }
          }
          # test for local improvement
          if (s.diff > next.score.diff[node]) {
            next.score.diff[node] <- s.diff
            next.pert[node] <- par
          }
        }
      }
    }

    best.node <- which.max(next.score.diff)

    # no possible improvements given the tabu list
    if (next.score.diff[best.node] == -Inf) {
      break
    }

    # update current graph and scores
    curr.g[next.pert[best.node], best.node] <- 1L - curr.g[next.pert[best.node], best.node] # flip
    if (curr.g[best.node, next.pert[best.node]] == 1L) # need to reverse
      {
        curr.g[best.node, next.pert[best.node]] <- 0L
        # cat(node.sizes,'\t',best.node-1L,'\t',which(curr.g[,best.node]!=0)-1L,'\t',next.pert[best.node]-1L,'\t',which(curr.g[,next.pert[best.node]]!=0)-1L,'\n')
        curr.score.nodes[best.node] <- .Call("bnstruct_score_node", data, node.sizes, best.node - 1L,
          which(curr.g[, best.node] != 0) - 1L, scoring.func, ess,
          PACKAGE = "bnstruct"
        )
        curr.score.nodes[next.pert[best.node]] <- .Call("bnstruct_score_node", data, node.sizes, next.pert[best.node] - 1L,
          which(curr.g[, next.pert[best.node]] != 0) - 1L, scoring.func, ess,
          PACKAGE = "bnstruct"
        )
      } else {
      curr.score.nodes[best.node] <- curr.score.nodes[best.node] + next.score.diff[best.node]
    }

    # print(c(tabu.pt,best.node,next.pert[best.node],sum(curr.score.nodes)))

    if (global.best.score < sum(curr.score.nodes)) # check for global best
      {
        wm.count <- 0
        global.best.g <- curr.g
        global.best.score <- sum(curr.score.nodes)
      } else {
      wm.count <- wm.count + 1
    }

    # update tabu list
    tabu[, , tabu.pt] <- curr.g
    tabu.pt <- (tabu.pt) %% tabu.tenure + 1
    # print(curr.g)
  }

  return(global.best.g)
}


mmpc <- function(data, node.sizes, cont.nodes = NULL, chi.th = 0.05,
                 layering = NULL, layer.struct = NULL,
                 max.fanin = length(node.sizes) - 1, mandatory.edges = NULL) {
  n.nodes <- ncol(data)
  n.cases <- nrow(data)
  min.counts <- 5

  if (length(max.fanin) == 1) {
    max.fanin <- rep(max.fanin, n.nodes)
  }

  # just to be sure
  storage.mode(node.sizes) <- "integer"

  # quantize data of continuous nodes
  levels <- rep(0, n.nodes)
  levels[cont.nodes] <- node.sizes[cont.nodes]

  # data <- quantize.with.na.matrix( data, levels )
  # data <- quantize.matrix( data, levels )
  out.data <- quantize.matrix(data, levels)
  data <- out.data$quant
  # quantiles(bn) <- out.data$quantiles
  # quantiles(dataset) <- out.data$quantiles

  # default values for layering
  if (is.null(layering)) {
    layering <- rep.int(1, n.nodes)
    if (!is.null(layer.struct)) {
      stop("Argument layer.struct without layering\n")
    }
  }
  n.layers <- length(unique(layering))
  if (is.null(layer.struct)) {
    layer.struct <- matrix(1, n.layers, n.layers)
    if (n.layers > 1) {
      layer.struct[lower.tri(layer.struct)] <- 0
      layer.struct[1, 1] <- 0 # default: no edges between nodes at level 1
    }
  }
  # print(n.layers)
  # print(layering)
  # print(layer.struct)

  # apply layering
  layer.mat <- matrix(1, n.nodes, n.nodes)
  for (i in 1:n.layers) {
    for (j in 1:n.layers) {
      layer.mat[layering == i, layering == j] <- layer.struct[i, j]
    }
  }
  diag(layer.mat) <- 0

  cpc.mat <- layer.mat | t(layer.mat) # constrain the cpc search
  allowed <- cpc.mat

  # mandatory edges
  if (is.null(mandatory.edges)) {
    medges <- matrix(0, n.nodes, n.nodes)
  } else {
    medges <- mandatory.edges | t(mandatory.edges)
  }

  # print(cpc.mat)

  storage.mode(data) <- "integer"

  # forward addition of nodes
  for (i in 1:n.nodes)
  {
    cpc.mat[i, ] <- mmpc.fwd(data, node.sizes, allowed, i, chi.th, min.counts, max.fanin)
    # cat("cpcMat ",i,": ",cpc.mat[i,],"\n")
    allowed[, i] <- allowed[, i] & t(cpc.mat[i, ])
  }

  # print(cpc.mat)

  # backwards removal of nodes
  for (i in 1:n.nodes)
  {
    cpc.mat[i, ] <- mmpc.bwd(data, node.sizes, cpc.mat[i, ], i, chi.th, min.counts, max.fanin)
    # cat("cpcMat ",i,": ",cpc.mat[i,],"\n")
    cpc.mat[, i] <- cpc.mat[, i] & t(cpc.mat[i, ])
  }

  # print(cpc.mat)

  # symmetry enformcement
  cpc.mat <- cpc.mat * t(cpc.mat)

  # further filter with layering
  cpc.mat <- cpc.mat * layer.mat

  # force mandatory edges
  cpc.mat <- cpc.mat | medges
  storage.mode(cpc.mat) <- "integer"

  return(cpc.mat)
}

mmpc.fwd <- function(data, node.sizes, allowed, x, chi.th, min.counts, max.fanin) {
  # cat("\n",x,":\n");
  n.nodes <- ncol(data)
  n.cases <- nrow(data)

  # test without conditioning
  # print(chi.th)
  minAssoc <- rep(0, n.nodes)
  for (y in 1:n.nodes) {
    if (allowed[x, y]) {
      minAssoc[y] <- g2(data, node.sizes, x, y, chi.th, min.counts = min.counts)
    }
  }
  allowed[x, minAssoc == 0] <- 0 # remove already independent nodes

  m <- max(minAssoc)
  if (m == 0) {
    return(rep(0, n.nodes))
  }

  m.ind <- which.max(minAssoc)
  cpc <- m.ind
  # cat(cpc,",\t",which( !allowed[x,] ),"\n")
  allowed[x, m.ind] <- 0
  minAssoc[m.ind] <- 0

  while (sum(allowed[x, ]) > 0) {
    # try to add one node to the cpc
    for (y in 1:n.nodes)
    {
      if (allowed[x, y] && length(cpc) > 0) {
        # condition on all possible combinations of cpc elements,
        # from smaller to larger
        n <- as.integer(length(cpc))
        s <- sort(node.sizes[cpc], index.return = T)$ix # useful for early stopping
        for (zsize in seq_len(min(c(n, max.fanin[y]))))
        {
          # test for possible early stopping, no indep test can be performed
          if (prod(node.sizes[c(x, y, cpc[s[1:zsize]])]) > n.cases / min.counts) {
            break
          }

          comb <- as.integer(1:zsize)
          while (comb[1] != 0) # check for the end of combinations
          {
            assoc <- g2(data, node.sizes, x, y, chi.th, cpc[comb], min.counts)
            if (assoc < minAssoc[y]) {
              minAssoc[y] <- assoc
            }
            if (assoc == 0) # do not try with y anymore
              {
                allowed[x, y] <- 0
                break
              }
            comb <- .Call("bnstruct_next_comb", comb, n, PACKAGE = "bnstruct")
          }
          if (allowed[x, y] == 0) { # came out from the break
            break
          }
        }
      }
    }
    # cat(minAssoc,"\n")
    m <- max(minAssoc)
    if (m == 0) {
      break
    }

    m.ind <- which.max(minAssoc)
    minAssoc[m.ind] <- 0
    cpc <- union(cpc, m.ind)
    # cat(cpc,",\t",which( !allowed[x,] ),"\n")
    allowed[x, m.ind] <- 0
  }

  cpc.vec <- rep(0, n.nodes)
  cpc.vec[cpc] <- 1
  return(cpc.vec)
}

mmpc.bwd <- function(data, node.sizes, cpc.vec, x, chi.th, min.counts, max.fanin) {
  #  cat("\n",x,":\n");
  #  cat(which(cpc.vec!=0),"\n");
  # not worth if cpc has less than 2 elmts
  if (sum(cpc.vec) < 2) {
    return(cpc.vec)
  }

  n.nodes <- ncol(data)
  n.cases <- nrow(data)

  for (y in 1:n.nodes)
  {
    if (cpc.vec[y] != 0) {
      # condition on all possible combinations of cpc elements,
      # from smaller to larger
      cpc <- setdiff(which(cpc.vec > 0), y)
      n <- as.integer(length(cpc))
      s <- sort(node.sizes[cpc], index.return = T)$ix # useful for early stopping
      for (zsize in seq_len(min(c(n, max.fanin[y]))))
      {
        # test for possible early stopping, no indep test can be performed
        if (prod(node.sizes[c(x, y, cpc[s[1:zsize]])]) > n.cases / min.counts) {
          break
        }

        comb <- as.integer(1:zsize)
        while (comb[1] != 0) # check for the end of combinations
        {
          assoc <- g2(data, node.sizes, x, y, chi.th, cpc[comb], min.counts)
          if (assoc == 0) # do not try with y anymore
            {
              cpc.vec[y] <- 0
              #            cat(which(cpc.vec!=0),"\n");
              break
            }
          comb <- .Call("bnstruct_next_comb", comb, n, PACKAGE = "bnstruct")
        }
        if (cpc.vec[y] == 0) { # came out from the break
          break
        }
      }
    }
  }

  return(cpc.vec)
}

g2 <- function(data, sizes, x, y, chi.th = 0.05, z = c(), min.counts = 5) {
  # if less than min.counts (default 5) counts per cell on average, cannot exclude dependence
  if (dim(data)[1] / prod(sizes[c(x, y, z)]) < min.counts) {
    return(chi.th)
  }

  # tab <- compute.counts( data[,c(x,y,z)], sizes[c(x,y,z)] )
  #   tab <- .Call("compute_counts", data = data[,c(x,y,z)], node_sizes = sizes[c(x,y,z)],
  #                PACKAGE = "bnstruct" )
  #
  #   dim(tab) <- c(sizes[x], sizes[y], prod(sizes[z]))
  #   sz <- array( rep(apply(tab, 3, sum), each=sizes[x]*sizes[y]), dim(tab) )
  #   syz <- array( rep(apply(tab, c(2,3), sum), each=sizes[x]), dim(tab) )
  #   sxz <- aperm( array(rep(apply(tab, c(1,3), sum), each=sizes[y]),
  #                       c(dim(tab)[2],dim(tab)[1],dim(tab)[3])), c(2,1,3) )
  #
  #   s <- tab * log((tab*sz)/(syz*sxz)) # there can be NaNs for cells with count 0
  #   s[is.na(s)] <- 0
  #
  #   print(tab)
  #   print(sz)
  #   print(syz)
  #   print(sxz)
  #   print(s)
  #
  #   df <- sum(tab!=0) + sum(sz!=0)/(sizes[x]*sizes[y]) -
  #     (sum(sxz!=0)/sizes[y] + sum(syz!=0)/sizes[x])
  #
  #   # sanity check
  #   if( df < 0 )
  #   {
  #     cat("Negative df\n")
  #     return( chi.th )
  #   }
  #
  #   print(c(2*sum(s),df))
  #
  #   return( max(pchisq(2*sum(s),df) - 1+chi.th, 0) )

  stat <- .Call("bnstruct_g2_stat",
    data = data[, c(x, y, z)], node_sizes = sizes[c(x, y, z)],
    PACKAGE = "bnstruct"
  )

  # print(stat)

  # sanity check
  if (stat[2] < 0) {
    # cat("Negative df\n")
    return(chi.th)
  }

  return(max(pchisq(stat[1], stat[2]) - 1 + chi.th, 0))
}

hc_w_disc <- function(data, node.sizes, scoring.func = 0, cpc, cont.nodes = c(), ess = 1, tabu.tenure = 100,
                      max.parents = length(node.sizes) - 1,
                      init.net = NULL, wm.max = 15, layering = NULL, layer.struct = NULL,
                      mandatory.edges = NULL, L.card = NULL, max.disc.cycles = 10, approx.parents = 0, n.pre.binning = 0, greedy.heuristic = 0) {
  if (is.null(L.card)) {
    L.card <- -1L
  } else {
    L.card <- as.integer(L.card)
  }
  n.nodes <- ncol(data)
  n.cases <- nrow(data)
  # just to be sure
  storage.mode(node.sizes) <- "integer"
  storage.mode(approx.parents) <- "integer"
  storage.mode(max.disc.cycles) <- "integer"
  storage.mode(cont.nodes) <- "integer"
  storage.mode(n.pre.binning) <- "integer"
  storage.mode(greedy.heuristic) <- "integer"
  # store original data
  original.data <- data
  # Initialize lgamma cache vector
  lgamma_max <- min(n.cases + 100, 1e6)
  lgamma_cache_vec <- lgamma(c(0:(lgamma_max)))
  # quantize data of continuous nodes
  levels <- rep(0, n.nodes)
  levels[cont.nodes] <- node.sizes[cont.nodes]

  # data <- quantize.with.na.matrix( data, levels )
  # data <- quantize.matrix( data, levels )
  curr.discretization <- rep(list(NA), n.nodes)
  out.data <- quantize.matrix(data, levels)
  data <- out.data$quant
  curr.discretization[cont.nodes] <- lapply(out.data$quantiles[cont.nodes], function(x) x[2:(length(x) - 1)])
  # quantiles(bn) <- out.data$quantiles
  # quantiles(dataset) <- out.data$quantiles
  storage.mode(data) <- "integer"
  storage.mode(original.data) <- "double"
  # apply mandatory edges (if present)
  m.edges <- matrix(0L, n.nodes, n.nodes)
  if (!is.null(mandatory.edges)) {
    m.edges <- m.edges | mandatory.edges
  }
  storage.mode(m.edges) <- "integer"

  # start with init.net if not NULL, otherwise with empty matrix
  if (!is.null(init.net)) {
    curr.g <- init.net
    # just to be sure
    storage.mode(curr.g) <- "integer"
    # add init edges out of cpc
    cpc <- cpc | init.net | t(init.net)
  } else {
    curr.g <- matrix(0L, n.nodes, n.nodes)
  } # integers!
  curr.g <- curr.g | m.edges

  # apply layering
  if (!is.null(layer.struct) && is.null(layering)) {
    stop("layer.struct provided without layering.\n")
  }
  if (!is.null(layering) && length(unique(layering)) > 1) {
    n.layers <- length(unique(layering))

    if (is.null(layer.struct)) {
      layer.struct <- matrix(1L, n.layers, n.layers)
      layer.struct[lower.tri(layer.struct)] <- 0
      layer.struct[1, 1] <- 0
    }

    layers <- matrix(1L, n.nodes, n.nodes)
    for (i in 1:n.layers) {
      for (j in 1:n.layers) {
        layers[which(layering == i), which(layering == j)] <- layer.struct[i, j]
      }
    }
    diag(layers) <- 0

    # keep only edges allowed by both the CPC / initial network and the layering
    cpc <- cpc & layers

    if (sum(m.edges | t(m.edges)) > 0L && sum(curr.g & layers) == 0L) {
      stop("the mandatory edges provided are inconsistent with the given layering.\n")
    }
  }
  # end apply layering

  # parents and children lists

  out.par_child <- adj_to_parents_children(curr.g)
  curr.parents <- out.par_child$parents
  curr.children <- out.par_child$children
  rm(out.par_child)

  curr.score.nodes <- array(0, n.nodes)
  for (i in 1L:n.nodes) {
    curr.score.nodes[i] <- .Call("bnstruct_score_node", data, node.sizes, i - 1L, which(curr.g[, i] != 0) - 1L,
      scoring.func, ess,
      PACKAGE = "bnstruct"
    )
  }

  # global best solution
  global.best.g <- curr.g
  global.best.score <- sum(curr.score.nodes)
  global.best.discretization <- curr.discretization

  # tabu list
  tabu <- array(0L, c(n.nodes, n.nodes, tabu.tenure))
  tabu.pt <- 1

  # worsening moves
  wm.count <- 0
  # wm.max <- 15
  while (wm.count < wm.max) {
    next.score.diff <- rep(-Inf, n.nodes)
    next.pert <- rep(-1L, n.nodes)
    # print(tabu[,,1:10])
    # try all possible perturbations
    for (node in 1L:n.nodes)
    {
      for (par in 1L:n.nodes)
      {
        if (cpc[par, node] == 1L && (m.edges[par, node] == 0 && m.edges[node, par] == 0)) {
          next.g <- curr.g
          s.diff <- -Inf

          if (curr.g[par, node] == 1L) # edge removal
            {
              next.g[par, node] <- 0L
              if (.Call("bnstruct_is_acyclic", next.g, PACKAGE = "bnstruct") & !.Call("bnstruct_in_tabu", next.g, tabu, PACKAGE = "bnstruct")) {
                #               cat(node.sizes,'\t',node-1L,'\t',which(next.g[,node]!=0)-1L,'\n')
                s.diff <- .Call("bnstruct_score_node", data, node.sizes, node - 1L, which(next.g[, node] != 0) - 1L,
                  scoring.func, ess,
                  PACKAGE = "bnstruct"
                ) - curr.score.nodes[node]
              }
            }
          # edge addition
          # check also if there is room for one more parent
          else if (curr.g[node, par] == 0L && sum(curr.g[, node]) < max.parents) {
            next.g[par, node] <- 1L
            # print(c(node,par,.Call("is_acyclic", next.g, PACKAGE = "bnstruct"),!.Call("in_tabu", next.g, tabu, PACKAGE = "bnstruct")))
            if (.Call("bnstruct_is_acyclic", next.g, PACKAGE = "bnstruct") & !.Call("bnstruct_in_tabu", next.g, tabu, PACKAGE = "bnstruct")) {
              # print("here\n");

              #               cat(node.sizes,'\t',node-1L,'\t',which(next.g[,node]!=0)-1L,'\n')
              s.diff <- .Call("bnstruct_score_node", data, node.sizes, node - 1L, which(next.g[, node] != 0) - 1L,
                scoring.func, ess,
                PACKAGE = "bnstruct"
              ) - curr.score.nodes[node]
            }
          }
          # edge reversal
          # check also if there is room for one more parent
          else if (sum(curr.g[, node]) < max.parents) {
            next.g[par, node] <- 1L
            next.g[node, par] <- 0L
            if (.Call("bnstruct_is_acyclic", next.g, PACKAGE = "bnstruct") & !.Call("bnstruct_in_tabu", next.g, tabu, PACKAGE = "bnstruct")) {
              #               cat(node.sizes,'\t',node-1L,'\t',which(next.g[,node]!=0)-1L,'\t',par-1L,'\t',which(next.g[,par]!=0)-1L,'\n')
              s.diff <- .Call("bnstruct_score_node", data, node.sizes, node - 1L, which(next.g[, node] != 0) - 1L,
                scoring.func, ess,
                PACKAGE = "bnstruct"
              ) +
                .Call("bnstruct_score_node", data, node.sizes, par - 1L, which(next.g[, par] != 0) - 1L,
                  scoring.func, ess,
                  PACKAGE = "bnstruct"
                ) -
                (curr.score.nodes[node] + curr.score.nodes[par])
            }
          }
          # test for local improvement
          if (s.diff > next.score.diff[node]) {
            next.score.diff[node] <- s.diff
            next.pert[node] <- par
          }
        }
      }
    }

    best.node <- which.max(next.score.diff)
    # browser()
    # no possible improvements given the tabu list
    if (next.score.diff[best.node] == -Inf) {
      break
    }

    # save old graph
    old.g <- curr.g

    # update current graph and scores
    curr.g[next.pert[best.node], best.node] <- 1L - curr.g[next.pert[best.node], best.node] # flip
    if (curr.g[next.pert[best.node], best.node] == 1L) {
      curr.parents[[best.node]] <- c(curr.parents[[best.node]], next.pert[best.node])
      curr.children[[next.pert[best.node]]] <- c(curr.children[[next.pert[best.node]]], best.node)
    } else {
      curr.parents[[best.node]] <- setdiff(curr.parents[[best.node]], next.pert[best.node])
      curr.children[[next.pert[best.node]]] <- setdiff(curr.children[[next.pert[best.node]]], best.node)
    }
    if (curr.g[best.node, next.pert[best.node]] == 1L) # need to reverse
      {
        curr.g[best.node, next.pert[best.node]] <- 0L
        curr.parents[[next.pert[best.node]]] <- setdiff(curr.parents[[next.pert[best.node]]], best.node)
        curr.children[[best.node]] <- setdiff(curr.children[[best.node]], next.pert[best.node])
        # cat(node.sizes,'\t',best.node-1L,'\t',which(curr.g[,best.node]!=0)-1L,'\t',next.pert[best.node]-1L,'\t',which(curr.g[,next.pert[best.node]]!=0)-1L,'\n')
        curr.score.nodes[best.node] <- .Call("bnstruct_score_node", data, node.sizes, best.node - 1L,
          which(curr.g[, best.node] != 0) - 1L, scoring.func, ess,
          PACKAGE = "bnstruct"
        )
        curr.score.nodes[next.pert[best.node]] <- .Call("bnstruct_score_node", data, node.sizes, next.pert[best.node] - 1L,
          which(curr.g[, next.pert[best.node]] != 0) - 1L, scoring.func, ess,
          PACKAGE = "bnstruct"
        )
      } else {
      curr.score.nodes[best.node] <- curr.score.nodes[best.node] + next.score.diff[best.node]
    }

    # detect continuous nodes affected by the perturbation and run discretization routine
    affected <- detect_affected_continuous_vars(old.g, curr.g, cont.nodes)
    n.cont_affected <- length(affected)
    if (n.cont_affected > 0) {
      top_order <- topological.sort(curr.g)
      affected_order <- as.integer(rev(top_order[top_order %in% affected]))

      out.disc <- .Call("bnstruct_dvbn_discretize_all", original.data[, affected_order, drop = F], data, n.nodes, n.cases, n.cont_affected, affected_order - 1L, node.sizes, lapply(curr.parents, function(x) as.integer(x - 1L)), lapply(curr.children, function(x) as.integer(x - 1L)),
        max.disc.cycles, approx.parents, L.card, n.pre.binning, greedy.heuristic, lgamma_cache_vec,
        PACKAGE = "bnstruct"
      )
      curr.discretization[affected_order] <- out.disc
      data[, affected] <- apply_discretization(original.data[, affected, drop = F], curr.discretization[affected])
      node.sizes[affected] <- sapply(curr.discretization[affected], length) + 1L

      # recalculate scores for affected nodes
      nodes_to_rescore <- unique(unlist(c(curr.children[affected], affected)))
      for (node_idx in nodes_to_rescore) {
        curr.score.nodes[node_idx] <- .Call("bnstruct_score_node", data, node.sizes, node_idx - 1L,
          which(curr.g[, node_idx] != 0) - 1L, scoring.func, ess,
          PACKAGE = "bnstruct"
        )
      }
    }
    # print(c(tabu.pt,best.node,next.pert[best.node],sum(curr.score.nodes)))

    if (global.best.score < sum(curr.score.nodes)) # check for global best
      {
        wm.count <- 0
        global.best.g <- curr.g
        global.best.score <- sum(curr.score.nodes)
        global.best.discretization <- curr.discretization
      } else {
      wm.count <- wm.count + 1
    }

    # update tabu list
    tabu[, , tabu.pt] <- curr.g
    tabu.pt <- (tabu.pt) %% tabu.tenure + 1
    # print(curr.g)
  }

  return(list(dag = global.best.g, discretization = global.best.discretization))
}

hc_block_disc <- function(data, node.sizes, scoring.func = 0, cpc, cont.nodes = c(), ess = 1, tabu.tenure = 100,
                          max.parents = length(node.sizes) - 1,
                          init.net = NULL, wm.max = 15, layering = NULL, layer.struct = NULL,
                          mandatory.edges = NULL, L.card = NULL, max.disc.cycles = 10, approx.parents = 0,
                          n.pre.binning = 0, greedy.heuristic = 0, max.outer.iter = 10,
                          disc.tol = 1e-8, outer.trace = FALSE) {
  if (is.null(L.card)) {
    L.card <- -1L
  } else {
    L.card <- as.integer(L.card)
  }
  n.nodes <- ncol(data)
  n.cases <- nrow(data)

  # just to be sure
  storage.mode(node.sizes) <- "integer"
  storage.mode(approx.parents) <- "integer"
  storage.mode(max.disc.cycles) <- "integer"
  storage.mode(cont.nodes) <- "integer"
  storage.mode(n.pre.binning) <- "integer"
  storage.mode(greedy.heuristic) <- "integer"
  storage.mode(max.outer.iter) <- "integer"

  original.data <- data
  storage.mode(original.data) <- "double"

  # Initialize lgamma cache vector
  lgamma_max <- min(n.cases + 100, 1e6)
  lgamma_cache_vec <- lgamma(c(0:(lgamma_max)))

  # Initial quantization from provided cardinalities
  levels <- rep(0, n.nodes)
  levels[cont.nodes] <- node.sizes[cont.nodes]
  out.data <- quantize.matrix(data, levels)
  data <- out.data$quant
  storage.mode(data) <- "integer"

  curr.discretization <- rep(list(NA), n.nodes)
  if (length(cont.nodes) > 0L) {
    curr.discretization[cont.nodes] <- lapply(out.data$quantiles[cont.nodes], function(x) x[2:(length(x) - 1L)])
    node.sizes[cont.nodes] <- as.integer(sapply(curr.discretization[cont.nodes], length) + 1L)
  }

  if (!is.null(init.net)) {
    curr.g <- init.net
  } else {
    curr.g <- matrix(0L, n.nodes, n.nodes)
  }
  storage.mode(curr.g) <- "integer"

  same_discretization <- function(old.disc, new.disc, nodes, tol) {
    if (length(nodes) == 0L) {
      return(TRUE)
    }
    for (v in nodes) {
      old.v <- old.disc[[v]]
      new.v <- new.disc[[v]]
      if (all(is.na(old.v)) && all(is.na(new.v))) {
        next
      }
      if (length(old.v) != length(new.v)) {
        return(FALSE)
      }
      if (length(old.v) > 0L && any(abs(old.v - new.v) > tol)) {
        return(FALSE)
      }
    }
    TRUE
  }

  max_discretization_diff <- function(old.disc, new.disc, nodes) {
    if (length(nodes) == 0L) {
      return(0)
    }
    max.diff <- 0
    for (v in nodes) {
      old.v <- old.disc[[v]]
      new.v <- new.disc[[v]]
      if (all(is.na(old.v)) && all(is.na(new.v))) {
        next
      }
      if (length(old.v) != length(new.v)) {
        return(Inf)
      }
      if (length(old.v) > 0L) {
        max.diff <- max(max.diff, max(abs(old.v - new.v)))
      }
    }
    max.diff
  }

  total_score <- function(curr.data, curr.sizes, curr.g, curr.scoring.func, curr.ess) {
    node.scores <- numeric(ncol(curr.data))
    for (i in seq_len(ncol(curr.data))) {
      node.scores[i] <- .Call(
        "bnstruct_score_node",
        curr.data,
        curr.sizes,
        i - 1L,
        which(curr.g[, i] != 0) - 1L,
        curr.scoring.func,
        curr.ess,
        PACKAGE = "bnstruct"
      )
    }
    sum(node.scores)
  }

  iter <- 0L
  exit.reason <- ""
  history <- data.frame(
    iter = integer(0),
    score = numeric(0),
    edge.count = integer(0),
    dag.diff = integer(0),
    disc.max.diff = numeric(0),
    dag.stable = logical(0),
    disc.stable = logical(0),
    stringsAsFactors = FALSE
  )

  repeat {
    iter <- iter + 1L
    prev.g <- curr.g
    prev.discretization <- curr.discretization

    # Full structure learning phase on fixed discretized data
    curr.g <- hc(
      data = data,
      node.sizes = node.sizes,
      scoring.func = scoring.func,
      cpc = cpc,
      cont.nodes = integer(0),
      ess = ess,
      tabu.tenure = tabu.tenure,
      max.parents = max.parents,
      init.net = curr.g,
      wm.max = wm.max,
      layering = layering,
      layer.struct = layer.struct,
      mandatory.edges = mandatory.edges
    )
    storage.mode(curr.g) <- "integer"

    # Full re-discretization phase
    if (length(cont.nodes) > 0L) {
      out.par_child <- adj_to_parents_children(curr.g)
      curr.parents <- out.par_child$parents
      curr.children <- out.par_child$children

      top_order <- topological.sort(curr.g)
      affected_order <- as.integer(rev(top_order[top_order %in% cont.nodes]))

      out.disc <- .Call(
        "bnstruct_dvbn_discretize_all",
        original.data[, affected_order, drop = F],
        data,
        n.nodes,
        n.cases,
        length(affected_order),
        affected_order - 1L,
        node.sizes,
        lapply(curr.parents, function(x) as.integer(x - 1L)),
        lapply(curr.children, function(x) as.integer(x - 1L)),
        max.disc.cycles,
        approx.parents,
        L.card,
        n.pre.binning,
        greedy.heuristic,
        lgamma_cache_vec,
        PACKAGE = "bnstruct"
      )

      curr.discretization[affected_order] <- out.disc
      data[, cont.nodes] <- apply_discretization(original.data[, cont.nodes, drop = F], curr.discretization[cont.nodes])
      node.sizes[cont.nodes] <- as.integer(sapply(curr.discretization[cont.nodes], length) + 1L)
    }

    dag.stable <- identical(prev.g, curr.g)
    dag.diff <- sum(prev.g != curr.g)
    disc.max.diff <- max_discretization_diff(prev.discretization, curr.discretization, cont.nodes)
    disc.stable <- same_discretization(prev.discretization, curr.discretization, cont.nodes, disc.tol)
    curr.score <- total_score(data, node.sizes, curr.g, scoring.func, ess)
    history <- rbind(
      history,
      data.frame(
        iter = iter,
        score = curr.score,
        edge.count = sum(curr.g != 0L),
        dag.diff = dag.diff,
        disc.max.diff = disc.max.diff,
        dag.stable = dag.stable,
        disc.stable = disc.stable,
        stringsAsFactors = FALSE
      )
    )

    if (isTRUE(outer.trace)) {
      message(sprintf(
        "hc_block_disc outer iteration %d: score=%.6f edges=%d dag.diff=%d disc.max.diff=%g dag.stable=%s disc.stable=%s",
        iter,
        curr.score,
        sum(curr.g != 0L),
        dag.diff,
        disc.max.diff,
        dag.stable,
        disc.stable
      ))
    }

    if (dag.stable && disc.stable) {
      exit.reason <- "converged"
      break
    }
    if (iter >= max.outer.iter) {
      exit.reason <- "max.outer.iter"
      break
    }
  }

  if (isTRUE(outer.trace)) {
    message(sprintf("hc_block_disc exit: reason=%s iterations=%d", exit.reason, iter))
  }

  list(
    dag = curr.g,
    discretization = curr.discretization,
    diagnostics = list(
      exit.reason = exit.reason,
      iterations = iter,
      history = history
    )
  )
}

K2 <- function(data, node.sizes, scoring.func = 0, cpc, cont.nodes = c(), ess = 1,
               max.parents = length(node.sizes) - 1,
               init.net = NULL, wm.max = 15, layering = NULL, layer.struct = NULL,
               mandatory.edges = NULL, n.starts = 1) {
  n.nodes <- ncol(data)
  n.cases <- nrow(data)
  # just to be sure
  storage.mode(node.sizes) <- "integer"

  # quantize data of continuous nodes
  levels <- rep(0, n.nodes)
  levels[cont.nodes] <- node.sizes[cont.nodes]

  # data <- quantize.with.na.matrix( data, levels )
  # data <- quantize.matrix( data, levels )
  out.data <- quantize.matrix(data, levels)
  data <- out.data$quant

  # apply mandatory edges (if present)
  m.edges <- matrix(0L, n.nodes, n.nodes)
  if (!is.null(mandatory.edges)) {
    m.edges <- m.edges | mandatory.edges
  }
  storage.mode(m.edges) <- "integer"

  # start with init.net if not NULL, otherwise with empty matrix
  if (!is.null(init.net)) {
    curr.g <- init.net
    # just to be sure
    storage.mode(curr.g) <- "integer"
    # add init edges out of cpc
    cpc <- cpc | init.net | t(init.net)
  } else {
    curr.g <- matrix(0L, n.nodes, n.nodes)
  } # integers!
  curr.g <- curr.g | m.edges

  # apply layering
  if (!is.null(layer.struct) && is.null(layering)) {
    stop("layer.struct provided without layering.\n")
  }
  if (!is.null(layering) && length(unique(layering)) > 1) {
    n.layers <- length(unique(layering))

    if (is.null(layer.struct)) {
      layer.struct <- matrix(1L, n.layers, n.layers)
      layer.struct[lower.tri(layer.struct)] <- 0
      layer.struct[1, 1] <- 0
    }

    layers <- matrix(1L, n.nodes, n.nodes)
    for (i in 1:n.layers) {
      for (j in 1:n.layers) {
        layers[which(layering == i), which(layering == j)] <- layer.struct[i, j]
      }
    }
    diag(layers) <- 0

    # keep only edges allowed by both the CPC / initial network and the layering
    cpc <- cpc & layers
    if (sum(m.edges | t(m.edges)) > 0L && sum(curr.g & layers) == 0L) {
      stop("the mandatory edges provided are inconsistent with the given layering.\n")
    }
  }
  # end apply layering

  # generate random orders compatible with the partially ordered graph
  p.order_graph <- matrix(0L, n.nodes, n.nodes)
  for (i in 1:n.nodes) {
    for (j in 1:n.nodes) {
      if (cpc[i, j] != cpc[j, i]) {
        if (cpc[i, j] == 1L) {
          p.order_graph[i, j] <- 1L
        } else {
          p.order_graph[j, i] <- 1L
        }
      }
    }
  }
  p.order_graph <- p.order_graph | m.edges
  # generate random orders consistent with the partially ordered graph
  orders <- lapply(rep(list(p.order_graph), n.starts * 10), rand.topological.sort)
  # Keep only n.starts unique orders if available
  orders <- unique(orders)
  if (length(orders) > n.starts) {
    orders <- orders[1:n.starts]
  } else {
    message(sprintf("Only %d unique orders could be generated. Requested n.starts were %d. Running with available orders.", length(orders), n.starts))
  }
  # end generate random orders

  # initialize scores
  curr.score.nodes <- array(0, n.nodes)
  for (i in 1L:n.nodes) {
    curr.score.nodes[i] <- .Call("bnstruct_score_node", data, node.sizes, i - 1L, which(curr.g[, i] != 0) - 1L,
      scoring.func, ess,
      PACKAGE = "bnstruct"
    )
  }

  # initialize global best solution
  global.best.g <- curr.g
  global.best.score <- sum(curr.score.nodes)

  curr.g.preK2 <- curr.g
  curr.score.nodes.preK2 <- curr.score.nodes

  # K2 search for each order
  for (start in seq_along(orders)) {
    order <- orders[[start]]
    # reset to pre-K2 state
    curr.g <- curr.g.preK2
    curr.score.nodes <- curr.score.nodes.preK2
    # Main K2 loop
    for (i in 1L:n.nodes) {
      node <- order[i]
      while (TRUE) {
        best.par <- -1L
        best.s.diff <- 0L
        for (par in order[seq_len(i - 1L)]) {
          if (curr.g[par, node] == 0 && cpc[par, node] == 1L && (m.edges[par, node] == 0 && m.edges[node, par] == 0) && sum(curr.g[, node]) < max.parents) {
            next.g <- curr.g
            next.g[par, node] <- 1L
            s.diff <- .Call("bnstruct_score_node", data, node.sizes, node - 1L, which(next.g[, node] != 0) - 1L,
              scoring.func, ess,
              PACKAGE = "bnstruct"
            ) - curr.score.nodes[node]
            if (s.diff > best.s.diff) {
              best.s.diff <- s.diff
              best.par <- par
            }
          }
        }
        if (best.par != -1L) {
          curr.g[best.par, node] <- 1L
          curr.score.nodes[node] <- curr.score.nodes[node] + best.s.diff
        } else {
          break
        }
      }
    } # end K2 for the given order

    # check for global best
    if (global.best.score < sum(curr.score.nodes)) # check for global best
      {
        global.best.g <- curr.g
        global.best.score <- sum(curr.score.nodes)
      }
  } # end loop over orders

  return(global.best.g)
}

K2_w_disc <- function(data, node.sizes, scoring.func = 0, cpc, cont.nodes = c(), ess = 1,
                      max.parents = length(node.sizes) - 1,
                      init.net = NULL, wm.max = 15, layering = NULL, layer.struct = NULL,
                      mandatory.edges = NULL, n.starts = 1, max.disc.cycles = 5, approx.parents = 0, L.card = NULL, n.pre.binning = 0, greedy.heuristic = 1L) {
  if (is.null(L.card)) {
    L.card <- -1L
  } else {
    L.card <- as.integer(L.card)
  }
  storage.mode(approx.parents) <- "integer"
  storage.mode(max.disc.cycles) <- "integer"
  storage.mode(cont.nodes) <- "integer"
  storage.mode(n.pre.binning) <- "integer"
  storage.mode(greedy.heuristic) <- "integer"
  n.nodes <- ncol(data)
  n.cases <- nrow(data)
  # just to be sure
  storage.mode(node.sizes) <- "integer"
  # Initialize lgamma cache vector
  lgamma_max <- min(n.cases + 100, 1e6)
  lgamma_cache_vec <- lgamma(c(0:(lgamma_max)))

  # quantize data of continuous nodes
  # Store original continuous data before quantization
  cont_data <- data[, cont.nodes, drop = F]
  storage.mode(cont_data) <- "double"

  levels <- rep(0, n.nodes)
  levels[cont.nodes] <- node.sizes[cont.nodes]

  out.data <- quantize.matrix(data, levels)
  data <- out.data$quant
  storage.mode(data) <- "integer"

  # apply mandatory edges (if present)
  m.edges <- matrix(0L, n.nodes, n.nodes)
  if (!is.null(mandatory.edges)) {
    m.edges <- m.edges | mandatory.edges
  }
  storage.mode(m.edges) <- "integer"

  # start with init.net if not NULL, otherwise with empty matrix
  if (!is.null(init.net)) {
    curr.g <- init.net
    # just to be sure
    storage.mode(curr.g) <- "integer"
    # add init edges out of cpc
    cpc <- cpc | init.net | t(init.net)
  } else {
    curr.g <- matrix(0L, n.nodes, n.nodes)
  } # integers!
  curr.g <- curr.g | m.edges

  # apply layering
  if (!is.null(layer.struct) && is.null(layering)) {
    stop("layer.struct provided without layering.\n")
  }
  if (!is.null(layering) && length(unique(layering)) > 1) {
    n.layers <- length(unique(layering))

    if (is.null(layer.struct)) {
      layer.struct <- matrix(1L, n.layers, n.layers)
      layer.struct[lower.tri(layer.struct)] <- 0
      layer.struct[1, 1] <- 0
    }

    layers <- matrix(1L, n.nodes, n.nodes)
    for (i in 1:n.layers) {
      for (j in 1:n.layers) {
        layers[which(layering == i), which(layering == j)] <- layer.struct[i, j]
      }
    }
    diag(layers) <- 0

    # keep only edges allowed by both the CPC / initial network and the layering
    cpc <- cpc & layers
    if (sum(m.edges | t(m.edges)) > 0L && sum(curr.g & layers) == 0L) {
      stop("the mandatory edges provided are inconsistent with the given layering.\n")
    }
  }
  # end apply layering

  # generate random orders compatible with the partially ordered graph
  p.order_graph <- matrix(0L, n.nodes, n.nodes)
  for (i in 1:n.nodes) {
    for (j in 1:n.nodes) {
      if (cpc[i, j] != cpc[j, i]) {
        if (cpc[i, j] == 1L) {
          p.order_graph[i, j] <- 1L
        } else {
          p.order_graph[j, i] <- 1L
        }
      }
    }
  }
  p.order_graph <- p.order_graph | m.edges
  # generate random orders consistent with the partially ordered graph
  orders <- lapply(rep(list(p.order_graph), n.starts * 10), rand.topological.sort)
  # Keep only n.starts unique orders if available
  orders <- unique(orders)
  if (length(orders) > n.starts) {
    orders <- orders[1:n.starts]
  } else {
    message(sprintf("Only %d unique orders could be generated. Requested n.starts were %d. Running with available orders.", length(orders), n.starts))
  }
  # end generate random orders

  # initialize scores
  curr.score.nodes <- array(0, n.nodes)
  for (i in 1L:n.nodes) {
    curr.score.nodes[i] <- .Call("bnstruct_score_node", data, node.sizes, i - 1L, which(curr.g[, i] != 0) - 1L,
      scoring.func, ess,
      PACKAGE = "bnstruct"
    )
  }

  # initialize global best solution
  global.best.g <- curr.g
  global.best.score <- sum(curr.score.nodes)
  global.best.disc <- out.data$quantiles

  curr.g.preK2 <- curr.g
  curr.score.nodes.preK2 <- curr.score.nodes
  out.par_child <- adj_to_parents_children(curr.g)

  # K2 search for each order
  for (start in seq_along(orders)) {
    order <- orders[[start]]

    # Reset to pre-K2 state
    curr.g <- curr.g.preK2
    curr.score.nodes <- curr.score.nodes.preK2
    curr.parents <- out.par_child$parents
    curr.children <- out.par_child$children
    # Reset entire data matrix
    data <- out.data$quant
    storage.mode(data) <- "integer"
    curr.disc <- out.data$quantiles
    curr.node.sizes <- as.integer(node.sizes) # Explicit copy and ensure integer type
    # Main K2 loop
    for (i in 1L:n.nodes) {
      node <- order[i]
      while (TRUE) {
        best.par <- -1L
        best.s.diff <- 0L # Initialize to 0 so we only accept improvements

        # --- STEP 1: SCAN CANDIDATES (Fixed Discretization) ---
        # We assume the current 'data' is the "truth" for this comparison
        for (par in order[seq_len(i - 1L)]) {
          # Check validity (CPC, mandatory, max parents)
          if (curr.g[par, node] == 0 && cpc[par, node] == 1L &&
            (m.edges[par, node] == 0 && m.edges[node, par] == 0) &&
            sum(curr.g[, node]) < max.parents) {
            # Temporarily form the parent set for scoring
            current_pars <- curr.parents[[node]]
            temp_pars <- c(current_pars, par)

            # Compute score with Candidate Parent using CURRENT data
            s.new <- .Call("bnstruct_score_node", data, curr.node.sizes, node - 1L,
              as.integer(temp_pars - 1L), scoring.func, ess,
              PACKAGE = "bnstruct"
            )

            s.diff <- s.new - curr.score.nodes[node]

            if (s.diff > best.s.diff) {
              best.s.diff <- s.diff
              best.par <- par
            }
          }
        }

        # --- STEP 2: COMMIT & RE-DISCRETIZE ---
        if (best.par != -1L) {
          # 1. Update Structure
          prev.g <- curr.g # Save previous graph for 'detect_affected'
          curr.g[best.par, node] <- 1L
          curr.parents[[node]] <- c(curr.parents[[node]], best.par)
          curr.children[[best.par]] <- c(curr.children[[best.par]], node)

          # 2. Alternating Step: Re-Discretize
          # We check if this structural change affects any continuous variables
          affected <- detect_affected_continuous_vars(prev.g, curr.g, cont.nodes)
          n.cont_affected <- length(affected)

          if (n.cont_affected > 0) {
            top_order <- topological.sort(curr.g)
            affected_order <- as.integer(rev(top_order[top_order %in% affected]))
            affected_order_cont_data <- match(affected_order, cont.nodes)

            # CALL C DISCRETIZER
            out.disc <- .Call(
              "bnstruct_dvbn_discretize_all",
              cont_data[, affected_order_cont_data, drop = F],
              data, n.nodes, n.cases, n.cont_affected, affected_order - 1L, curr.node.sizes,
              lapply(curr.parents, function(x) as.integer(x - 1L)),
              lapply(curr.children, function(x) as.integer(x - 1L)),
              max.disc.cycles, approx.parents, L.card, n.pre.binning,
              greedy.heuristic, lgamma_cache_vec
            )

            # Update Data
            data[, affected_order] <- apply_discretization(cont_data[, affected_order_cont_data, drop = F], out.disc)
            storage.mode(data) <- "integer"
            curr.disc[affected_order] <- out.disc
            # Update Node Sizes
            curr.node.sizes[affected_order] <- sapply(curr.disc[affected_order], length) + 1L
            storage.mode(curr.node.sizes) <- "integer"

            # 3. Update Scores
            # Because discretization changed, we muast re-calculate the score for affected nodes and their children
            children_of_affected <- unlist(curr.children[affected])
            nodes_to_update <- unique(c(affected, children_of_affected, node))
            for (n in nodes_to_update) {
              curr.score.nodes[n] <- .Call("bnstruct_score_node", data, curr.node.sizes, n - 1L,
                as.integer(curr.parents[[n]] - 1L), scoring.func, ess,
                PACKAGE = "bnstruct"
              )
            }
          } else {
            # No discretization change, use incremental score update
            curr.score.nodes[node] <- curr.score.nodes[node] + best.s.diff
          }
        } else {
          break # No parent improves score
        }
      }
    }

    # check for global best
    if (global.best.score < sum(curr.score.nodes)) {
      global.best.g <- curr.g
      global.best.score <- sum(curr.score.nodes)
      global.best.disc <- curr.disc
    }
  }

  return(list(dag = global.best.g, discretization = global.best.disc))
}
