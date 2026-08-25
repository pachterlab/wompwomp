# Pure R port of the NeighborNet circular-ordering algorithm (Bryant and
# Moulton, 2004; Huson and Bryant, 2006), following the implementation in the
# Python 'splitspy' package (David J. Bryant and Daniel H. Huson,
# GPL-licensed), which wompwomp previously called via reticulate.
#
# Nodes are represented as plain integer ids (0 = the linked-list head
# sentinel; 1..n = the original taxa; n+1.. = internal nodes created during
# agglomeration), with per-node fields (nbr/ch1/ch2/nxt/prv/Rx/Sx) stored as
# plain vectors indexed by id+1, mutated in place via an environment. An
# earlier version used one R6 (environment) object per node with `$`-based
# field access; R6 dispatch overhead in the O(n^3) inner loop made it ~15x
# slower than the Python reference for larger n (e.g. ~15s vs ~1s at n=150).
# Plain vector indexing avoids that overhead entirely while preserving the
# exact same algorithm and output.

# NA_integer_ means "no pointer" (used for nbr/ch1/ch2/nxt/prv). Node id 0 is
# the head sentinel and is never itself a valid nbr/ch1/ch2 target, but it can
# legitimately appear as the `prv` of the first real node, so a value of 0 is
# not usable as a "no pointer" sentinel.
.nnet_setup_nodes <- function(n_tax, max_nodes) {
    size <- max_nodes + 1L
    nodes <- new.env(parent = emptyenv())
    nodes$nbr <- rep(NA_integer_, size)
    nodes$ch1 <- rep(NA_integer_, size)
    nodes$ch2 <- rep(NA_integer_, size)
    nodes$nxt <- rep(NA_integer_, size)
    nodes$prv <- rep(NA_integer_, size)
    nodes$Rx <- numeric(size)
    nodes$Sx <- numeric(size)

    head_nxt <- NA_integer_
    for (i in n_tax:1) {
        nodes$nxt[i + 1L] <- head_nxt
        head_nxt <- i
    }
    nodes$nxt[0L + 1L] <- head_nxt

    node <- 0L
    while (!is.na(nodes$nxt[node + 1L])) {
        nxt_node <- nodes$nxt[node + 1L]
        nodes$prv[nxt_node + 1L] <- node
        node <- nxt_node
    }

    # Reusable scratch buffer for .nnet_active_ids(); active-node count never
    # exceeds n_tax, so this is always large enough without reallocating.
    nodes$scratch_ids <- integer(n_tax)

    nodes
}

.nnet_setup_matrix <- function(n, dist_matrix) {
    max_number_of_nodes <- max(3, 3 * n - 5)
    mat <- matrix(0, nrow = max_number_of_nodes, ncol = max_number_of_nodes)
    dm <- as.matrix(dist_matrix)
    mat[1:n, 1:n] <- dm[1:n, 1:n]
    mat
}

.nnet_join2way <- function(nodes, x, y) {
    nodes$nbr[x + 1L] <- y
    nodes$nbr[y + 1L] <- x
    invisible(NULL)
}

# All currently-active node ids (the nxt-chain from the head), as a plain
# vector, so callers can do vectorized matrix reads/writes instead of a
# scalar loop over the linked list.
.nnet_active_ids <- function(nodes) {
    buf <- nodes$scratch_ids
    cnt <- 0L
    p <- nodes$nxt[0L + 1L]
    while (!is.na(p)) {
        cnt <- cnt + 1L
        buf[cnt] <- p
        p <- nodes$nxt[p + 1L]
    }
    buf[seq_len(cnt)]
}

.nnet_join3way <- function(nodes, mat_env, joins_env, x, y, z, num_nodes) {
    u <- num_nodes + 1L
    v <- num_nodes + 2L

    nodes$ch1[u + 1L] <- x
    nodes$ch2[u + 1L] <- y
    nodes$ch1[v + 1L] <- y
    nodes$ch2[v + 1L] <- z

    u_nxt <- nodes$nxt[x + 1L]
    u_prv <- nodes$prv[x + 1L]
    nodes$nxt[u + 1L] <- u_nxt
    nodes$prv[u + 1L] <- u_prv
    if (!is.na(u_nxt)) nodes$prv[u_nxt + 1L] <- u
    if (!is.na(u_prv)) nodes$nxt[u_prv + 1L] <- u

    v_nxt <- nodes$nxt[z + 1L]
    v_prv <- nodes$prv[z + 1L]
    nodes$nxt[v + 1L] <- v_nxt
    nodes$prv[v + 1L] <- v_prv
    if (!is.na(v_nxt)) nodes$prv[v_nxt + 1L] <- v
    if (!is.na(v_prv)) nodes$nxt[v_prv + 1L] <- v

    y_nxt <- nodes$nxt[y + 1L]
    y_prv <- nodes$prv[y + 1L]
    if (!is.na(y_nxt)) nodes$prv[y_nxt + 1L] <- y_prv
    if (!is.na(y_prv)) nodes$nxt[y_prv + 1L] <- y_nxt

    nodes$nbr[u + 1L] <- v
    nodes$nbr[v + 1L] <- u

    # Vectorized distance-matrix row update: collect the active ids with one
    # cheap linked-list traversal (no matrix access), then update all of row
    # u/v in a single vectorized pass instead of one scalar mat[i,j] access
    # per node. This is the dominant cost of the whole algorithm (profiling
    # showed ~50% of total runtime here at n=150), since it's O(active nodes)
    # per join and there are O(n) joins.
    p_vec <- .nnet_active_ids(nodes)
    if (length(p_vec) > 0L) {
        mp_x <- mat_env$m[x, p_vec]
        mp_y <- mat_env$m[y, p_vec]
        mp_z <- mat_env$m[z, p_vec]

        val_u <- (2.0 / 3.0) * mp_x + mp_y / 3.0
        mat_env$m[u, p_vec] <- val_u
        mat_env$m[p_vec, u] <- val_u

        val_v <- (2.0 / 3.0) * mp_z + mp_y / 3.0
        mat_env$m[v, p_vec] <- val_v
        mat_env$m[p_vec, v] <- val_v
    }
    mat_env$m[u, u] <- 0.0
    mat_env$m[v, v] <- 0.0

    joins_env$top <- joins_env$top + 1L
    joins_env$stack[joins_env$top] <- u

    u
}

.nnet_join4way <- function(nodes, mat_env, joins_env, x2, x, y, y2, num_nodes) {
    u <- .nnet_join3way(nodes, mat_env, joins_env, x2, x, y, num_nodes)
    num_nodes <- num_nodes + 2L

    .nnet_join3way(nodes, mat_env, joins_env, u, nodes$nbr[u + 1L], y2, num_nodes)
    num_nodes <- num_nodes + 2L

    num_nodes
}

.nnet_compute_rx <- function(nodes, mat_env, z, c_x, c_y) {
    r_x <- 0.0
    c_x_nbr <- nodes$nbr[c_x + 1L]
    c_y_nbr <- nodes$nbr[c_y + 1L]
    p <- nodes$nxt[0L + 1L]
    while (!is.na(p)) {
        if (p == c_x || (!is.na(c_x_nbr) && p == c_x_nbr) ||
            p == c_y || (!is.na(c_y_nbr) && p == c_y_nbr) ||
            is.na(nodes$nbr[p + 1L])) {
            r_x <- r_x + mat_env$m[z, p]
        } else {
            r_x <- r_x + mat_env$m[z, p] / 2.0
        }
        p <- nodes$nxt[p + 1L]
    }
    r_x
}

.nnet_join_nodes <- function(nodes, mat_env, joins_env, n) {
    num_nodes <- n
    num_active <- n
    num_clusters <- n

    while (num_active > 3) {
        if (num_active == 4 && num_clusters == 2) {
            p <- nodes$nxt[0L + 1L]
            p_nxt <- nodes$nxt[p + 1L]
            p_nbr <- nodes$nbr[p + 1L]
            q <- if (is.na(p_nbr) || p_nxt != p_nbr) p_nxt else nodes$nxt[p_nxt + 1L]
            q_nbr <- nodes$nbr[q + 1L]

            if (mat_env$m[p, q] + mat_env$m[p_nbr, q_nbr] < mat_env$m[p, q_nbr] + mat_env$m[p_nbr, q]) {
                .nnet_join3way(nodes, mat_env, joins_env, p, q, q_nbr, num_nodes)
            } else {
                .nnet_join3way(nodes, mat_env, joins_env, p, q_nbr, q, num_nodes)
            }
            num_nodes <- num_nodes + 2L
            break
        }

        p <- nodes$nxt[0L + 1L]
        while (!is.na(p)) {
            nodes$Sx[p + 1L] <- 0.0
            p <- nodes$nxt[p + 1L]
        }

        p <- nodes$nxt[0L + 1L]
        while (!is.na(p)) {
            p_nbr <- nodes$nbr[p + 1L]
            if (is.na(p_nbr) || p_nbr > p) {
                q <- nodes$nxt[p + 1L]
                while (!is.na(q)) {
                    q_nbr <- nodes$nbr[q + 1L]
                    if (is.na(q_nbr) || (q_nbr > q && q_nbr != p)) {
                        if (is.na(p_nbr) && is.na(q_nbr)) {
                            d_pq <- mat_env$m[p, q]
                        } else if (!is.na(p_nbr) && is.na(q_nbr)) {
                            d_pq <- (mat_env$m[p, q] + mat_env$m[p_nbr, q]) / 2.0
                        } else if (is.na(p_nbr) && !is.na(q_nbr)) {
                            d_pq <- (mat_env$m[p, q] + mat_env$m[p, q_nbr]) / 2.0
                        } else {
                            d_pq <- (mat_env$m[p, q] + mat_env$m[p, q_nbr] +
                                mat_env$m[p_nbr, q] + mat_env$m[p_nbr, q_nbr]) / 4.0
                        }
                        nodes$Sx[p + 1L] <- nodes$Sx[p + 1L] + d_pq
                        if (!is.na(p_nbr)) nodes$Sx[p_nbr + 1L] <- nodes$Sx[p_nbr + 1L] + d_pq
                        nodes$Sx[q + 1L] <- nodes$Sx[q + 1L] + d_pq
                        if (!is.na(q_nbr)) nodes$Sx[q_nbr + 1L] <- nodes$Sx[q_nbr + 1L] + d_pq
                    }
                    q <- nodes$nxt[q + 1L]
                }
            }
            p <- nodes$nxt[p + 1L]
        }

        c_x <- NA_integer_
        c_y <- NA_integer_
        best <- 0.0

        p <- nodes$nxt[0L + 1L]
        while (!is.na(p)) {
            p_nbr <- nodes$nbr[p + 1L]
            if (!is.na(p_nbr) && (p_nbr < p)) {
                p <- nodes$nxt[p + 1L]
                next
            }
            q <- nodes$nxt[0L + 1L]
            while (!is.na(q)) {
                if (q == p) break
                q_nbr <- nodes$nbr[q + 1L]
                if (!is.na(q_nbr) && (q_nbr < q)) {
                    q <- nodes$nxt[q + 1L]
                    next
                }
                if (!is.na(q_nbr) && q_nbr == p) {
                    q <- nodes$nxt[q + 1L]
                    next
                }
                if (is.na(p_nbr) && is.na(q_nbr)) {
                    d_pq <- mat_env$m[p, q]
                } else if (!is.na(p_nbr) && is.na(q_nbr)) {
                    d_pq <- (mat_env$m[p, q] + mat_env$m[p_nbr, q]) / 2.0
                } else if (is.na(p_nbr) && !is.na(q_nbr)) {
                    d_pq <- (mat_env$m[p, q] + mat_env$m[p, q_nbr]) / 2.0
                } else {
                    d_pq <- (mat_env$m[p, q] + mat_env$m[p, q_nbr] +
                        mat_env$m[p_nbr, q] + mat_env$m[p_nbr, q_nbr]) / 4.0
                }

                q_pq <- (num_clusters - 2.0) * d_pq - nodes$Sx[p + 1L] - nodes$Sx[q + 1L]

                if ((is.na(c_x) || (q_pq < best)) && (is.na(p_nbr) || p_nbr != q)) {
                    c_x <- p
                    c_y <- q
                    best <- q_pq
                }
                q <- nodes$nxt[q + 1L]
            }
            p <- nodes$nxt[p + 1L]
        }

        x <- c_x
        y <- c_y

        c_x_nbr <- nodes$nbr[c_x + 1L]
        c_y_nbr <- nodes$nbr[c_y + 1L]

        if (!is.na(c_x_nbr) || !is.na(c_y_nbr)) {
            nodes$Rx[c_x + 1L] <- .nnet_compute_rx(nodes, mat_env, c_x, c_x, c_y)
            if (!is.na(c_x_nbr)) nodes$Rx[c_x_nbr + 1L] <- .nnet_compute_rx(nodes, mat_env, c_x_nbr, c_x, c_y)
            nodes$Rx[c_y + 1L] <- .nnet_compute_rx(nodes, mat_env, c_y, c_x, c_y)
            if (!is.na(c_y_nbr)) nodes$Rx[c_y_nbr + 1L] <- .nnet_compute_rx(nodes, mat_env, c_y_nbr, c_x, c_y)
        }

        m <- num_clusters
        if (!is.na(c_x_nbr)) m <- m + 1
        if (!is.na(c_y_nbr)) m <- m + 1

        best <- (m - 2.0) * mat_env$m[c_x, c_y] - nodes$Rx[c_x + 1L] - nodes$Rx[c_y + 1L]
        if (!is.na(c_x_nbr)) {
            q_pq <- (m - 2.0) * mat_env$m[c_x_nbr, c_y] - nodes$Rx[c_x_nbr + 1L] - nodes$Rx[c_y + 1L]
            if (q_pq < best) {
                x <- c_x_nbr
                y <- c_y
                best <- q_pq
            }
        }

        if (!is.na(c_y_nbr)) {
            q_pq <- (m - 2.0) * mat_env$m[c_x, c_y_nbr] - nodes$Rx[c_x + 1L] - nodes$Rx[c_y_nbr + 1L]
            if (q_pq < best) {
                x <- c_x
                y <- c_y_nbr
                best <- q_pq
            }
        }

        if (!is.na(c_x_nbr) && !is.na(c_y_nbr)) {
            q_pq <- (m - 2.0) * mat_env$m[c_x_nbr, c_y_nbr] - nodes$Rx[c_x_nbr + 1L] - nodes$Rx[c_y_nbr + 1L]
            if (q_pq < best) {
                x <- c_x_nbr
                y <- c_y_nbr
            }
        }

        x_nbr <- nodes$nbr[x + 1L]
        y_nbr <- nodes$nbr[y + 1L]

        if (is.na(x_nbr) && is.na(y_nbr)) {
            .nnet_join2way(nodes, x, y)
            num_clusters <- num_clusters - 1L
        } else if (is.na(x_nbr)) {
            .nnet_join3way(nodes, mat_env, joins_env, x, y, y_nbr, num_nodes)
            num_nodes <- num_nodes + 2L
            num_active <- num_active - 1L
            num_clusters <- num_clusters - 1L
        } else if (is.na(y_nbr) || (num_active == 4)) {
            .nnet_join3way(nodes, mat_env, joins_env, y, x, x_nbr, num_nodes)
            num_nodes <- num_nodes + 2L
            num_active <- num_active - 1L
            num_clusters <- num_clusters - 1L
        } else {
            num_nodes <- .nnet_join4way(nodes, mat_env, joins_env, x_nbr, x, y, y_nbr, num_nodes)
            num_active <- num_active - 2L
            num_clusters <- num_clusters - 1L
        }
    }

    invisible(NULL)
}

.nnet_expand_nodes <- function(nodes, joins_env, n) {
    x <- nodes$nxt[0L + 1L]
    y <- nodes$nxt[x + 1L]
    z <- nodes$nxt[y + 1L]
    nodes$nxt[z + 1L] <- x
    nodes$prv[x + 1L] <- z

    while (joins_env$top > 0L) {
        u <- joins_env$stack[joins_env$top]
        joins_env$top <- joins_env$top - 1L

        v <- nodes$nbr[u + 1L]
        x <- nodes$ch1[u + 1L]
        y <- nodes$ch2[u + 1L]
        z <- nodes$ch2[v + 1L]

        if (v != nodes$nxt[u + 1L]) {
            tmp <- u
            u <- v
            v <- tmp
            tmp <- x
            x <- z
            z <- tmp
        }

        nodes$prv[x + 1L] <- nodes$prv[u + 1L]
        nodes$nxt[nodes$prv[x + 1L] + 1L] <- x
        nodes$nxt[x + 1L] <- y
        nodes$prv[y + 1L] <- x
        nodes$nxt[y + 1L] <- z
        nodes$prv[z + 1L] <- y
        nodes$nxt[z + 1L] <- nodes$nxt[v + 1L]
        nodes$prv[nodes$nxt[z + 1L] + 1L] <- z
    }

    while (x != 1L) {
        x <- nodes$nxt[x + 1L]
    }

    cycle <- integer(n + 1L)
    a <- x
    i <- 1L
    repeat {
        i <- i + 1L
        cycle[i] <- a
        a <- nodes$nxt[a + 1L]
        if (a == x) break
    }

    cycle
}

# cycle is a 0-prefixed integer vector: cycle[1] == 0, cycle[2:length(cycle)] are node ids.
.nnet_normalize_cycle <- function(cycle) {
    L <- length(cycle)
    pos_of_1 <- 1
    for (i in 1:(L - 1)) {
        if (cycle[i + 1] == 1) {
            pos_of_1 <- i
            break
        }
    }

    last <- L - 1
    pos_prev <- if (pos_of_1 == 1) last else pos_of_1 - 1
    pos_next <- if (pos_of_1 == last) 1 else pos_of_1 + 1

    if (cycle[pos_prev + 1] > cycle[pos_next + 1]) {
        if (pos_of_1 == 1) {
            return(cycle)
        }
        result <- integer(L)
        i <- pos_of_1
        j <- 1L
        while (j < L) {
            j <- j + 1L
            result[j] <- cycle[i + 1]
            i <- if (i < last) i + 1 else 1
        }
        return(result)
    }

    result <- integer(L)
    i <- pos_of_1
    j <- 1L
    while (j < L) {
        j <- j + 1L
        result[j] <- cycle[i + 1]
        i <- if (i > 1) i - 1 else last
    }
    result
}

# Computes a NeighborNet circular ordering ("cycle") of `labels` from the n x n
# distance matrix `mat`. Returns an integer vector of length n + 1: element 1 is
# always 0 (dropped by R's 1-based indexing, e.g. labels[cycle]), and the
# remaining n elements are a permutation of 1:n giving the circular order.
neighbor_net_cycle <- function(labels, mat) {
    n <- length(labels)

    if (n <= 3) {
        return(0:n)
    }

    max_number_of_nodes <- max(3, 3 * n - 5)

    nodes <- .nnet_setup_nodes(n, max_number_of_nodes)

    mat_env <- new.env(parent = emptyenv())
    mat_env$m <- .nnet_setup_matrix(n, mat)

    joins_env <- new.env(parent = emptyenv())
    joins_env$stack <- integer(max_number_of_nodes)
    joins_env$top <- 0L

    .nnet_join_nodes(nodes, mat_env, joins_env, n)

    cycle <- .nnet_expand_nodes(nodes, joins_env, n)

    .nnet_normalize_cycle(cycle)
}
