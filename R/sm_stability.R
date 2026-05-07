sm_path_subsets <- function(sinks, n.nodes) {
    bitmask <- 2^(seq_len(n.nodes) - 1)
    mask <- 2^n.nodes - 1
    subsets <- integer(n.nodes)

    for (i in seq_len(n.nodes)) {
        subsets[i] <- mask
        sink <- sinks[mask + 1]
        mask <- bitXor(mask, bitmask[sink])
    }

    subsets
}

sm_stability_bdeu <- function(
    data,
    node.sizes,
    ess,
    aflml,
    bps,
    sinks,
    F,
    ifm,
    bound.type = c("cons", "tight"),
    scope = c("path", "all")) {
    bound.type <- match.arg(bound.type)
    scope <- match.arg(scope)

    if (scope != "path") {
        stop("stability.scope = 'all' not yet implemented.")
    }

    bounds <- .Call(
        "bnstruct_bdeu_loo_bounds",
        data, node.sizes, ifm, ess,
        PACKAGE = "bnstruct"
    )

    lower_delta <- bounds$lower
    upper_delta <- bounds$upper

    if (bound.type == "cons") {
        Bg <- .Call(
            "bnstruct_compute_bg_cons",
            lower_delta, upper_delta, ifm,
            PACKAGE = "bnstruct"
        )
    } else {
        Bg <- .Call(
            "bnstruct_compute_bg_tight",
            aflml, bps, lower_delta, upper_delta, ifm,
            PACKAGE = "bnstruct"
        )
    }

    BF <- .Call(
        "bnstruct_compute_bf",
        Bg,
        PACKAGE = "bnstruct"
    )

    blocal <- pmax(abs(lower_delta), abs(upper_delta))

    n.nodes <- ncol(data)
    bitmask <- 2^(seq_len(n.nodes) - 1)

    subsets <- sm_path_subsets(sinks, n.nodes)

    sink_failures <- list()
    parent_failures <- list()

    for (mask in subsets) {
        nodes <- which(bitAnd(mask, bitmask) > 0)
        if (length(nodes) == 0) {
            next
        }

        i_star <- sinks[mask + 1]

        if (length(nodes) >= 2) {
            A_vals <- numeric(length(nodes))
            BA_vals <- numeric(length(nodes))

            for (idx in seq_along(nodes)) {
                i <- nodes[idx]
                U_mask <- bitXor(mask, bitmask[i])
                u_idx <- U_mask + 1
                A_vals[idx] <- F[u_idx] + aflml[i, bps[i, u_idx]]
                BA_vals[idx] <- BF[u_idx] + Bg[i, u_idx]
            }

            idx_star <- which(nodes == i_star)
            if (length(idx_star) == 1) {
                comp_idx <- setdiff(seq_along(nodes), idx_star)
                if (length(comp_idx) > 0) {
                    if (bound.type == "tight") {
                        risk_vals <- A_vals[comp_idx] + BA_vals[comp_idx]
                        risk_idx <- comp_idx[which.max(risk_vals)]

                        A_star <- A_vals[idx_star]
                        BA_star <- BA_vals[idx_star]
                        lhs <- A_star - BA_star
                        rhs <- A_vals[risk_idx] + BA_vals[risk_idx]

                        if (!is.finite(lhs) || !is.finite(rhs) || lhs <= rhs) {
                            sink_failures[[length(sink_failures) + 1]] <- list(
                                subset = mask,
                                chosen = i_star,
                                competitor = nodes[risk_idx],
                                A_chosen = A_star,
                                A_competitor = A_vals[risk_idx],
                                BA_chosen = BA_star,
                                BA_competitor = BA_vals[risk_idx],
                                lhs = lhs,
                                rhs = rhs,
                                gap = A_star - A_vals[risk_idx],
                                bound = BA_star + BA_vals[risk_idx],
                                type = "sink"
                            )
                        }
                    } else {
                        comp_scores <- A_vals[comp_idx]
                        if (any(is.finite(comp_scores))) {
                            best_comp_idx <- comp_idx[which.max(comp_scores)]
                            A_star <- A_vals[idx_star]
                            A_second <- A_vals[best_comp_idx]
                            BA_star <- BA_vals[idx_star]
                            BA_second <- BA_vals[best_comp_idx]
                            gap <- A_star - A_second
                            bound <- BA_star + BA_second

                            if (!is.finite(gap) || !is.finite(bound) || gap <= bound) {
                                sink_failures[[length(sink_failures) + 1]] <- list(
                                    subset = mask,
                                    chosen = i_star,
                                    runner_up = nodes[best_comp_idx],
                                    gap = gap,
                                    bound = bound,
                                    type = "sink"
                                )
                            }
                        }
                    }
                }
            }
        }

        U_mask <- bitXor(mask, bitmask[i_star])
        u_idx <- U_mask + 1
        p_star_idx <- bps[i_star, u_idx]
        best_score <- aflml[i_star, p_star_idx]

        if (is.finite(best_score)) {
            if (bound.type == "tight") {
                max_risk <- -Inf
                p_risk <- NA_integer_

                p_mask <- U_mask
                repeat {
                    p_idx <- p_mask + 1
                    if (ifm[i_star, p_idx] && p_idx != p_star_idx) {
                        risk <- aflml[i_star, p_idx] + upper_delta[i_star, p_idx]
                        if (risk > max_risk) {
                            max_risk <- risk
                            p_risk <- p_idx
                        }
                    }

                    if (p_mask == 0) {
                        break
                    }
                    p_mask <- bitAnd(p_mask - 1, U_mask)
                }

                if (is.finite(max_risk)) {
                    lhs <- best_score + lower_delta[i_star, p_star_idx]
                    rhs <- max_risk
                    parent_gap <- best_score - aflml[i_star, p_risk]
                    parent_bound <- upper_delta[i_star, p_risk] - lower_delta[i_star, p_star_idx]

                    if (!is.finite(lhs) || !is.finite(rhs) || lhs <= rhs) {
                        parent_failures[[length(parent_failures) + 1]] <- list(
                            subset = mask,
                            node = i_star,
                            best_parent_set = p_star_idx,
                            competitor_parent_set = p_risk,
                            lhs = lhs,
                            rhs = rhs,
                            parent_gap = parent_gap,
                            bound = parent_bound,
                            type = "parent"
                        )
                    }
                }
            } else {
                max_other_score <- -Inf
                max_other_blocal <- -Inf
                blocal_star <- blocal[i_star, p_star_idx]

                p_mask <- U_mask
                repeat {
                    p_idx <- p_mask + 1
                    if (ifm[i_star, p_idx]) {
                        score <- aflml[i_star, p_idx]
                        bl <- blocal[i_star, p_idx]
                        if (p_idx != p_star_idx) {
                            if (score > max_other_score) {
                                max_other_score <- score
                            }
                            if (bl > max_other_blocal) {
                                max_other_blocal <- bl
                            }
                        }
                    }

                    if (p_mask == 0) {
                        break
                    }
                    p_mask <- bitAnd(p_mask - 1, U_mask)
                }

                if (is.finite(max_other_score)) {
                    parent_gap <- best_score - max_other_score
                    parent_bound <- blocal_star + max_other_blocal

                    if (!is.finite(parent_gap) || !is.finite(parent_bound) || parent_gap <= parent_bound) {
                        parent_failures[[length(parent_failures) + 1]] <- list(
                            subset = mask,
                            node = i_star,
                            best_parent_set = p_star_idx,
                            parent_gap = parent_gap,
                            bound = parent_bound,
                            type = "parent"
                        )
                    }
                }
            }
        }
    }

    stable_order <- length(sink_failures) == 0
    stable_dag <- stable_order && length(parent_failures) == 0

    list(
        stable_order = stable_order,
        stable_dag = stable_dag,
        sink_failures = sink_failures,
        parent_failures = parent_failures,
        witnesses = c(sink_failures, parent_failures),
        bound_type = bound.type,
        scope = scope
    )
}
