ARI <- function(true_cl, pred_cl) {
  res <- aricode::sortPairs(true_cl, pred_cl)

  n <- length(true_cl)
  K <- length(res$levels$c1)

  stot <- sum(choose(res$nij, 2), na.rm = TRUE)
  srow <- sum(choose(res$ni., 2), na.rm = TRUE)
  scol <- sum(choose(res$n.j, 2), na.rm = TRUE)

  expectedIndex <- (srow * scol) / (choose(n, 2))
  maximumIndex <- (srow + scol) / 2

  a <- (stot - expectedIndex + 0.5*n) / (maximumIndex - expectedIndex + n)

  return(list(ARI=a, stot=stot, srow=srow, scol=scol))
}

load_ret <- function() {
  exp_cases <- read_csv("../experimental_cases.csv")
  K.max <- 100
  tt.max <- 10
  methods <- c("xmeans", "gmeans", "dipmeans", "pgmeans", "smlsom", "mml_em")

  df <- data.frame()
  for (i in 1:nrow(exp_cases)) {
    cat(i, date(), "\n")

    ind <- exp_cases[i,]$index
    n <- exp_cases[i,]$n

    Dataname <- sprintf("Gaussians%03d_n=%d.RData", ind, n)

    load(str_c("ground_truth/uniform/", Dataname))

    for (method in methods) {
      if (method == "dipmeans" & n == 27000)
        next

      Valname <- sprintf("Gaussians%03d_n=%d_%s.RData", ind, n, method)

      load(str_c("uniform/values/", Valname))

      df_evals["n"] <- n
      df_evals["cARI"] <- 0
      df_evals["stot"] <- 0
      df_evals["srow"] <- 0
      df_evals["scol"] <- 0

      for (k in 1:K.max) {
        Modelname <- sprintf("Gaussians%03d_%03d_n=%d_%s.RData", ind, k, n, method)

        load(str_c("uniform/models/", Modelname))

        true_cl <- label_array[,k]

        v <- df_evals %>% filter(dataset==k)

        m_pos <- which(!sapply(Models, is.null))

        if (length(m_pos) == 0) next

        vals <- lapply(Models[m_pos], function(m) ARI(true_cl, m$classes))

        df_vals <- bind_rows(vals) %>% rename(cARI=ARI)

        v_pos <- which(v$K!=0)

        pos <- (k-1)*tt.max + v_pos

        df_evals[pos, "cARI"] <- df_vals$cARI
        df_evals[pos, "stot"] <- df_vals$stot
        df_evals[pos, "srow"] <- df_vals$srow
        df_evals[pos, "scol"] <- df_vals$scol
      }
      df <- bind_rows(df, df_evals)
    }
  }
  return(df)
}

contrast <- function(df, type=c("treatment","sum")) {
  type <- match.arg(type)

  if (type == "treatment") {
    contr <- matrix(c(1, 0, 0, 0, 0,
                      0, 1, 0, 0, 0,
                      0, 0, 1, 0, 0,
                      0, 0, 0, 1, 0), nrow=5)
    dimnames(contr) <- list(methods, methods[-length(methods)])
    contrasts(df$method) <- contr

    contr <- matrix(c(0, 1, 0, 0,
                      0, 0, 1, 0,
                      0, 0, 0, 1), nrow=4)
    dimnames(contr) <- list(1:4, 2:4)
    contrasts(df$cov_type) <- contr


    contr <- matrix(c(1,0,0,0,0,1), nrow=3)
  } else if (type == "sum") {
    contr <- matrix(c(1, 0, 0, 0, 0, -1,
                      0, 1, 0, 0, 0, -1,
                      0, 0, 1, 0, 0, -1,
                      0, 0, 0, 1, 0, -1,
                      0, 0, 0, 0, 1, -1
    ), nrow=6)
    dimnames(contr) <- list(methods, methods[-length(methods)])
    contrasts(df$method) <- contr

    contr <- matrix(c(-1, 1, 0, 0,
                      -1, 0, 1, 0,
                      -1, 0, 0, 1), nrow=4)
    dimnames(contr) <- list(1:4, 2:4)
    contrasts(df$cov_type) <- contr

    ##contr <- matrix(c(1, -1, 0,
    ##                  0, -1, 1), nrow=3)
    contr <- matrix(c(-1, 1, 0,
                      -1, 0, 1), nrow=3)
  }

  ### dimensionality
  #dimnames(contr) <- list(c(2, 6, 18), c(2, 18))
  dimnames(contr) <- list(c(2, 6, 18), c(6, 18))
  contrasts(df$p) <- contr

  ### cluster number
  #dimnames(contr) <- list(c(3, 6, 12), c(3, 12))
  dimnames(contr) <- list(c(3, 6, 12), c(6, 12))
  contrasts(df$K) <- contr

  ### overlap
  #dimnames(contr) <- list(c(0.01, 0.05, 0.1), c(0.01, 0.1))
  dimnames(contr) <- list(c(0.01, 0.05, 0.1), c(0.05, 0.1))
  contrasts(df$omega) <- contr

  ### sample number
  #dimnames(contr) <- list(c(3000, 9000, 27000), c(3000, 27000))
  dimnames(contr) <- list(c(3000, 9000, 27000), c(9000, 27000))
  contrasts(df$n) <- contr

  return(df)
}

get_rcoef <- function(object)
{
  xl <- object$xlevels
  if (!length(xl))
    return(as.list(coef(object)))
  Terms <- terms(object)
  tl <- attr(Terms, "term.labels")
  int <- attr(Terms, "intercept")
  facs <- attr(Terms, "factors")[-1, , drop = FALSE]
  Terms <- delete.response(Terms)
  mf <- object$model %||% model.frame(object)
  vars <- dimnames(facs)[[1]]
  xtlv <- lapply(mf[, vars, drop = FALSE], levels)
  nxl <- pmax(lengths(xtlv), 1L)
  lterms <- apply(facs, 2L, function(x) prod(nxl[x > 0]))
  nl <- sum(lterms)
  args <- sapply(vars, function(i) if (nxl[i] == 1)
    rep.int(1, nl)
    else factor(rep.int(xtlv[[i]][1L], nl), levels = xtlv[[i]]),
    simplify = FALSE)
  dummy <- do.call(data.frame, args)
  names(dummy) <- vars
  pos <- 0L
  rn <- rep.int(tl, lterms)
  rnn <- character(nl)
  for (j in tl) {
    i <- vars[facs[, j] > 0]
    ifac <- i[nxl[i] > 1]
    lt.j <- lterms[[j]]
    if (length(ifac) == 0L) {
      rnn[pos + 1L] <- j
    }
    else {
      p.j <- pos + seq_len(lt.j)
      if (length(ifac) == 1L) {
        dummy[p.j, ifac] <- x.i <- xtlv[[ifac]]
        rnn[p.j] <- as.character(x.i)
      }
      else {
        tmp <- expand.grid(xtlv[ifac], KEEP.OUT.ATTRS = FALSE)
        dummy[p.j, ifac] <- tmp
        rnn[p.j] <- apply(as.matrix(tmp), 1L, paste,
                          collapse = ":")
      }
    }
    pos <- pos + lt.j
  }
  attr(dummy, "terms") <- attr(mf, "terms")
  lcontr <- object$contrasts
  lci <- vapply(dummy, is.factor, NA)
  lcontr <- lcontr[names(lci)[lci]]
  mm <- model.matrix(Terms, dummy, lcontr, xl)
  if (anyNA(mm)) {
    warning("some terms will have NAs due to the limits of the method")
    mm[is.na(mm)] <- NA
  }
  coef <- object$coefficients
  coef[is.na(coef)] <- 0
  asgn <- attr(mm, "assign")

  V <- vcov(object)
  V[is.na(V)] <- 0

  coef_mat <- data.frame(matrix(0, nrow(mm), 4))
  colnames(coef_mat) <- c("factor", "level", "est", "se")
  coef_mat[,1] <- rn
  coef_mat[,2] <- rnn

  for (j in seq_along(tl)) {
    keep <- which(asgn == j)
    cf <- coef[keep]

    ij <- rn == tl[j]

    v <- V[keep, keep, drop=F]
    Cmat <- mm[ij, keep, drop=F]

    coef_mat[ij, "est"] <- setNames(drop(Cmat %*% cf), rnn[ij])
    coef_mat[ij, "se"] <- setNames(sqrt(diag(Cmat %*% v %*% t(Cmat))), rnn[ij])
  }

  coef_mat <- rbind(data.frame(factor="Intercept", level="Intercept",
                               est=as.numeric(coef[int]), se=sqrt(V[int, int])),
                    coef_mat)
  return(coef_mat)
}

plot.interval_main <- function(df_coef, xlim=c(-0.8, 0.8),
                               cex.axis=1, line=3.5) {

  tl <- unique(df_coef$factor)
  n_labs <- table(df_coef$factor)[tl]
  len_labs <- nrow(df_coef)

  levs <- df_coef$level

  plot(0,0, xlim=xlim, ylim=c(-0.5, len_labs-0.5),
       xlab="Coefficient", ylab="", yaxt="n", type="n")
  abline(v=0, lty=2)
  abline(h=len_labs - cumsum(n_labs) - 0.5, lty=2, col="gray", lwd=0.5)

  axis(side=2,
       at=(len_labs-1):0,
       labels=levs,
       las=2,
       cex.axis=0.5*cex.axis)

  axis(side=2,
       at=len_labs - cumsum(n_labs) - 0.5,
       labels=names(n_labs),
       las=2,
       line=line,
       cex.axis=0.6*cex.axis)

  for (i in 1:len_labs) {
    est <- df_coef[i, "est"]
    se <- df_coef[i, "se"]

    int <- c(est-2*se, est+2*se)

    vpos <- len_labs - i

    segments(x0=est-2*se, y0=vpos, x1=est+2*se, y1=vpos, lwd=0.5)
    col <- ifelse(int[1] < 0 & int[2] > 0, 1, 2)
    points(est, vpos, pch=20, cex=0.5, col=col)
  }
}

plot.interval_inter <- function(df_coef, xlim=c(-0.8, 0.8)) {

  intervals <- matrix(0, nrow(df_coef), 2)
  intervals[,1] <- df_coef$est - 2*df_coef$se
  intervals[,2] <- df_coef$est + 2*df_coef$se

  df_coef <- df_coef[!(intervals[,1] <= 0 & intervals[,2] >= 0), ]

  tl <- unique(df_coef$factor)
  n_labs <- table(df_coef$factor)[tl]
  len_labs <- nrow(df_coef)

  plot(0,0, xlim=xlim, ylim=c(0, len_labs),
       xlab="Coefficient", ylab="", yaxt="n", type="n")
  abline(v=0, lty=2)
  abline(h=len_labs - cumsum(n_labs) - 0.5, lty=2, col="gray", lwd=0.5)

  axis(side=2,
       at=len_labs - cumsum(n_labs) - 0.5,
       labels=names(n_labs),
       las=2,
       line=0.5,
       cex.axis=0.5)

  for (i in 1:len_labs) {
    est <- df_coef[i, "est"]
    se <- df_coef[i, "se"]

    int <- c(est-2*se, est+2*se)

    vpos <- len_labs - i

    if (int[1] <= 0 & int[2] >= 0) {
      next
    } else {
      segments(x0=est-2*se, y0=vpos, x1=est+2*se, y1=vpos, lwd=0.5)
      points(est, vpos, pch=20, cex=0.5, col=2)
    }
  }
}

plot_radar.base <- function(labels, trans, cex.text=0.7) {
  symbols(0, 0, circles=1, inches=F, asp=1,
          xlim=c(-1, 1), ylim=c(-1, 1), fg="gray",
          ann=F, bty="n", xaxt="n", yaxt="n")
  symbols(0, 0, circles=0.5, inches=F, asp=1, add=T, lty=2, fg="gray")

  cex.text <- 0.7
  for (i in 1:nrow(trans)) {
    adj <- ifelse(i==1, 0.5, 0.5 * (sign(trans[i,1])+1 ))
    segments(0, 0, trans[i,1], trans[i,2], col="gray", lwd=0.5)
    text(trans[i,1]+0.05*sign(trans[i,1]), trans[i,2]+0.05*sign(trans[i,2]),
         labels=labels[i], adj=adj, cex=cex.text)
  }
}

plot_radar.line <- function(values, trans, col, value_name="ARI_mn", error_name="ARI_sd") {
  coordinates <- values[[value_name]] * trans
  errors <- values[[error_name]] / 10
  symbols(coordinates, circles=errors*2, inches=F, asp=1, add=T, lty=1, fg=col, lwd=0.35)
  points(rbind(coordinates, coordinates[1,]), type="o", cex=0.3, pch=20, col=col)
}

plot_radar.omega <- function(df, labels, trans, cols, methods=methods, omegas=c(0.01, 0.05, 0.1)) {
  par(mfrow=c(1,3), mar=c(0, 0, 1, 3))
  for (f in 1:length(omegas)) {
    omega <- omegas[f]

    df_sub <- df %>% filter(omega == !!omega)

    plot_radar.base(labels, trans)

    mtext(text=bquote(omega==.(omega)), side=3, line=-2)

    if (f == 1)
      legend("topright", legend=methods, col=cols, lty=1)

    if (f == 3)
      legend("topright", legend=str_c("ARI=", c(1, 0.5)), lty=1:2, col="gray")

    for (j in 1:length(methods)) {
      values <- df_sub %>% filter(method==!!methods[j]) %>% arrange(p, n)
      if (methods[j] == "dipmeans")
        plot_radar.line(values, trans[-length(labels),], col=cols[j])
      else
        plot_radar.line(values, trans, col=cols[j])
    }
  }
}
