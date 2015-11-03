#PBIB test modified
##How to use it?
##Step1 : Run order.group2
##STEP2: Run PBIB.test.modif (a modified agricolae PBIB function)

library(agricolae)
library(nlme) 

order.group2 <- 
  function (trt, means, N, MSerror, Tprob, std.err, parameter = 1, 
            snk = 0, DFerror = NULL, alpha = NULL, sdtdif = NULL, vartau = NULL, 
            console) 
  {
    replications <- N
    n <- length(means)
    z <- data.frame(trt, means, N = N, std.err, replications)
    letras <- c(letters[1:26], LETTERS[1:26], 1:120)
    w <- z[order(z[, 2], decreasing = TRUE), ]
    M <- rep("", n)
    k <- 1
    j <- 1
    k <- 1
    cambio <- n
    cambio1 <- 0
    chequeo = 0
    M[1] <- letras[k]
    N <- w$N
    while (j < n) {
      chequeo <- chequeo + 1
      if (chequeo > n) 
        break
      for (i in j:n) {
        nx <- abs(i - j) + 1
        if (nx == 1) 
          nx = 2
        if (snk == 1) 
          Tprob <- qtukey(p = 1 - alpha, nmeans = nx, df = DFerror)
        if (snk == 2) 
          Tprob <- qtukey(p = (1 - alpha)^(nx - 1), nmeans = nx, 
                          df = DFerror)
        if (!is.null(vartau)) {
          ii0 <- w[i, 1]
          jj0 <- w[j, 1]
          if (snk == 3 | snk == 4) 
            sdtdif <- sqrt(vartau[ii0, ii0] + vartau[jj0, 
                                                     jj0] - 2 * vartau[ii0, jj0])
          if (snk == 5 | snk == 6) 
            sdtdif <- sqrt(vartau[ii0, jj0])
          if (snk == 3 | snk == 5) 
            Tprob <- qt(p = (1 - alpha/2), df = DFerror)
          if (snk == 4 | snk == 6) 
            Tprob <- qtukey(p = (1 - alpha), nmeans = n, 
                            df = DFerror)
        }
        if (is.null(sdtdif)) 
          minimo <- Tprob * sqrt(parameter * MSerror * 
                                   (1/N[i] + 1/N[j]))
        if (!is.null(sdtdif)) 
          minimo <- Tprob * sdtdif
        s <- abs(w[i, 2] - w[j, 2]) <= minimo
        if (s) {
          if (lastC(M[i]) != letras[k]) 
            M[i] <- paste(M[i], letras[k], sep = "")
        }
        else {
          k <- k + 1
          cambio <- i
          cambio1 <- 0
          ja <- j
          for (jj in cambio:n) M[jj] <- paste(M[jj], sep = "")
          M[cambio] <- paste(M[cambio], letras[k], sep = "")
          for (v in ja:cambio) {
            nx <- abs(v - cambio) + 1
            if (nx == 1) 
              nx = 2
            if (snk == 1) 
              Tprob <- qtukey(p = 1 - alpha, nmeans = nx, 
                              df = DFerror)
            if (snk == 2) 
              Tprob <- qtukey(p = (1 - alpha)^(nx - 1), 
                              nmeans = nx, df = DFerror)
            if (!is.null(vartau)) {
              ii0 <- w[i, 1]
              jj0 <- w[j, 1]
              if (snk == 3 | snk == 4) 
                sdtdif <- sqrt(vartau[ii0, ii0] + vartau[jj0, 
                                                         jj0] - 2 * vartau[ii0, jj0])
              if (snk == 5 | snk == 6) 
                sdtdif <- sqrt(vartau[ii0, jj0])
              if (snk == 3 | snk == 5) 
                Tprob <- qt(p = (1 - alpha/2), df = DFerror)
              if (snk == 4 | snk == 6) 
                Tprob <- qtukey(p = (1 - alpha), nmeans = n, 
                                df = DFerror)
            }
            if (is.null(sdtdif)) 
              minimo <- Tprob * sqrt(parameter * MSerror * 
                                       (1/N[i] + 1/N[j]))
            if (!is.null(sdtdif)) 
              minimo <- Tprob * sdtdif
            if (abs(w[v, 2] - w[cambio, 2]) > minimo) {
              j <- j + 1
              cambio1 <- 1
            }
            else break
          }
          break
        }
      }
      if (cambio1 == 0) 
        j <- j + 1
    }
    w <- data.frame(w, stat = M)
    trt <- as.character(w$trt)
    means <- as.numeric(w$means)
    means1 <- signif(means, 4)
    N <- as.numeric(w$N)
    std.err <- as.numeric(w$std.err)
    replications <- w$replications
    cmax <- max(nchar(trt))
    trt <- paste(trt, "                                      ")
    trt <- substr(trt, 1, cmax)
    output1 <- cbind(trt, means)
    rownames(output1) <- M
    for (i in 1:n) {
      if (console) 
        cat(rownames(output1)[i], "\t", output1[i, 1], "\t", 
            means1[i], "\n")
    }
    output <- data.frame(trt, means, M, N = replications, std.err)
    invisible(output)
  }


###This a modified function of PBIB.test (agricolae package)
PBIB.test.modif <- function (block, trt, replication, y, k, method = c("REML", "ML", 
                                                                       "VC"), test = c("lsd", "tukey"), alpha = 0.05, console = FALSE, 
                             group = TRUE) 
{
  test <- match.arg(test)
  if (test == "lsd") 
    snk = 3
  if (test == "tukey") 
    snk = 4
  method <- match.arg(method)
  if (method == "REML") 
    nMethod <- "Residual (restricted) maximum likelihood"
  if (method == "ML") 
    nMethod <- "Maximum likelihood"
  if (method == "VC") 
    nMethod <- "Variances component model"
  name.y <- paste(deparse(substitute(y)))
  name.r <- paste(deparse(substitute(replication)))
  name.b <- paste(deparse(substitute(block)))
  name.trt <- paste(deparse(substitute(trt)))
  block.adj <- as.factor(block)
  trt.adj <- as.factor(trt)
  replication <- as.factor(replication)
  mean.trt <- as.matrix(by(y, trt, function(x) mean(x, na.rm = TRUE)))
  mi <- as.matrix(by(y, trt, function(x) min(x, na.rm = TRUE)))
  ma <- as.matrix(by(y, trt, function(x) max(x, na.rm = TRUE)))
  n.rep <- as.matrix(by(y, trt, function(x) length(na.omit(x))))
  sds <- as.matrix(by(y, trt, function(x) sd(x, na.rm = TRUE)))
  std <- sds
  indice <- rownames(mean.trt)
  ntr <- nlevels(trt.adj)
  r <- nlevels(replication)
  s <- ntr/k
  obs <- length(na.omit(y))
  if (method == "VC" & obs != r * ntr) {
    if (console) 
      cat("\nWarning.. incomplete repetition. Please you use method REML or ML\n")
    return()
  }
  modelo <- formula(paste(name.y, "~ replication + trt.adj+ block.adj%in%replication"))
  model <- lm(modelo)
  glerror <- df.residual(model)
  ANOVA <- anova(model)
  rownames(ANOVA)[2] <- name.trt
  CMerror <- as.numeric(deviance(model)/glerror)
  Mean <- mean(y, na.rm = TRUE)
  if (method == "VC") {
    rownames(ANOVA) <- c(name.r, paste(name.trt, ".unadj", 
                                       sep = ""), paste(name.b, "/", name.r, sep = ""), 
                         "Residual")
  }
  if (method == "REML" | method == "ML") {
    trt.adj <- as.factor(trt)
    if (method == "REML") {
      modlmer <- lme(y ~ 0 + trt.adj, random = ~1 | replication/block.adj, 
                     method = "REML")
      model <- lme(y ~ trt.adj, random = ~1 | replication/block.adj, 
                   method = "REML")
    }
    if (method == "ML") {
      modlmer <- lme(y ~ 0 + trt.adj, random = ~1 | replication/block.adj, 
                     method = "ML")
      model <- lme(y ~ trt.adj, random = ~1 | replication/block.adj, 
                   method = "ML")
    }
    Afm <- anova(model)
    VarRand <- matrix(rep(0, 3), nrow = 3)
    VarRand[c(1, 3), 1] <- as.numeric(VarCorr(model)[4:5, 
                                                     1])
    VarRand[2, 1] <- as.numeric(VarCorr(model)[2, 1])
    CMerror <- as.numeric(VarRand[3, 1])
    VarRand <- data.frame(VarRand)
    names(VarRand) <- "Variance"
    rownames(VarRand) <- c(paste(name.b, ":", name.r, sep = ""), 
                           name.r, "Residual")
    ANOVA <- ANOVA[c(2, 4), ]
    ANOVA[2, 3] <- CMerror
    ANOVA[1, 1] <- Afm[2, 1]
    ANOVA[1, 4] <- Afm[2, 3]
    ANOVA[1, 5] <- Afm[2, 4]
    ANOVA[1, 3] <- ANOVA[2, 3] * ANOVA[1, 4]
    ANOVA[, 2] <- ANOVA[, 1] * ANOVA[, 3]
    ANOVA[2, 4:5] <- NA
  }
  b <- s * r
  glt <- ntr - 1
  if (method == "VC") {
    SCt <- anova(model)[2, 2]
    Ee <- deviance(model)/glerror
    Eb <- anova(model)[3, 3]
    X <- rep(0, obs * ntr)
    dim(X) <- c(obs, ntr)
    for (i in 1:obs) {
      tr <- trt[i]
      X[i, tr] <- 1
    }
    R <- rep(0, obs * r)
    dim(R) <- c(obs, r)
    for (i in 1:obs) {
      rp <- replication[i]
      R[i, rp] <- 1
    }
    Z <- rep(0, obs * b)
    dim(Z) <- c(obs, b)
    for (i in 1:obs) {
      rb <- block[i]
      Z[i, rb] <- 1
    }
    N <- t(X) %*% Z
    In <- diag(1, obs)
    c0 <- t(Z) %*% (In - (1/r) * X %*% t(X)) %*% y
    Js <- diag(s)
    Ir <- diag(r)
    Jr <- matrix(1, r, r)
    Js <- matrix(1, s, s)
    Ib <- diag(b)
    Iv <- diag(ntr)
    q <- k - floor(k/s) * s
    if (q <= s/2) 
      g <- floor(k/s)
    if (q > s/2) 
      g <- floor(k/s) + 1
    phi <- r * (Eb - Ee)/((r - 1) * Ee)
    lambda <- 1/(r * k * (1/phi + 1) - k)
    W <- t(N) %*% N - k * Ib - g * kronecker((Jr - Ir), Js)
    inversa <- ginv(Ib - lambda * W)
    tauIntra <- t(X) %*% y/r - lambda * N %*% inversa %*% 
      c0
    vartau <- (Ee/r) * (Iv + lambda * N %*% inversa %*% t(N))
    dvar <- sqrt(diag(vartau))
  }
  if (method == "REML" | method == "ML") {
    tauIntra <- fixef(modlmer)
    vartau <- vcov(modlmer)
    DIA <- as.matrix(vartau)
    dvar <- sqrt(diag(DIA))
  }
  ntr0 <- ncol(vartau)
  if (ntr0 < ntr) {
    ntr <- ntr0
    n.rep <- na.omit(n.rep)
    std <- na.omit(std)
    mean.trt <- na.omit(mean.trt)
    mi <- na.omit(mi)
    ma <- na.omit(ma)
  }
  vardif <- matrix(0, ntr, ntr)
  for (i in 1:(ntr - 1)) {
    for (j in (i + 1):ntr) {
      vardif[i, j] <- vartau[i, i] + vartau[j, j] - 2 * 
        vartau[i, j]
      vardif[j, i] <- vardif[i, j]
    }
  }
  media <- mean(y, na.rm = TRUE)
  if (console) {
    cat("\nANALYSIS PBIB: ", name.y, "\n\nClass level information\n")
    cat(name.b, ":", b, "\n")
    cat(name.trt, ":", ntr)
    cat("\n\nNumber of observations: ", length(y), "\n\n")
    cat("Estimation Method: ", nMethod, "\n\n")
  }
  Fstat <- data.frame(c(AIC(model), BIC(model)))
  rownames(Fstat) <- c("AIC", "BIC")
  names(Fstat) <- "Fit Statistics"
  if (method == "REML" | method == "ML") {
    Fstat <- rbind(Fstat, model$logLik)
    rownames(Fstat)[3] <- "-2 Res Log Likelihood"
    if (console) {
      cat("Parameter Estimates\n")
      print(VarRand)
      cat("\n")
    }
  }
  if (console) {
    print(Fstat)
    cat("\n")
  }
  CV <- sqrt(CMerror) * 100/media
  design <- data.frame(. = c(ntr, k, b/r, r))
  rownames(design) <- c(name.trt, paste(name.b, "size"), paste(name.b, 
                                                               "/", name.r, sep = ""), name.r)
  E <- (ntr - 1) * (r - 1)/((ntr - 1) * (r - 1) + r * (s - 
                                                         1))
  if (console) {
    print(ANOVA)
    cat("\ncoefficient of variation:", round(CV, 1), "%\n")
    cat(name.y, "Means:", media, "\n")
    cat("\nParameters PBIB\n")
    print(design)
    cat("\nEfficiency factor", E, "\n")
    cat("\nComparison test", test, "\n")
  }
  parameters <- data.frame(treatments = ntr, blockSize = k, 
                           blocks = s, r = r)
  statistics <- data.frame(Efficiency = E, Mean = Mean, CV = CV)
  rownames(parameters) <- " "
  rownames(statistics) <- " "
  comb <- combn(ntr, 2)
  nn <- ncol(comb)
  dif <- rep(0, nn)
  stdt <- rep(0, nn)
  pvalue <- rep(0, nn)
  for (k in 1:nn) {
    i <- comb[1, k]
    j <- comb[2, k]
    dif[k] <- tauIntra[i] - tauIntra[j]
    stdt[k] <- sqrt(vartau[i, i] + vartau[j, j] - 2 * vartau[i, 
                                                             j])
    tc <- abs(dif[k])/stdt[k]
    if (test == "lsd") 
      pvalue[k] <- 2 * round(1 - pt(tc, glerror), 6)
    if (test == "tukey") 
      pvalue[k] <- round(1 - ptukey(tc, ntr, glerror), 
                         6)
  }
  tr.i <- comb[1, ]
  tr.j <- comb[2, ]
  groups <- NULL
  if (group) {
    groups <- order.group2(trt = 1:ntr, tauIntra, n.rep, MSerror = NULL, 
                           Tprob = NULL, std.err = dvar, parameter = 1, snk, 
                           DFerror = glerror, alpha, sdtdif = 1, vartau, console = FALSE)
    names(groups)[2] <- "mean.adj"
    rownames(groups) <- groups$trt
    indices <- as.numeric(as.character(groups$trt))
    groups$trt <- indice[indices]
    names(groups)[1] <- name.trt
    groups <- groups[, 1:3]
  }
  cat("\n<<< to see the objects: means, comparison and groups. >>>\n\n")
  comparison <- data.frame(Difference = dif, stderr = stdt, 
                           pvalue = pvalue)
  means <- data.frame(means = mean.trt, trt = 1:ntr, mean.adj = as.numeric(tauIntra), 
                      SE = dvar, r = n.rep, std, Min = mi, Max = ma)
  names(means)[1] <- name.y
  rownames(comparison) <- paste(indice[tr.i], indice[tr.j], 
                                sep = " - ")
  output <- list(ANOVA = ANOVA, method = nMethod, parameters = parameters, 
                 statistics = statistics, Fstat = Fstat, comparison = comparison, 
                 means = means, groups = groups, vartau = vartau)
  invisible(output)
}


output <- PBIB.test.modif(block = dbook$block,trt = dbook$Genotype,replication = dbook$replication,y = dbook$yield,k = 3,method = "REML",test = "lsd",alpha = 0.05,console = TRUE,group = TRUE)



