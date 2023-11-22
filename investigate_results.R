library(tidyverse)
library(grpreg)
library(gtsummary)
source("functions.R")

### =====================================================================
###              data load (1476000 records)
### =====================================================================

methods <- c("xmeans", "gmeans", "dipmeans", "pgmeans", "smlsom", "mml_em")

df <- read_csv("result_six_algorithms.csv")
## Columns
### method: six cluster search algorithms
### p:        Dimensionality (level: 2, 6, 18)
### n:        Sample size (level: 3000, 9000, 27000)
### omega:    Cluster overlap (level: 0.01, 0.05, 0.1)
### K:        The number of clusters (level: 3, 6, 12)
### cov_type: Covariance type (level: 1, 2, 3, 4)
### dataset:  Index of dataset (1, 2, ..., 100)
### rep:      Index of algorithm trial (1, 2, ..., 10)
### estK:     The estimated number of clusters by the algorithm
### elapsed:  Computation time of algorithm run
### valid:    Value of information criterion corresponding to each algorithm
### cARI:     cARI
### NMI:      Normalized Mutual Information
### ARI:      Adjusted Rand Index
### stot:     Value for calculating ARI (see aricode package)
### srow:     Value for calculating ARI (see aricode package)
### scol:     Value for calculating ARI (see aricode package)

## Choosing the best result from 10 runs
df_ret <- df %>%
    group_by(method, p, n, omega, K, cov_type, dataset) %>%
    reframe(
        minInd = ifelse(all(is.na(valid)), NA, which.min(valid)),
        valid=valid[minInd],
        estK=estK[minInd],
        NMI=NMI[minInd],
        ARI=ARI[minInd],
        cARI=cARI[minInd],
        elapsed=elapsed[minInd]
    )

df.sum <- df_ret %>%
    dplyr::select(method, p, n, omega, K, cov_type, dataset, estK, elapsed, cARI, NMI) %>%
    rename(ARI=cARI) %>%
    mutate(method = factor(method, levels=methods),
           diffK = estK - K,
           p = as.factor(p),
           n = as.factor(n),
           omega = as.factor(omega),
           K = as.factor(K),
           cov_type = as.factor(cov_type))

## sum-to-zero contrast
df.sum <- contrast(df.sum, type = "sum")

## computation time
df.sum %>%
    filter(estK!=0) %>%
    ggplot(aes(x=ARI, y=log10(elapsed), color=method)) +
    geom_point(alpha=0.6, shape=21, size=0.8) +
    facet_grid(K~p+n,
               labeller = labeller(K = label_both,
                                   p = label_both,
                                   n = label_both)) +
    guides(color = guide_legend(override.aes = list(size=2, alpha=1, shape=19))) +
    xlab("cARI") + ylab("log10(elapsed time[sec])") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ============================================================
##             Analysis Data (1476 records)
## ============================================================

## compression
df_fac <- df.sum %>%
    group_by(method, p, n, omega, K, cov_type) %>%
    summarize(
        n_valid=sum(!is.na(ARI)) / n(),
        ARI=mean(ARI, na.rm=T),
        estK=mean(estK, na.rm=T),
        diffK=mean(diffK, na.rm=T)) %>%
    ungroup()

## replications
df_fac %>%
    mutate(method_if = (method != "dipmeans")) %>%
    dplyr::select(method_if, p, n, omega, K, cov_type) %>%
    tbl_summary(by=method_if) %>%
    add_overall() %>%
    bold_labels()

## descriptive statistics
df_fac %>%
    mutate(diffK_r=diffK / (as.numeric(as.character(K))-1)) %>%
    dplyr::select(method, ARI, diffK_r) %>%
    pivot_longer(-method, names_to="variable") %>%
    group_by(method, variable) %>%
    summarize(
        Min=min(value),
        Q1=quantile(value, .25),
        Median=median(value),
        Q3=quantile(value, .75),
        Max=max(value),
        ) %>%
    ungroup() %>%
    gt(rowname_col = "method", groupname_col = "variable") %>%
    fmt_auto() # %>%
#    as_latex() %>%
#    cat()

### =======================================================================
###              Group Lasso analysis
### =======================================================================

y <- qnorm(df_fac$ARI) # probit transform

lm.s <- lm(y ~ method + p + n + omega + K + cov_type +
               method:(p + n + omega + K + cov_type)^3,
           weights = n_valid, data=df_fac)

summary(aov(lm.s)) # almost factors and interactions are significant in ANOVA

## =========== screening factor combinations using Group Lasso ===================
model <- lm.s

## extract design matrix
M <- model.matrix(model)
y <- model$fitted.values + model$residuals

## identify group of variables
tl <- attr(terms(model), "term.labels")
asgn <- attr(M, "assign")

## group lasso
fit <- grpreg::grpreg(M[,-1], y, asgn[-1])

## determine penalty parameter using EBIC
fit.selEBIC <- grpreg::select(fit, criterion = "EBIC")

## identify dropping variable groups
x <- fit.selEBIC$beta==0
zero_var <- tl[unique(asgn[x])]

## model update
chng <- paste("-", paste(zero_var, collapse=" - "))
refit <- update(model, paste("~.", chng))

## ============= extract ceofficients =====================
refit.s <- update(refit, ~.-n -method:K -method:cov_type)

object <- refit.s

tl <- attr(terms(object), "term.labels")
intr_order <- lengths(strsplit(tl, ":"))

print(tl)

refit.coef_mat <- get_rcoef(object)

## =====================================================
##               Visualization
## =====================================================

df_ord1 <- refit.coef_mat %>%
  filter(factor %in% tl[intr_order==2]) %>%
  mutate(method=factor(sapply(str_split(level, ":"), function(x) x[1]), levels=methods),
         item=factor(sapply(str_split(factor, ":"), function(x) x[2]), levels=item_names),
         value=as.numeric(sapply(str_split(level, ":"), function(x) x[2])))

df_ord2 <- refit.coef_mat %>%
  filter(factor %in% tl[intr_order==3]) %>%
  mutate(method=factor(sapply(str_split(level, ":"), function(x) x[1]), levels=methods),
         item1=factor(sapply(str_split(factor, ":"), function(x) x[2]), levels=item_names),
         item2=factor(sapply(str_split(factor, ":"), function(x) x[3]), levels=item_names),
         value1=as.numeric(sapply(str_split(level, ":"), function(x) x[2])),
         value2=as.numeric(sapply(str_split(level, ":"), function(x) x[3])))

## main effects
par(mar=c(4, 9, .1, .5))
n_len <- sum(refit.coef_mat$factor %in% tl[intr_order==1])
df_coef <- refit.coef_mat[1:(n_len+1),]
plot.interval_main(df_coef, xlim=c(-1, 1), cex.axis=1.4, line=4)

## two-factor interaction effects
df_ord1 %>%
  rename(coefficient=est) %>%
  ggplot(aes(x=value, y=coefficient, color=method)) +
  geom_point(size=1) +
  geom_line(lwd=0.5) +
  geom_errorbar(aes(ymin=coefficient-2*se, ymax=coefficient+2*se), width=0) +
  facet_grid(~item, scales="free_x") + ylim(-1, 1) +
  theme(text = element_text(size = 15))

## three-factor interaction effects
df_ord2 %>%
  filter(item1=="p") %>%
  rename(coefficient=est, value=value2) %>%
  filter(!(coefficient == 0 & se == 0)) %>%
  mutate(value1 = factor(value1, levels=c(2, 6, 18))) %>%
  ggplot(aes(x=value, y=coefficient, color=method)) +
  geom_point(size=1) + geom_line(lwd=0.5) +
  geom_errorbar(aes(ymin=coefficient-2*se, ymax=coefficient+2*se), width=0) +
  facet_grid(value1~item2, scales="free_x") +
  theme(text = element_text(size = 15))

## box plot
g1 <- df_fac %>%
  mutate(y=qnorm(ARI)) %>%
  ggplot(aes(x=n, y=y, color=method)) +
  geom_boxplot() + facet_grid(~p)
ggsave(filename="image/aggplot_method_p_n.png", plot=g1, width=9, height=3)

g2 <- df_fac %>%
  mutate(y=qnorm(ARI)) %>%
  ggplot(aes(x=K, y=y, color=method)) +
  geom_boxplot() + facet_grid(~p)
ggsave(filename="image/aggplot_method_p_K.png", plot=g2, width=9, height=3)

g3 <- df_fac %>%
  mutate(y=qnorm(ARI)) %>%
  ggplot(aes(x=cov_type, y=y, color=method)) +
  geom_boxplot() + facet_grid(~p)
ggsave(filename="image/aggplot_method_p_cov_type.png", plot=g3, width=9, height=3)
