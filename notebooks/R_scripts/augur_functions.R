# adjusted augur variance filter that uses dispersion per gene
select_variance_disp <- function (mat, var_quantile = 0.5, filter_negative_residuals = FALSE) {
    sds = rowSds(mat)
    sds[is.na(sds)] <- 0
    mat %<>% extract(sds > 0, )
    if (var_quantile < 1 | filter_negative_residuals) {
        means = rowMeans(mat)
        sds %<>% extract(. > 0)
        # # old method
        # cvs = means/sds
        # new calculation for dispersion
        cvs = sds^2 / means
        lower = quantile(cvs, 0.01)
        upper = quantile(cvs, 0.99)
        keep = between(cvs, lower, upper)
        cv0 = cvs[keep]
        mean0 = means[keep]
        if (any(mean0 < 0)) {
            model = loess(cv0 ~ mean0)
        }
        else {
            fit1 = loess(cv0 ~ mean0)
            fit2 = loess(cv0 ~ log(mean0))
            cox = coxtest(fit1, fit2)
            probs = cox$`Pr(>|z|)`
            if (probs[1] < probs[2]) {
                model = fit1
            }
            else {
                model = fit2
            }
        }
        genes = rownames(mat)[keep]
        residuals = setNames(model$residuals, genes)
        if (filter_negative_residuals == T) {
            genes = names(residuals)[residuals > 0]
        }
        else {
            genes = names(residuals)[residuals > quantile(residuals, 
                var_quantile)]
        }
        mat %<>% extract(genes, )
    }
    return(mat)
}
