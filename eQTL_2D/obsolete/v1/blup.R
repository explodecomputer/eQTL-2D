
calc.resid <- function(formula, dat, Ainv, alpha)
{
	X <- Matrix(model.matrix(formula, dat))
	Y <- Matrix(model.frame(formula, dat)[, 1])
	Z <- Matrix(model.matrix(dat[, all.vars(formula)[1]] ~ dat[,1]))

	xtx <- crossprod(X)
	xtz <- crossprod(X, Z)
	ztx <- crossprod(Z, X)
	ztzainv <- crossprod(Z) + alpha * Ainv

	nr <- nrow(xtx) + nrow(ztx)
	LHS <- Matrix(0, ncol = nr, nrow = nr)
	LHS[1:nrow(xtx), 1:ncol(xtx)] <- xtx
	LHS[nrow(xtx) + seq(1:nrow(ztx)), 1:ncol(ztx)] <- ztx
	LHS[1:nrow(xtz), ncol(xtx) + seq(1, ncol(xtz))] <- xtz
	LHS[nrow(xtz) + seq(1:nrow(ztzainv)), ncol(ztx) + seq(1, ncol(ztzainv))] <- ztzainv
	RHS <- Matrix(0, ncol = 1, nrow = nr)
	xty <- crossprod(X, Y)
	RHS[1:nrow(xty), 1] <- xty
	zty <- crossprod(Z, Y)
	RHS[nrow(xty) + seq(1, nrow(zty))] <- zty
	sol <- solve(LHS, RHS)

	return(Y - sol)

}


blup2 <- function (formula, ped, alpha, Ainv, trim = FALSE) 
{
    colnames(ped)[1:3] <- c("ID", "SIRE", "DAM")
    if (any(!all.vars(formula) %in% colnames(ped))) 
        stop("Formula has variables which are not present in the data.")
    ww <- match(all.vars(formula)[-1], colnames(ped))
    ped$b <- apply(ped, 1, function(x) !any(is.na(x[ww])))
    if (trim) 
        ped <- ped[trimPed(ped, ped$b), ]
    ped$ID <- factor(ped$ID)
    ped$SIRE <- factor(ped$SIRE)
    ped$DAM <- factor(ped$DAM)
    X <- Matrix(model.matrix(formula, ped))
    Y <- Matrix(model.frame(formula, ped)[, 1])
    Z <- Matrix(model.matrix(ped[, all.vars(formula)[1]] ~ ped$ID))
    xtx <- crossprod(X)
    xtz <- crossprod(X, Z)
    ztx <- crossprod(Z, X)
    ztzainv <- crossprod(Z) + alpha * Ainv
    nr <- nrow(xtx) + nrow(ztx)
    LHS <- Matrix(0, ncol = nr, nrow = nr)
    LHS[1:nrow(xtx), 1:ncol(xtx)] <- xtx
    LHS[nrow(xtx) + seq(1:nrow(ztx)), 1:ncol(ztx)] <- ztx
    LHS[1:nrow(xtz), ncol(xtx) + seq(1, ncol(xtz))] <- xtz
    LHS[nrow(xtz) + seq(1:nrow(ztzainv)), ncol(ztx) + seq(1, 
        ncol(ztzainv))] <- ztzainv
    RHS <- Matrix(0, ncol = 1, nrow = nr)
    xty <- crossprod(X, Y)
    RHS[1:nrow(xty), 1] <- xty
    zty <- crossprod(Z, Y)
    RHS[nrow(xty) + seq(1, nrow(zty))] <- zty
    sol <- solve(LHS, RHS)
    row.names(sol) <- c(colnames(X), as.character(ped$ID))
    as.numeric(Y - sol[-c(1:ncol(X)), 1])
}

hsqs <- seq(0.1, 0.9, 0.1)
blups <- matrix(0, nrow(test), length(hsqs))
for(i in 1:length(hsqs))
{
	blups[,i] <- blup2(test)
	
	
	
}
