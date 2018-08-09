#' Title
#'
#' @param sp
#' @param dfvec
#'
#' @return
#' @export
#'
#' @examples

indsp.func <- function(sp, dfvec = c(4, 7, 10, 15, 20, 33)){
        # indsp.func:
        # Calculates the data matrices of indices for different df.
        # Arguments are the species data frame sp,
        # with column headers marked "site", "year", and "count",
        # and a list of required degrees of freedom (could be just one) for
        # a selection of GAM fits.
        #
        # All rows with NA's (missing values) must be removed from the data before
        # applying this function.
        # Result is a matrix with columns giving indices for each of the Nyears
        # years, with df as marked at the top of the column.
        #
        # example: species data frame is called cb (for common bird):
        # at command line:   indcb <- indsp.func(cb, c(4, 7, 10, 15, 20, 33))
        #
        # NOTE: SOME PARTS OF THIS FUNCTION MAY BE SPECIFIC TO THE GAM
        # FITTED ON count ~ as.factor(site) + s(year).
        # IF USING A DIFFERENT FORMULA, NEED TO MODIFY CODE: ESPECIALLY SEE
        # NOTE MARKED *** BELOW.
        #
        if(length(sp$count[is.na(sp$count)]) > 0) stop(
                "no missing data allowed")

        # fit.func fits the GAM using the MGCV library and extracts the indices for a given
        # value of df:
        fit.func <- function(dfval)
        {
                gam.df <- mgcv::gam(count ~ as.factor(site) + s(year, fx=TRUE, k = dfval + 1),
                              family = poisson(link = log), data = sp)
                pred.df <- mgcv::predict.gam(gam.df,type="terms")
                srow <- length(pred.df[1,])
                # *** For old MGCV versions, change line above to
                #     srow <- length(pred.df[,1])
                # srow is the row containing the smooth term in the predict object
                #
                # *** IF CHANGING THE ORDER OF VARIABLES IN FORMULA, OR FITTING GAM
                # WITH COVARIATES, NEED TO CHECK WHETHER IT IS STILL THE SAME ROW REQUIRED ***
                #
                Nentries <- length(sp$count)
                reqd.entries <- (1:Nentries)[!duplicated(sp$year)] #
                # reqd.entries gives the entries corresponding to the first occurrence of each
                # year in the data frame sp.  The years are not in order, but given by
                # years.entries.  Both are then ordered before being returned to ind.df.
                years.entries <- sp$year[!duplicated(sp$year)]
                years.pred <- pred.df[reqd.entries, srow]
                # *** For old MGCV versions, change line above to
                #     years.pred <- pred.df[srow, reqd.entries]
                years.ord <- sort(years.entries)
                pred.ord <- years.pred[order(years.entries)]
                ind.df <- exp(pred.ord)/exp(pred.ord[1])        #

                # indices are scaled around the base year: chosen here as year 1.
                cat("df ", dfval, " complete \n")
                ind.df
        }
        ind.all <- sapply(dfvec, fit.func)
        dimnames(ind.all) <- list(character(0), paste("df", as.character(dfvec),
                                                      sep = ""))
        ind.all
}

#' Title
#'
#' @param sp
#' @param defree
#'
#' @return
#' @export
#'
#' @examples
bootstrap.func <- function(sp, defree = 10)
{
        # bootstrap.func:
        # Outputs the index curve for a single bootstrap replicate.
        # Arguments are the species data frame sp,
        # with column headers marked "site", "year", and "count";
        # and the required degrees of freedom for the fit, defree.
        #
        # example: species data frame is called cb (for common bird):
        # call:  bootstrap.func(cb, 10)
        # but generally this function is only called from within the outer
        # bootstrap routine, outer.boot.func.
        #
        # NOTE: SOME PARTS OF THIS FUNCTION MAY BE SPECIFIC TO THE GAM
        # FITTED ON count ~ as.factor(site) + s(year).
        # IF USING A DIFFERENT FORMULA, NEED TO MODIFY CODE: ESPECIALLY SEE
        # NOTE MARKED *** BELOW.
        #
        # first obtain the resampled data:
        uniq.site <- unique(sp$site)
        Nentries <- length(sp$count)
        Nyears <- length(unique(sp$year))
        Nsites <- length(uniq.site)	#
        # take sample of sites to be included in the resample
        # (bootstrap across sites):
        sam <- sample(uniq.site, replace = T)	#

        # elements.func lists the rows of the data frame sp that need to
        # be included in the resample, for an individual site x in the sample sam:
        elements.func <- function(x)
        {
                (1:Nentries)[sp$site == x]
        }
        elements <- lapply(sam, elements.func)	#
        # elements is the set of rows of the data frame to be included in the
        # resample (with repetitions listed separately).  It is a ragged list.
        dr.site <- rep(1:Nsites, as.vector(unlist(lapply(elements, length))))	#
        # dr.site is the vector of site levels for the data resample: these should
        # go from 1 to Nsites even though some sites have been included in the
        # resample more than once, because we want to fit the model to the
        # repetitions as if they were different sites.  However, year and count
        # for the resample should be extracted directly from the original data
        # frame, sp.
        #
        elements <- as.vector(unlist(elements))
        data.resample <- data.frame(site = dr.site, year = sp$year[elements],
                                    count = sp$count[elements])	#

        # fit GAM on this df using MGCV library:
        dr.gam <- mgcv::gam(count~as.factor(site) + s(year, fx=TRUE, k=defree+1), family =
                              poisson(link = log), data = data.resample)	#
        # dr.gam is the GAM on the data resample.
        #
        # Now make out a new data frame consisting of all years (in case some were
        # omitted from data.resample) and all sites in data.resample.
        # The predict feature doesn't seem to work unless all sites are included.

        year.base <- min(sp$year)-1
        dr.new <- expand.grid(year=year.base+1:Nyears, site=sort(unique(dr.site)))
        result <- mgcv::predict.gam(dr.gam, newdata=dr.new, type = "terms")
        srow <- length(result[1,])
        # *** For old MGCV versions, change line above to
        #     srow <- length(result[,1])
        #
        # srow is the row containing the smooth term in the predict object
        #
        # *** IF CHANGING THE ORDER OF VARIABLES IN FORMULA, OR FITTING GAM
        # WITH COVARIATES, NEED TO CHECK WHETHER IT IS STILL THE SAME ROW REQUIRED ***
        #
        pred.allyears <-  result[1:Nyears, srow]
        # *** For old MGCV versions, change line above to
        #     pred.allyears <-  result[srow, 1:Nyears]
        indices <- exp(pred.allyears)/exp(pred.allyears[1])        #
        indices


}

#' Title
#'
#' @param sp
#' @param defree
#' @param Nrep
#' @param indexfile
#'
#' @return
#' @export
#'
#' @examples
outer.boot.func <- function(sp, defree = 10, Nrep = 399)
{
        # outer.boot.func:
        # outer function for bootstrapping.
        # Calls the function bootstrap.func that will compute a
        # bootstrap resample, fit the GAM and find the indices.
        #
        # sp is the dataframe with columns headed "site", "year", "count",
        # and no missing values.
        # Nrep is number of bootstrap replicates; defree is degrees of freedom
        # for the GAM fit.  Indexfile is the name of a file in the working
        # directory (see below): it should be in quotes (" ").
        #
        # Example:
        #   cb.bootind.399.10 <- outer.boot.func(cb, 10, 399, "cb.bootind.399.10")
        #
        # NOTE:
        # The file indexfile is a safeguard against premature termination of the
        # bootstrap routine: all results are written to this file as the
        # function proceeds.  If the program crashes, the early results can be
        # salvaged by:
        #   sp.bootind.crash <- matrix(scan("indexfile"), nrow = k, byrow = T)
        # where k is the number of replicates completed before the crash.
        # The remaining replicates can then be made up by restarting the
        # oter.boot.func function.
        #
        # ** IMPORTANT NOTE **
        # When restarting any Splus simulation after a crash, be sure to reset the
        # random seed before resuming: otherwise the previous results (before the
        # crash) will just be repeated.  The easiest way to reset the seed is
        # simply to generate some new random numbers: e.g. runif(100).
        # So (e.g.) if the procedure crashed after replicate 50 out of 399, the
        # full run and salvage operation would be something like:
        #    cb.bootind.399.10 <- outer.boot.func(cb, 10, 399,"cb.results")
        #        (this crashes and the object cb.bootind.399.10 is never formed)
        #    cb.bootind.crash50 <- matrix(scan("cb.results"), nrow=50, byrow=T)
        #    runif(100)
        #    cb.bootind.last349 <- outer.boot.func(cb, 10, 349, "cb.results2")
        #    cb.bootind.399.10 <- rbind(cb.bootind.crash50, cb.bootind.last349)
        #
        #
        # START OF FUNCTION
        # some housekeeping jobs first:

        # if(!is.character(indexfile)) stop(
        #         "filenames have to be character strings")	#
        # the dataset sp is assumed to have NO MISSING ENTRIES:
        if(length(sp$count[is.na(sp$count)]) > 0) stop(
                "no missing data allowed")	#
        # start bootstrap loop:
        replicate.func <- function(i)
        {
                print(i)
                rep.indices <- bootstrap.func(sp, defree = defree)
                # write(rep.indices, file = indexfile, ncolumns = 4, append = T)
                rep.indices
        }
        index.matrix <- sapply(1:Nrep, replicate.func)
        matrix(unlist(index.matrix), byrow = T, nrow = Nrep)
}

#' Title
#'
#' @param x
#' @param h
#' @param d
#' @param interval
#'
#' @return
#' @export
#'
#' @examples
diffop.func <- function(x, h, d, interval)
{
        # diffop.func:
        # Calculates the difference operator on vector x to give an
        # approximate second derivative curve for x.
        #
        # h is the desired window size (see below). Corresponds to the quantity r
        #    in eqns 9, 10, and 11 of the Ecology paper.
        #    Usually h will be chosen to be the smallest possible (h=1) but try
        #    experimenting with h=2 and h=3 from within function sp.plot,
        #    to see if it has much effect on the changepoints.
        # d is the order of the difference operator: can be 2, 4, or 6,
        #    representing the second derivative estimates obtained from 2nd
        #    differences, 4th differences, or 6th differences (eqns 9, 10, and 11
        #    in the Ecology paper).  d=6 should give the most accurate results,
        #    but again try experimenting with d=4 and d=2.
        #"interval" is the time elapsed between points of x, which will be treated
        #    as 1 unit:
        #    e.g. if x consists of records taken every 2 years, then 2 years is
        #    equal to 1 unit, and interval=2.  Usually interval=1.
        #
        # The window size h is measured relative to the interval variable,
        # so if the window extends h units then it covers (h*interval) of
        # actual time: e.g. if interval = 2 years, and h = 3 units, then actual
        # window size is 2*3=6 years).
        #   [The reason for this setup is because the function works with the
        #   quantity x[t+h], not x(t+h) --- i.e. h is the index of a vector
        #   element, not a measurement of time.  Note however that to obtain the
        #   2nd derivative, the h^2 on the denominator of eqns 9, 10, and 11 is
        #   supposed to be a measurement of time; so in this function we have to
        #   divide by (h*interval)^2 (see penultimate line of code where the
        #   adjustment by interval^2 is made), instead of just by h^2.]
        #
        # NOTE:
        # The quantity x corresponds to the vector of abundance indices
        # in the Ecology paper: x=(I(1), I(2), ...).
        #
        # EXAMPLE:
        # this function is usually called from within sp.plot via the
        # function indmat.derivmat.func, but here are some experimentations
        # with the values of d and h:
        #     plot(1:34, diffop.func(indcb[,"df10"], 1, 6, 1), type="l")
        #     abline(h=0)
        #     lines(1:34, diffop.func(indcb[,"df10"], 1, 4, 1), col=2)
        #     lines(1:34, diffop.func(indcb[,"df10"], 1, 2, 1), col=3)
        #     lines(1:34, diffop.func(indcb[,"df10"], 2, 6, 1), col=4)
        #     lines(1:34, diffop.func(indcb[,"df10"], 3, 6, 1), col=5)
        #     etc.
        #
        # first some housekeeping
        if(is.na(match(d, c(2, 4, 6)))) stop("d must be 2, 4 or 6.")	#
        # now work out the vector of window sizes:
        # at endpoints it is not possible to use d=4 or d=6, or h>1, because
        # there are insufficiently many points at the end of the series to
        # accommodate the larger windows required.  Thus whatever values were
        # specified for d and h, the endpoints of the series must have d=2 and
        # h=1.  This function maintains h at 1 until sufficiently far from the
        # endpoints to allow d to increase (if necessary) to its specified value,
        # after which h is also allowed to increase (if necessary) to its
        # specified value.  The actual window sizes used for every time point
        # in the series are stored in vector hvec:
        N <- length(x)	# x is recorded at N equally spaced points
        hvec <- 1:(N/2)
        hvec <- (hvec - 1) %/% (d/2)
        hvec[hvec > h] <- h
        if(N %% 2)
                hvec <- c(hvec, max(hvec), rev(hvec))	# (N odd)
        else hvec <- c(hvec, rev(hvec))
        hvec[hvec == 0] <- 1	#
        # similarly, dvec stores the actual values (2, 4, or 6) for the
        # difference operator d along the series:
        dvec <- rep(d, N)
        if(d == 6) {
                dvec[c(2, N - 1)] <- 2
                dvec[c(3, N - 2)] <- 4
        }
        else if(d == 4) {
                dvec[c(2, N - 1)] <- 2
        }
        dvec[c(1, N)] <- 0	#
        # Now code the various difference operators approximating 2nd derivative,
        # corresponding to d=2, 4, or 6.
        # xinds is the set of indices of x for which the derivative is to be
        # calculated (usually all those with d taking a certain value).
        # hvals are the h-values corresponding to those indices in xinds.
        diff.op2 <- function(xinds, hvals)
        {
                (x[xinds + hvals] - 2 * x[xinds] + x[xinds - hvals])/hvals^2
        }
        diff.op4 <- function(xinds, hvals)
        {
                ( - x[xinds + 2 * hvals] + 16 * x[xinds + hvals] - 30 * x[xinds
                                                                          ] + 16 * x[xinds - hvals] - x[xinds - 2 * hvals])/(12 *
                                                                                                                                     hvals^2)
        }
        diff.op6 <- function(xinds, hvals)
        {
                (2 * x[xinds + 3 * hvals] - 27 * x[xinds + 2 * hvals] + 270 * x[
                        xinds + hvals] - 490 * x[xinds] + 270 * x[xinds - hvals
                                                                  ] - 27 * x[xinds - 2 * hvals] + 2 * x[xinds - 3 * hvals
                                                                                                        ])/(180 * hvals^2)
        }
        ##
        secderiv <- numeric(N)
        secderiv[0] <- 0
        secderiv[N] <- 0	#
        # apply the relevant difference operators throughout the series to get the
        # estimate secderiv of the second derivative vector:
        secderiv[dvec == 2] <- diff.op2((1:N)[dvec == 2], hvec[dvec == 2])
        secderiv[dvec == 4] <- diff.op4((1:N)[dvec == 4], hvec[dvec == 4])
        secderiv[dvec == 6] <- diff.op6((1:N)[dvec == 6], hvec[dvec == 6])	#
        # correct for interval^2 as described above:
        secderiv <- secderiv/(interval^2)
        secderiv
}

#' Title
#'
#' @param sp.bootind
#' @param h
#' @param d
#' @param interval
#'
#' @return
#' @export
#'
#' @examples
indmat.derivmat.func <- function(sp.bootind, h, d, interval)
{
        # indmat.derivmat.func:
        # Calculates the matrix of 2nd derivative estimates from the
        # matrix sp.bootind of bootstrapped abundance indices.
        #
        # sp.bootind is a matrix of the form sp.bootind.Nrep.df
        # (e.g. cb.bootind.399.10).
        # Each row of the output matrix is the vector of 2nd derivative estimates
        # for the abundance index curve in the corresponding row of the index
        # matrix (i.e. the curve corresponding to a particular bootstrap
        # replicate).
        # See diffop.func for a full explanation of arguments  h, d, and interval.
        # interval is the time lapse between consecutive index points: e.g. if
        #    the index curve measures abundance every 2 years, then interval=2.
        #    Usually interval=1.
        # h is related to window size for the second derivative estimates: usually
        #    will be chosen to be the smallest possible (h=1) but try
        #    experimenting with h=2 and h=3 to see if the results are any
        #    different.
        # d is the order of the difference operator: can be 2, 4, or 6,
        #    representing the second derivative estimates obtained from 2nd
        #    differences, 4th differences, or 6th differences (eqns 9, 10, and 11
        #    in the Ecology paper).  d=6 should give the most accurate results.
        #
        # The output from this function would normally be named
        # sp.bootderiv.Nrep.df, although this function is usually only called
        # from within sp.plot.
        # example:
        # cb.bootderiv.399.10 <- indmat.derivmat.func(cb.bootind.399.10, 1, 6, 1)
        #
        #
        t(apply(sp.bootind, 1, diffop.func, h = h, d = d, interval = interval))
}


#' Title
#'
#' @param index.matrix
#' @param conf
#'
#' @return
#' @export
#'
#' @examples
index.ci.func <- function(index.matrix, conf = 0.95)
{
        # index.ci.func:
        # Given a matrix of bootstrapped index values of the form
        # index.matrix = sp.bootind.Nrep.df, this function extracts
        # the lower and upper (100*conf)% confidence limits for the abundance
        # indices.
        #
        # conf should be chosen so that (Nrep+1)*(1-conf)/2 is an integer:
        # otherwise it will not be possible to find integers l and u such that the
        # l'th and u'th ordered values give us exactly a conf% confidence interval.
        #
        # syntax:  sp.ind.ci <- index.ci.func(sp.bootind.Nrep.df, 0.95)
        # example:  cb.ind.ci <- index.ci.func(cb.bootind.399.10, 0.95)
        #
        # Usually called from within the function sp.plot.
        #
        Nyears <- ncol(index.matrix)
        Nrep <- nrow(index.matrix)	#
        # Nrep is the number of bootstrap replicates in the index matrix.
        # Find alp: alp=(1-conf)/2:
        alp <- (1 - conf)/2	#
        # Interpolation is not ideal, so discard some rows of the bootstrap
        # matrix if (Nrep+1)*alp is not an integer:
        if(abs((Nrep + 1) * alp - round((Nrep + 1) * alp)) > 1e-05)
                stop("need to discard rows of index.matrix, or change conf, so that (Nrep+1)*(1-conf)/2 is an integer."
                )

        lowerpt <- (Nrep + 1) * alp
        upperpt <- (Nrep + 1) * (1 - alp)

        # The confidence interval goes from the (Nrep+1)*alp 'th ordered
        # bootstrap value (low) to the (Nrep+1)*(1-alp) 'th ordered bootstrap
        # value (high).
        inner.func <- function(yr)
        {
                sort(index.matrix[, yr])[c(lowerpt, upperpt)]
        }
        index.ci <- sapply(seq(1, Nyears), inner.func)
        dimnames(index.ci) <- list(c("lower", "upper"), NULL)
        index.ci
}


#' Title
#'
#' @param deriv.matrix
#' @param conf
#'
#' @return
#' @export
#'
#' @examples
deriv.ci.func <- function(deriv.matrix, conf = 0.95)
{
        # deriv.ci.func:
        # Given a matrix of bootstrapped 2nd derivative values (of the form
        # deriv.matrix = sp.bootderiv.Nrep.df), this function extracts
        # the lower and upper (100*conf)% confidence limits for the second
        # derivatives.
        #
        # conf should be chosen so that (Nrep+1)*(1-conf)/2 is an integer:
        # otherwise it will not be possible to find integers l and u such that the
        # l'th and u'th ordered values give us exactly a conf% confidence interval.
        #
        # syntax:  sp.deriv.ci <- deriv.ci.func(sp.bootderiv.Nrep.df, 0.95)
        # example:  cb.deriv.ci <- deriv.ci.func(cb.bootderiv.399.10, 0.95)
        #
        # However, this function is rarely called directly: almost always called
        # from within the function sp.plot.
        #
        Nyears <- ncol(deriv.matrix)
        Nrep <- nrow(deriv.matrix)	#
        # Nrep is the number of bootstrap replicates in the 2nd derivative
        # matrix.
        # Find alp: alp=(1-conf)/2:
        alp <- (1 - conf)/2	#
        # Interpolation is not ideal, so discard some rows of the bootstrap
        # matrix if (Nrep+1)*alp is not an integer:
        if(abs((Nrep + 1) * alp - round((Nrep + 1) * alp)) > 1e-05)
                stop("need to discard rows of deriv.matrix, or change\n\tconf, so that (Nrep+1)*(1-conf)/2 is an integer."
                )
        lowerpt <- (Nrep + 1) * alp
        upperpt <- (Nrep + 1) * (1 - alp)

        # The confidence interval goes from the (Nrep+1)*alp 'th ordered
        # bootstrap value (low) to the (Nrep+1)*(1-alp) 'th ordered bootstrap
        # value (high).
        inner.func <- function(yr)
        {
                sort(deriv.matrix[, yr])[c(lowerpt, upperpt)]
        }
        deriv.ci <- sapply(seq(1, Nyears), inner.func)
        dimnames(deriv.ci) <- list(c("lower", "upper"), NULL)
        deriv.ci
}

#' Title
#'
#' @param sp
#' @param defree
#' @param sp.bootind
#' @param h
#' @param d
#' @param interval
#' @param conf
#' @param add
#' @param col
#' @param cex
#'
#' @return
#' @export
#'
#' @examples
sp.plot <- function(sp, defree, sp.bootind, h, d, interval, conf = 0.95,
                    add = F, col = 2, cex = 2)
{
        # sp.plot:
        # This function wraps everything up and performs some of the
        # "invisible" calculations (e.g. the second derivative estimates and
        # changepoints).
        # It plots the estimated abundance indices and confidence intervals, and
        # marks significant upturns and downturns (i.e. "changepoints") on the
        # curve in open and closed circles.
        #
        # "sp" is the name of the species data frame, and must be a character
        #    string: e.g. "cb".
        #    It is assumed that the R object "indsp" (e.g. indcb) containing
        #    the GAM fits to the original data for the required df has been
        #    obtained previously (using function indsp.func), and is available
        #    BY THAT NAME on the permanent database.
        #    If it is not there, then form it first at the command line by typing:
        #    indsp <- indsp.func(sp, dfvec=defree)
        #    (e.g. indcb <- indsp.func(cb, defree) where defree is replaced by
        #    the required value).
        # "defree" is the df of the required GAM fit: e.g. if defree=10 then
        #    the indices indsp[,"df10"] will be extracted for the main curve.
        # "sp.bootind" is a matrix of bootstrapped indices, e.g.
        #    cb.bootind.399.10.  It should be the relevant object for the GAM
        #    with specified df (i.e. if df=10 then sp.bootind should be
        #    sp.bootind.Nrep.10 for some Nrep).
        #
        # See function diffop.func for explanation of "h", "d", and "interval".
        # "conf" gives confidence level for the ABUNDANCE INDICES confidence
        #    intervals: *not* the confidence level for the changepoints.
        #    The changepoints are automatically calculated at the 5% significance
        #    level. This can be changed if required by manually editing the value
        #    "deriv.conf" in the function body.
        # "add" determines whether to add the plot to a current plot, or
        #    whether to start a new plot.  The default is to start a new plot.
        # "col": if add=T, col determines what colour the new plot (superposed on
        #    the existing plot) should be plotted in.
        # "cex" determines the size of the text and changepoints on the plot:
        #    increase it for larger text, decrease it for smaller text.
        #
        # EXAMPLE:
        # # plot for 10 df:
        #      sp.plot("cb", 10, cb.bootind.399.10, 1, 6, 1, conf=0.95)
        # # then superpose the plot for 20 df:
        #      sp.plot("cb", 20, cb.bootind.399.20, 1, 6, 1, conf=0.95, add=T)
        #
        deriv.conf <- 0.95
        par(font = 3)	#
        # housekeeping:
        if(!is.character(sp))
                stop("first argument must be a character string.")
        dfname <- paste("df", defree, sep = "")
        sp.ind.ci <- index.ci.func(sp.bootind, conf = conf)
        Nyears <- length(sp.ind.ci[1,  ])
        yearvec <- seq(1, Nyears)	#
        # Might want to change yearvec to the actual years of the survey: e.g.
        # yearvec <- seq(1962, 1995) instead of seq(1, 34).
        #
        lower <- min(sp.ind.ci[1,  ])
        upper <- max(sp.ind.ci[2,  ])
        indsp <- eval(parse(text = paste("ind", sp, sep = "")))
        if(!add) {
                plot(yearvec, indsp[, dfname], type = "l", ylim = c(lower,
                                                                    upper), xlab = "", ylab = "", cex = cex, lwd = 2, bty
                     = "l", font = 3)
                textint <- 0.2/1.7 * (upper - lower)	#
                # textint is for beautification and nothing else
                text(min(yearvec) - 1, upper + textint, "index", cex = cex,
                     font = 3)
                text(max(yearvec), lower - textint, "year", cex = cex, font = 3
                )
        }
        if(add) {
                par(col = col)
                lines(yearvec, indsp[, dfname], lwd = 2)
        }
        lines(yearvec, sp.ind.ci[1,  ], lty = 4, lwd = 2)
        lines(yearvec, sp.ind.ci[2,  ], lty = 4, lwd = 2)	#
        # Now calculate matrix of bootstrapped second derivatives:
        sp.bootderiv <- indmat.derivmat.func(sp.bootind = sp.bootind, h = h, d
                                             = d, interval = interval)
        sp.deriv.ci <- deriv.ci.func(sp.bootderiv, conf = deriv.conf)	#
        # Calculate the changepoints:
        # the lower changepoints are those where the *upturn* of the curve
        # is statistically significant (i.e. the lower confidence limit is
        # greater than 0);
        # the upper changepoints are those where the *downturn* of the curve
        # is statistically significant (i.e. the upper confidence limit is
        # less than 0);
        changepts.lower <- seq(1, Nyears)[sp.deriv.ci["lower",  ] > 0]
        changepts.upper <- seq(1, Nyears)[sp.deriv.ci["upper",  ] < 0]
        if((length(changepts.lower > 0) & (!is.na(changepts.lower[1])))) {
                points(changepts.lower + min(yearvec) - 1, indsp[
                        changepts.lower, dfname], pch = 1, cex = cex)
                changepts.lower <- changepts.lower + min(yearvec) - 1
        }
        else changepts.lower <- "none"
        if((length(changepts.upper > 0) & (!is.na(changepts.upper[1])))) {
                points(changepts.upper + min(yearvec) - 1, indsp[
                        changepts.upper, dfname], pch = 16, cex = cex)
                changepts.upper <- changepts.upper + min(yearvec) - 1
        }
        else changepts.upper <- "none"
        par(col = 1)
        list(upturns = changepts.lower, downturns = changepts.upper)
}

