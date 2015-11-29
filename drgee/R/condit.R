condit <-
    function(y, x.cent, id, param.names) {

        ## browser()
        ## Assuming that the observations are sorted by id

        ## Find the total number of observations
        n.obs <- length(id)

        ## Find the number of parameters
        n.params <- ncol(x.cent)
        
        ## Find cluster sums
        clust.size <- ave(y, id, FUN=length)
        
        ## Find cluster sums
        ysum <- ave(y, id, FUN=sum)

        ## Find cluster with more than half being cases
        inv <- ifelse(ysum > 0.5*clust.size, 1, 0)

        ## Find discordant clusters
        disc = ifelse(0 < ysum & ysum < clust.size, 1, 0)

        ## Find the indices corresponding to discordant clusters
        disc.idx <- which(disc == 1)

        ## Find the cluster ids to use
        clusterid.disc <- unique(id[disc.idx])

        ## Find start indices for the discordant clusters
        min.idx.disc <- match( clusterid.disc, id )
        
        ## Find cluster sums for the discordant clusters
        ysum.disc <- ysum[min.idx.disc]

        ## Find cluster sizes for the discordant clusters
        clust.size.disc <- clust.size[min.idx.disc]

        ## Find clusters with more than half being cases
        ## for the discordant clusters
        inv.disc <- inv[min.idx.disc]

        ## Fit using clogit
        fit.clogit <- try(survival::clogit(y~x.cent + strata(id), method =
        "exact", subset = disc.idx) )

        if (class(fit.clogit)[1] == "try-error") {

            beta.hat <- rep(NA, n.params)

            eq.res <- rep(NA, n.obs)

            d.eq.res <- matrix( rep(NA, n.obs*n.params), nrow = n.obs )

            naive.var = NULL
            
        } else {
                
            beta.hat <- coef(fit.clogit)

            resids <- .Call("conditRes", beta.hat, ysum.disc,
                            clust.size.disc, min.idx.disc, inv.disc,
                            y, x.cent, PACKAGE = "drgee")

            eq.res <- as.vector(resids$res)

            d.eq.res <- resids$dres

            naive.var = vcov(fit.clogit)
                
        }


        
        names(beta.hat) <- param.names
        colnames(d.eq.res) <- param.names

        return( list(coefficients = beta.hat,
                     res = eq.res,
                     d.res = d.eq.res,
                     eq.x = x.cent,
                     optim.object = NULL,
                     naive.var = naive.var) )

    }
