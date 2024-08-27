#' Estimate hyperparameters for Duration Relational Event Model (DREM)
#' 
#' @description Function to estimate the hyperparameters for duration relational event model (DREM) using grid search
#' 
#'  
#' @param effects_start formula object for remstats, used to compute the start statistics
#' @param effects_end formula object for remstats, used to compute the end statistics
#' @param edgelist data.frame with columns (start_time, sender, receiver, end_time)
#' @param psi_start_candidates numeric vector of candidate values for psi_start. default value is 1
#' @param psi_end_candidates numeric vector of candidate values for psi_end. default value is 1
#' @param memory if "full" then no memory effects are incorporated. If "decay" then decay memory effects with specified half life
#' @param half_life_candidates numeric vector of candidate values for half life parameters
#' @param dur_undirected TRUE if riskset for the end DREM model needs to be undirected. See details
#' @param reh_undirected TRUE if riskset for both start and end DREM models needs to be undirected
#' @param strip_return logical, if TRUE then strip the heavy elements from a glm output object
#' @param save_dir character, local directory where to save fitted candidate model files
#' 
#' @return list with element \code{loglik}, a matrix or array of loglikelihood of fitted DREM candidate models and \code{mle}, a vector of candidate values (psi_start,psi_end(,half_life)) with maximum likelihood
#' @details
#' \code{dur_undirected} is set to \code{TRUE} if riskset for the duration model needs to be undirected. i.e if dyad A->B is in an event, the undirected dyad (AB==BA) is at risk to end the event. This argument can be used when it is not directly observable whether the sender or receiver ended the event
#' 
#' A list of available effects for the start and end models of DREM can be obtained with \code{\link[remstats:tie_effects]{remstats::tie_effects()}} and
#' for a list of undirected effects \code{\link[remstats:tie_effects]{remstats::tie_effects(directed = FALSE)}}
#' @examples
#' 
#' # Define effects for the start and end model of DREM
#' effects_start <- ~ 1 + remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
#' effects_end <- ~ 1 + remstats::outdegreeSender(scaling = "std")
#'
#' # Fit a DREM model
#' drem::dremstimate.grid(effects_start, effects_end, dat$edgelist)
#' 
#' @export
dremstimate.grid <- function(
    effects_start,
    effects_end,
    edgelist,    
    psi_start_candidates = 1,
    psi_end_candidates = 1,
    memory = c("full","decay"),
    half_life_candidates = NA,
    dur_undirected = FALSE,
    reh_undirected = FALSE,
    strip_return = TRUE,
    save_dir = NULL){
        
    if(any(edgelist$end_time < edgelist$start_time, na.rm = T)) {
        stop("End time of an event cannot be before it's start time.")
    }
    if(anyNA(edgelist$end_time)) {
        stop("Missing event end time")
    }
    
    memory <- match.arg(memory)
   
    if(memory=="decay" & any(is.na(half_life_candidates))){
        stop("Incorrect half life supplied for decay memory.")
    }

    if(memory == "decay"){
        loglik <- array(data = NA, dim = c(length(psi_start_candidates), length(psi_end_candidates),length(half_life_candidates)))
    }else{
        loglik <- array(data = NA, dim = c(length(psi_start_candidates), length(psi_end_candidates)))
    }

    # Progress bar
    pb <- txtProgressBar(min = 0, max = length(psi_start_candidates)*length(psi_end_candidates), style = 3)
    
    for(j in 1:length(psi_start_candidates)){
        for(k in 1:length(psi_end_candidates)){
            if(memory == "decay"){
                for(l in 1:length(half_life_candidates)){
                     fit <- dremstimate(effects_start, effects_end, edgelist,
                        psi_start = psi_start_candidates[j],
                        psi_end = psi_end_candidates[k],
                        memory = "decay",
                        half_life = half_life_candidates[l],
                        dur_undirected = dur_undirected,
                        reh_undirected = reh_undirected,
                        strip_return = strip_return)

                        loglik[j,k,l] = as.numeric(logLik(fit))
                        if(!is.null(save_dir)){
                            file_path = paste0(save_dir,"drem_fit_psi_start=",psi_start_candidates[j],"_psi_end=",psi_end_candidates[k],"_half_life=",half_life_candidates[l],".rdata")

                            save(fit,file = file_path)
                        }
                }
            }else{
                fit <- dremstimate(effects_start, effects_end, edgelist,
                        psi_start = psi_start_candidates[j],
                        psi_end = psi_end_candidates[k],
                        memory = "full",
                        dur_undirected = dur_undirected,
                        reh_undirected = reh_undirected,
                        strip_return = strip_return)
                loglik[j,k] = as.numeric(logLik(fit))

                if(!is.null(save_dir)){
                   file_path = paste0(save_dir,"drem_fit_psi_start=",psi_start_candidates[j],"_psi_end=",psi_end_candidates[k],".rdata")

                    save(fit,file = file_path)
                }
            }    
                   
            # Update progress
            setTxtProgressBar(pb, (j-1)*length(psi_start_candidates) + k)
        }    
    }

    dimnames(loglik)[[1]] = psi_start_candidates
    dimnames(loglik)[[2]] = psi_end_candidates
    
    indices <- which(loglik == max(loglik), arr.ind = TRUE)
    mle = c(dimnames(loglik)[[1]][indices[1]],dimnames(loglik)[[2]][indices[2]])

    if(memory=="decay"){
        dimnames(loglik)[[3]] = half_life_candidates
        mle = c(dimnames(loglik)[[1]][indices[1]],dimnames(loglik)[[2]][indices[2]],dimnames(loglik)[[3]][indices[3]])
    }
    
    return(list(loglik = loglik , mle = mle))
}


#' Estimate Duration Relational Event Model (DREM)
#' 
#' @description Function to estimate duration relational event model (DREM)
#' 
#'  
#' @param effects_start formula object for remstats, used to compute the start statistics
#' @param effects_end formula object for remstats, used to compute the end statistics
#' @param edgelist data.frame with columns (start_time, sender, receiver, end_time)
#' @param psi_start numeric value of psi parameter for start DREM model
#' @param psi_end numeric value of psi parameter for end DREM model
#' @param memory if "full" then no memory effects are incorporated. If "decay" then decay memory effects with specified half life
#' @param half_life numeric value of half life parameter for decay memory
#' @param dur_undirected TRUE if riskset for the end DREM model needs to be undirected. See details
#' @param reh_undirected TRUE if riskset for both start and end DREM models needs to be undirected
#' @param strip_return Logical, if TRUE then strip the heavy elements from a glm output object
#' @return fitted drem model object
#' 
#' @details
#' \code{dur_undirected} is set to \code{TRUE} if riskset for the duration model needs to be undirected. i.e if dyad A->B is in an event, the undirected dyad (AB==BA) is at risk to end the event. This argument can be used when it is not directly observable whether the sender or receiver ended the event.
#' 
#' A list of available effects for the start and end models of DREM can be obtained with \code{\link[remstats:tie_effects]{remstats::tie_effects()}} and
#' for a list of undirected effects \code{\link[remstats:tie_effects]{remstats::tie_effects(directed = FALSE)}}
#' 
#' The observation period is assumed to start at t = 0. However, the start time and end times must not be 0. If needed, add small random noise to the first start time
#' @examples
#' 
#' # Define effects for the start and end model of DREM
#' effects_start <- ~ 1 + remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
#' effects_end <- ~ 1 + remstats::outdegreeSender(scaling = "std")
#'
#' # Fit a DREM model
#' drem::dremstimate(effects_start, effects_end, edgelist)
#' 
#' @export
dremstimate <- function(
    effects_start,
    effects_end,
    edgelist,    
    psi_start = 1,
    psi_end = 1,
    memory = c("full","decay"),
    half_life = NA,
    dur_undirected = FALSE,
    reh_undirected = FALSE,
    strip_return = TRUE){

    if(any(edgelist$end_time < edgelist$start_time, na.rm = T)) {
        stop("End time of an event cannot be before it's start time.")
    }
    if(anyNA(edgelist$end_time)) {
        stop("Missing event end time")
    }
    
    memory <- match.arg(memory)
    if(memory=="decay" & is.na(half_life)){
        warning("No half life parameter supplied for decay memory. All events will have equal weight.")
        memory = "full"
    }
    
    stats_list = dremstats(effects_start, effects_end, edgelist,
                    psi_start = psi_start,
                    psi_end = psi_end,
                    memory = memory,
                    half_life = half_life,
                    dur_undirected = dur_undirected,
                    reh_undirected = reh_undirected)
    
    stat_stack <- drem_statstack(edgelist = edgelist, stats_list$start_stats, stats_list$end_stats, dur_undirected = dur_undirected,reh_undirected)   
    stat_names = colnames(stat_stack)[6:length(colnames(stat_stack))]
    
    fit <- glm(create_glm_formula(stat_names),
            family = poisson(),
            data = stat_stack)
    
    #Make the output smaller
    if(strip_return){
        fit <- strip_glm(fit)
    }

    return(fit)
    
}

# function to compute statistics for the duration rem using remstats 
# The stats generated from this function still require processing to be used to estimate the model
# 
# @param effects_start formula object for remstats, used to compute the start statistics
# @param effects_end formula object for remstats, used to compute the end statistics
# @param edgelist data.frame with columns in order: (start_time, sender, receiver, end_time)
# @param psi_start numeric value of psi parameter for start DREM model
# @param psi_end numeric value of psi parameter for end DREM model
# @param memory if "full" then no memory effects are incorporated. If "decay" then decay memory effects with specified half life
# @param half_life numeric value of half life parameter for decay memory
# @param dur_undirected TRUE if riskset for the end DREM model needs to be undirected. See details
# @param reh_undirected TRUE if riskset for both start and end DREM models needs to be undirected
# 
# @details
# dur_undirected is set to TRUE if riskset for the duration model needs to be undirected. i.e if dyad A->B is in an event, the undirected dyad (AB==BA) is at risk to end the event. This argument can be used when it is not directly observable whether the sender or receiver ended the event
# 
#@export 
dremstats <- function(effects_start, 
    effects_end, 
    edgelist, 
    psi_start = 0,
    psi_end = 0,
    memory = c("full","decay"),
    half_life = NA,
    dur_undirected = FALSE,
    reh_undirected = FALSE){
    
    colnames(edgelist)[1:4] = c("start_time","sender","receiver","end_time")
    
    #### start model
    edgelist$weight <- (edgelist$end_time - edgelist$start_time)^psi_start
    
    #split the edgelist into start and end event types
    dur.edgelist = edgelist[rep(seq_len(nrow(edgelist)), each = 2), ]
    dur.edgelist$type = rep(c("start","end"),nrow(edgelist))
    dur.edgelist$time = ifelse(dur.edgelist$type=="start",dur.edgelist$start_time,dur.edgelist$end_time)
    
    #memory
    memory <- match.arg(memory)
    memory <- ifelse(is.na(half_life), "full", "decay")
    
    #remify
    suppressWarnings({
        reh_dir <- remify::remify(dur.edgelist[,c("time","sender","receiver","weight","type")], directed = !reh_undirected, types = c("start","end") )
    })
    
    start_stats <- remstats::tomstats(effects = effects_start, reh = reh_dir, memory, half_life)
    dimnames(start_stats)[[3]] = paste0(dimnames(start_stats)[[3]],".start")
    
    #### end model
    edgelist$weight <- (edgelist$end_time - edgelist$start_time)^psi_end
    dur.edgelist = edgelist[rep(seq_len(nrow(edgelist)), each = 2), ]
    dur.edgelist$type = rep(c("start","end"),nrow(edgelist))    
    dur.edgelist$time = ifelse(dur.edgelist$type=="start",dur.edgelist$start_time,dur.edgelist$end_time)

    #end model can be either directed or undirected
    suppressWarnings({         
        reh_end <- remify::remify(dur.edgelist[,c("time","sender","receiver","weight","type")], types = c("start","end"), directed = !dur_undirected)
    })    
 
    end_stats <- remstats::tomstats(effects = effects_end, reh = reh_end, memory, half_life)

    dimnames(end_stats)[[3]] = paste0(dimnames(end_stats)[[3]],".end")

    return(list(start_stats = start_stats, end_stats = end_stats))
}


#'@keywords internal
create_glm_formula <- function(stat_names){
    return(paste0(" obs ~ -1 + offset(logtimediff) + ",paste(stat_names,collapse = " + ")))
}



# function to arrange remstats array into a stat stack for glm for the duration rem model
# 
# @param edgelist data.frame with columns in order (start_time, sender, receiver, end_time)
# @param start_stats remstats object for start stats
# @param end_stats remstats object for end stats
# @param dur_undirected (default=FALSE) boolean for if the end model should be undirected
# 
# @details
# start_stats has dimension (M,2D,P1) where M is the number of unique time points (over start and end of an event combined)
# D = total number of dyads in the network = N*(N-1)
# P1 = number of start statistics specified
# 
# end_stats has dimension (M,2D,P2) where P2 = number of end statistics specified
# 
# This function allows multiple start and end events to occur simultaneously, but only once per unique dyad in an interval
# @keywords internal
#'@export 
drem_statstack <- function(edgelist, start_stats, end_stats, dur_undirected = FALSE,reh_undirected = FALSE){
    
    colnames(edgelist)[1:4] = c("start_time","sender","receiver","end_time")
    rs = attr(start_stats,"riskset")
    #start and end event type dyad id
    if(!reh_undirected){
        edgelist$start_dyad <- apply(edgelist,1,function(x){
        return(which(rs[,1]== x[2] & rs[,2]==x[3] & rs[,3] == "start"))
        })
        edgelist$end_dyad <- apply(edgelist,1,function(x){
            return(which(rs[,1]== x[2] & rs[,2]==x[3] & rs[,3] == "end"))
        })        
    }else{
        edgelist$start_dyad <- apply(edgelist,1,function(x){
        ind1 = which(rs[,1]== x[2] & rs[,2]==x[3] & rs[,3] == "start")
        ind2 = which(rs[,2]== x[2] & rs[,1]==x[3] & rs[,3] == "start")
        return(max(ind1,ind2,na.rm = TRUE))
        })
        edgelist$end_dyad <- apply(edgelist,1,function(x){
            ind1 = which(rs[,1]== x[2] & rs[,2]==x[3] & rs[,3] == "end")
            ind2 = which(rs[,2]== x[2] & rs[,1]==x[3] & rs[,3] == "end")
            return(max(ind1,ind2,na.rm = TRUE))
        })        
    }

    if(dur_undirected & !reh_undirected){
        dur.rs = attr(end_stats,"riskset")
        edgelist$end_dyad_und <- apply(edgelist,1,function(x){
            ind1 = which(dur.rs[,2]== as.numeric(x[2]) & dur.rs[,1]==as.numeric(x[3]) & dur.rs[,3] == "end")
            ind2 = which(dur.rs[,1]== as.numeric(x[2]) & dur.rs[,2]==as.numeric(x[3]) & dur.rs[,3] == "end")
            return(max(ind1,ind2,na.rm = TRUE))
        })
    }
    
    P_start = dim(start_stats)[3]
    P_end  = dim(end_stats)[3]
    
    start_dyads = which(rs[,"type"]== "start")
    
    unique_times <- sort(unique(c(edgelist$start_time, edgelist$end_time)))
    
    #number of intervals
    M <- length(unique_times)
    
    logtimediff  <- log(unique_times - c(0,unique_times[-M]))
    
    stat_stack <- do.call(rbind.data.frame, lapply(1:M, function(i){
        stack_row = data.frame(matrix(0,nrow = 0,ncol = (5 + dim(start_stats)[3]+ dim(end_stats)[3])))
        
        ###### type 1::State: an event endped
        # end(d) is observed to end
        #Also allows instantaneous events because start_time <= t
        if(i > 1){
            observed_events = subset(edgelist,start_time <= unique_times[i] & end_time == unique_times[i])
            if(dur_undirected & !reh_undirected){
                dyads_observed = unique(observed_events[,"end_dyad_und"])
            }else{
                dyads_observed = unique(observed_events[,"end_dyad"])    
            }
            if(length(dyads_observed) > 0){
                stack_row = rbind(stack_row,as.data.frame(do.call(rbind,lapply(dyads_observed, function(d){
                    return(c(1, d, i, logtimediff[i],1, rep(0,P_start),unname(end_stats[i,d,])))
                }))))
            }    
        }

        ###### type 2::State: event ongoing
        # end(d) is at risk to end
        #all events active right now but didnt start in the interval
        active_events = subset(edgelist,start_time < unique_times[i] & end_time > unique_times[i])
        #end dyads currently in an event 
        if(dur_undirected& !reh_undirected){
            dyads_at_risk_end = unique(active_events[,"end_dyad_und"])
        }else{
            dyads_at_risk_end = unique(active_events[,"end_dyad"])
        }
        
        #add the active events with end type to riskset
        if(length(dyads_at_risk_end>0)){
            stack_row = rbind(stack_row,as.data.frame(do.call(rbind,lapply(dyads_at_risk_end, function(d){
                return(c(0, d, i, logtimediff[i] ,2,rep(0,P_start),unname(end_stats[i,d,])))
            }))))
        }
    
        ###### type 3::State: event just started
        #start(d) is observed
        #all events started in this interval, may also end in the same time => point events
        active_events = subset(edgelist,start_time == unique_times[i] & end_time >= unique_times[i])
        #end dyads currently in an event 
        dyads_at_risk_end = unique(active_events[,"start_dyad"])
        
        #add the active events with end type to riskset
        if(length(dyads_at_risk_end>0)){
            stack_row = rbind(stack_row,as.data.frame(do.call(rbind,lapply(dyads_at_risk_end, function(d){
                return(c(1, d, i, logtimediff[i] ,3, unname(start_stats[i,d,]),rep(0,P_end)))
            }))))
        }

        ###### type 4::State: event at risk to start
        # Set difference of all dyads that can start and the dyads that are currently active
        # start time can be current time for active event i.e event just started
        # end time cannot be current time i.e event just ended because then it would be at risk to start
        active_events = subset(edgelist,start_time <= unique_times[i] & end_time > unique_times[i])
        active_dyads = unique(active_events["start_dyad"])
        dyads_at_risk_start = setdiff(start_dyads,active_dyads)
        #add the active events with start type to riskset
        if(length(dyads_at_risk_start)>0){
            stack_row = rbind(stack_row,as.data.frame(do.call(rbind,lapply(dyads_at_risk_start, function(d){
                return(c(0, d, i, logtimediff[i] ,4,unname(start_stats[i,d,]),rep(0,P_end)))
            }))))            
        }
        return(stack_row)
    }))
    
    colnames(stat_stack) = c(c("obs","tie","i","logtimediff","type"),dimnames(start_stats)[[3]],dimnames(end_stats)[[3]])
    return(stat_stack)
}


strip_glm = function(cm) {
    cm$y = c()
    cm$model = c()
    
    cm$residuals = c()
    cm$fitted.values = c()
    cm$effects = c()
    #cm$qr$qr = c()
    #cm$linear.predictors = c()
    cm$weights = c()
    # cm$prior.weights = c()
    cm$data = c()
    
    
    cm$family$variance = c()
    #cm$family$dev.resids = c()
    cm$family$aic = c()
    cm$family$validmu = c()
    cm$family$simulate = c()
    attr(cm$terms,".Environment") = c()
    attr(cm$formula,".Environment") = c()
    
    cm
}