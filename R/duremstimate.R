#' Estimate hyperparameters for Duration Relational Event Model (DuREM)
#' 
#' @description Function to estimate the hyperparameters for duration relational event model (DuREM) using grid search
#' 
#'  
#' @param start_effects formula object for remstats, used to compute the start statistics
#' @param end_effects formula object for remstats, used to compute the end statistics
#' @param edgelist data.frame with columns (start_time, sender, receiver, end_time)
#' @param psi_start_candidates numeric vector of candidate values for psi_start. default value is 1
#' @param psi_end_candidates numeric vector of candidate values for psi_end. default value is 1
#' @param start_undirected Logical. If `TRUE`, the risk set for start DuREM model is undirected. Default is `FALSE`.
#' @param end_undirected Logical. If `TRUE`, the risk set for the end DuREM model is undirected (see Details). Default is `FALSE`.
#' @param memory if "full" then no memory effects are incorporated. If "decay" then decay memory effects with specified half life
#' @param half_life_candidates numeric vector of candidate values for half life parameters
#' @param start_undirected TRUE if riskset for both start and end DuREM models needs to be undirected
#' @param strip_return logical, if TRUE then strip the heavy elements from a glm output object
#' @param save_dir character, local directory where to save fitted candidate model files
#' 
#' @return list with element \code{loglik}, a matrix or array of loglikelihood of fitted DuREM candidate models and \code{mle}, a vector of candidate values (psi_start,psi_end(,half_life)) with maximum likelihood
#' @details
#' \code{end_undirected} is set to \code{TRUE} if riskset for the duration model needs to be undirected. i.e if dyad A->B is in an event, the undirected dyad (AB==BA) is at risk to end the event. This argument can be used when it is not directly observable whether the sender or receiver ended the event
#' 
#' A list of available effects for the start and end models of DuREM can be obtained with \code{\link[remstats:tie_effects]{remstats::tie_effects()}} and
#' for a list of undirected effects \code{\link[remstats:tie_effects]{remstats::tie_effects(directed = FALSE)}}
#' @examples
#' 
#' # Define effects for the start and end model of DuREM
#' start_effects <- ~ 1 + remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
#' end_effects <- ~ 1 + remstats::outdegreeSender(scaling = "std")
#'
#' # Fit a DuREM model
#' durem::duremstimate.grid(start_effects, end_effects, dat$edgelist)
#' 
#' @export
duremstimate.grid <- function(
    start_effects,
    end_effects,
    edgelist,    
    psi_start_candidates = 1,
    psi_end_candidates = 1,
    start_undirected = FALSE,
    end_undirected = FALSE,
    memory = c("full","decay"),
    half_life_candidates = NA,    
    engaged_stat = FALSE,
    engaged_directed = FALSE,
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
                     fit <- duremstimate(start_effects, end_effects, edgelist,
                        psi_start = psi_start_candidates[j],
                        psi_end = psi_end_candidates[k],
                        memory = "decay",
                        half_life = half_life_candidates[l],
                        end_undirected = end_undirected,
                        start_undirected = start_undirected,
                        engaged_stat = engaged_stat,
                        engaged_directed = engaged_directed,
                        strip_return = strip_return)

                        loglik[j,k,l] = as.numeric(logLik(fit))
                        if(!is.null(save_dir)){
                            file_path = paste0(save_dir,"durem_fit_psi_start=",psi_start_candidates[j],"_psi_end=",psi_end_candidates[k],"_half_life=",half_life_candidates[l],".rdata")

                            save(fit,file = file_path)
                        }
                }
            }else{
                fit <- duremstimate(start_effects, end_effects, edgelist,
                        psi_start = psi_start_candidates[j],
                        psi_end = psi_end_candidates[k],
                        memory = "full",
                        end_undirected = end_undirected,
                        start_undirected = start_undirected,
                        engaged_stat = engaged_stat,
                        engaged_directed = engaged_directed,
                        strip_return = strip_return)
                loglik[j,k] = as.numeric(logLik(fit))

                if(!is.null(save_dir)){
                   file_path = paste0(save_dir,"durem_fit_psi_start=",psi_start_candidates[j],"_psi_end=",psi_end_candidates[k],".rdata")

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


#' Estimate Duration Relational Event Model (DuREM)
#'
#' @description 
#' This function estimates the Duration Relational Event Model (DuREM).
#'
#' @param start_effects Formula object for `remstats`, used to compute the start statistics.
#' @param end_effects Formula object for `remstats`, used to compute the end statistics.
#' @param edgelist A `data.frame` with columns: `start_time`, `sender`, `receiver`, and `end_time`.
#' @param psi_start Numeric value of the psi parameter for the start DuREM model. Default is `1`.
#' @param psi_end Numeric value of the psi parameter for the end DuREM model. Default is `1`.
#' @param start_undirected Logical. If `TRUE`, the risk set for start DuREM model is undirected. Default is `FALSE`.
#' @param end_undirected Logical. If `TRUE`, the risk set for the end DuREM model is undirected (see Details). Default is `FALSE`.
#' @param memory Character string indicating memory effects. If `"full"`, no memory effects are incorporated. 
#' If `"decay"`, decay memory effects are applied with the specified half-life. Default is `"full"`.
#' @param half_life Numeric value of the half-life parameter for decay memory. Required if `memory = "decay"`. Default is `NA`.
#' @param engaged_stat Logical. If `TRUE`, includes statistics to account for engagement effects. Default is `FALSE`.
#' @param engaged_directed Logical. If `TRUE`, assumes engagement effects are directed. Default is `TRUE`.
#' @param strip_return Logical. If `TRUE`, strips the heavy elements from a `glm` output object. Default is `TRUE`.
#' 
#' @return 
#' A fitted model object.
#'
#' @details 
#' - The `end_undirected` parameter should be set to `TRUE` if the risk set for the duration model needs to be undirected. 
#'   This means that if a dyad A → B is part of an event, the undirected dyad (AB = BA) is considered at risk to end the event. 
#'   This is useful when it is not directly observable whether the sender or receiver ended the event.
#'
#' - A list of available effects for the start and end models of DuREM can be obtained with 
#'   \code{\link[remstats:tie_effects]{remstats::tie_effects()}}. For a list of undirected effects, use 
#'   \code{\link[remstats:tie_effects]{remstats::tie_effects(directed = FALSE)}}.
#'
#' - The observation period is assumed to start at t = 0. However, the start times and end times must not be 0. 
#'   If necessary, add small random noise to the first start time to avoid this issue.
#'
#' @examples 
#' # Define effects for the start and end models of DuREM
#' start_effects <- ~ 1 + remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
#' end_effects <- ~ 1 + remstats::outdegreeSender(scaling = "std")
#'
#' # Create an edgelist data.frame (example format)
#' edgelist <- data.frame(
#'   start_time = c(1, 2, 3),
#'   sender = c("A", "B", "C"),
#'   receiver = c("B", "C", "A"),
#'   end_time = c(4, 5, 6)
#' )
#'
#' # Fit a DuREM model
#' model <- duremstimate(start_effects, end_effects, edgelist)
#'
#' @export
duremstimate <- function(
    start_effects,
    end_effects,
    edgelist,    
    psi_start = 1,
    psi_end = 1,
    start_undirected = FALSE,
    end_undirected = FALSE,    
    memory = c("full","decay"),
    half_life = NA,    
    engaged_stat = FALSE,
    engaged_directed = TRUE,
    strip_return = TRUE){

    colnames(edgelist)[1:4] = c("start_time","sender","receiver","end_time")

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
    
    stats_list = duremstats(start_effects, end_effects, edgelist,
                    psi_start = psi_start,
                    psi_end = psi_end,
                    memory = memory,
                    half_life = half_life,
                    end_undirected = end_undirected,
                    start_undirected = start_undirected)
    
    stat_stack <- durem_statstack(edgelist = edgelist,
                    start_stats = stats_list$start_stats,
                    end_stats = stats_list$end_stats,
                    end_undirected = end_undirected,
                    start_undirected = start_undirected,
                    engaged_stat = engaged_stat,
                    engaged_directed = engaged_directed)
                    
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
# @param start_effects formula object for remstats, used to compute the start statistics
# @param end_effects formula object for remstats, used to compute the end statistics
# @param edgelist data.frame with columns in order: (start_time, sender, receiver, end_time)
# @param psi_start numeric value of psi parameter for start DuREM model
# @param psi_end numeric value of psi parameter for end DuREM model
# @param memory if "full" then no memory effects are incorporated. If "decay" then decay memory effects with specified half life
# @param half_life numeric value of half life parameter for decay memory
# @param end_undirected TRUE if riskset for the end DuREM model needs to be undirected. See details
# @param start_undirected TRUE if riskset for both start and end DuREM models needs to be undirected
# 
# @details
# end_undirected is set to TRUE if riskset for the duration model needs to be undirected. i.e if dyad A->B is in an event, the undirected dyad (AB==BA) is at risk to end the event. This argument can be used when it is not directly observable whether the sender or receiver ended the event
#  
duremstats <- function(start_effects, 
    end_effects, 
    edgelist, 
    psi_start = 0,
    psi_end = 0,
    memory = c("full","decay"),
    half_life = NA,
    end_undirected = FALSE,
    start_undirected = FALSE){
    
    #colnames(edgelist)[1:4] = c("start_time","sender","receiver","end_time")
    
    #### start model
    edgelist$weight <- (edgelist$end_time - edgelist$start_time+1)^psi_start
    
    #split the edgelist into start and end event types
    dur.edgelist = edgelist[rep(seq_len(nrow(edgelist)), each = 2), ]
    dur.edgelist$type = rep(c("start","end"),nrow(edgelist))
    dur.edgelist$time = ifelse(dur.edgelist$type=="start",dur.edgelist$start_time,dur.edgelist$end_time)
    
    #memory
    memory <- match.arg(memory)
    #memory <- ifelse(is.na(half_life), "full", "decay")
    
    #remify
    suppressWarnings({
        reh_dir <- remify::remify(dur.edgelist[,c("time","sender","receiver","weight","type")], directed = !start_undirected, types = c("start","end"),model="tie")
    })
    
    start_stats <- remstats::tomstats(effects = start_effects, reh = reh_dir, memory= memory, memory_value = half_life)
    dimnames(start_stats)[[3]] = paste0(dimnames(start_stats)[[3]],".start")
    
    #### end model
    edgelist$weight <- (edgelist$end_time - edgelist$start_time+1)^psi_end
    dur.edgelist = edgelist[rep(seq_len(nrow(edgelist)), each = 2), ]
    dur.edgelist$type = rep(c("start","end"),nrow(edgelist))    
    dur.edgelist$time = ifelse(dur.edgelist$type=="start",dur.edgelist$start_time,dur.edgelist$end_time)

    #end model can be either directed or undirected
    suppressWarnings({         
        reh_end <- remify::remify(dur.edgelist[,c("time","sender","receiver","weight","type")], types = c("start","end"), directed = !end_undirected,model="tie")
    })    
 
    end_stats <- remstats::tomstats(effects = end_effects, reh = reh_end, memory = memory, memory_value = half_life)
    
    dimnames(end_stats)[[3]] = paste0(dimnames(end_stats)[[3]],".end")

    return(list(start_stats = start_stats, end_stats = end_stats))
}

#'@keywords internal
create_glm_formula <- function(stat_names){
    return(paste0(" obs ~ -1 + offset(logtimediff) + ",paste(stat_names,collapse = " + ")))
}


#' function to arrange remstats array into a stat stack for glm for the duration rem model
#' 
#' @param edgelist data.frame with columns in order (start_time, sender, receiver, end_time)
#' @param num_actors number of actors
#' @param start_stats remstats object for start stats
#' @param end_stats remstats object for end stats
#' @param end_undirected (default=FALSE) boolean for if the end model should be undirected
#' 
#' @details
#' start_stats has dimension (M,2D,P1) where M is the number of unique time points (over start and end of an event combined)
#' D = total number of dyads in the network = N*(N-1)
#' P1 = number of start statistics specified
#' 
#' end_stats has dimension (M,2D,P2) where P2 = number of end statistics specified
#' 
#' This function allows multiple start and end events to occur simultaneously, but only once per unique dyad in an interval
#' @keywords internal
#' @export
durem_statstack <- function(edgelist,
                    start_stats, end_stats,
                    end_undirected = FALSE,
                    start_undirected = FALSE,
                    engaged_stat = FALSE,
                    engaged_directed = TRUE){
    
    num_actors = length(unique(c(edgelist[,2],edgelist[,3])))

    #colnames(edgelist)[1:4] = c("start_time","sender","receiver","end_time")

    rs = attr(start_stats,"riskset")
    #start and end event type dyad id
    if(!start_undirected){
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

    if(end_undirected & !start_undirected){
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
    
    P_durem = 0
    if(engaged_stat){
        P_durem = P_durem + 2 #engaged actor for start and end model
        #if(engaged_directed){
        if(engaged_directed){ #one each for sender and recv
            P_durem = P_durem + 2
        }
    }
    
    stat_stack <- do.call(rbind.data.frame, lapply(1:M, function(i){
        
        stack_row = data.frame(matrix(0,nrow = 0,ncol = (5 + dim(start_stats)[3]+ dim(end_stats)[3] + P_durem)))
                  
        engaged_actors <- tabulate(as.numeric(unlist(edgelist[edgelist$start_time <= unique_times[i] & edgelist$end_time > unique_times[i], 2:3], use.names = FALSE)), nbins = num_actors)



        ###### type 1::State: an event ended
        # end(d) is observed to end
        #Also allows instantaneous events because start_time <= t
        if(i > 1){
            observed_events = subset(edgelist,start_time <= unique_times[i] & end_time == unique_times[i])
            if(end_undirected & !start_undirected){
                dyads_observed = unique(observed_events[,"end_dyad_und"])
            }else{
                dyads_observed = unique(observed_events[,"end_dyad"])    
            }
            if(length(dyads_observed) > 0){
                stack_row = rbind(stack_row,as.data.frame(do.call(rbind,lapply(dyads_observed, function(d){
                    if(P_durem == 0){
                        return(c(1, d, i, logtimediff[i],1, rep(0,P_start),unname(end_stats[i,d,])))
                    }else{
                        # if sender or receiver of d is currently involved in an event then stat_durem = c(0,n) otherwise c(0,0)
                        if(end_undirected & !start_undirected){
                            sender = as.numeric(dur.rs[d,1])
                            receiver = as.numeric(dur.rs[d,2])
                        }else{
                            sender = as.numeric(rs[d,1])
                            receiver = as.numeric(rs[d,2])
                        }
                        if(engaged_stat){
                            # Count engaged start overlaps
                            if(engaged_directed){
                                engaged_send = engaged_actors[sender]
                                engaged_recv = engaged_actors[receiver]
                                stat_durem = c(0,0, engaged_send, engaged_recv)
                            }else{
                                engaged_end = engaged_actors[sender] + engaged_actors[receiver] - 2                                
                            stat_durem = c(0, engaged_end)
                            }
                        }
                        
                        
                        return(c(1, d, i, logtimediff[i],1, rep(0,P_start),unname(end_stats[i,d,]),stat_durem))
                    }
                    
                }))))
                
            }    
        }

        ###### type 2::State: event ongoing
        # end(d) is at risk to end
        #all events active right now but didnt start in the interval
        active_events = subset(edgelist,start_time < unique_times[i] & end_time > unique_times[i])
        #end dyads currently in an event 
        if(end_undirected& !start_undirected){
            dyads_at_risk_end = unique(active_events[,"end_dyad_und"])
        }else{
            dyads_at_risk_end = unique(active_events[,"end_dyad"])
        }
        
        #add the active events with end type to riskset
        if(length(dyads_at_risk_end>0)){
            stack_row = rbind(stack_row,as.data.frame(do.call(rbind,lapply(dyads_at_risk_end, function(d){
                if(P_durem == 0){
                    return(c(0, d, i, logtimediff[i] ,2,rep(0,P_start),unname(end_stats[i,d,])))
                }else{                    
                        # if sender or receiver of d is currently involved in an event then stat_durem = c(n,0) otherwise c(0,0)
                        if(end_undirected & !start_undirected){
                            sender = as.numeric(dur.rs[d,1])
                            receiver = as.numeric(dur.rs[d,2])
                        }else{
                            sender = as.numeric(rs[d,1])
                            receiver = as.numeric(rs[d,2])
                        }
                        if(engaged_stat){
                            # Count engaged start overlaps
                            if(engaged_directed){
                                engaged_send = engaged_actors[sender] - 1 
                                engaged_recv = engaged_actors[receiver] -1

                                stat_durem = c(0,0, engaged_send, engaged_recv)
                            }else{
                                engaged_end = engaged_actors[sender] + engaged_actors[receiver] - 2                                
                                stat_durem = c(0, engaged_end)
                            }
                        }                    
                    return(c(0, d, i, logtimediff[i] ,2,rep(0,P_start),unname(end_stats[i,d,]),stat_durem))
                }
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
                if(P_durem==0){
                    return(c(1, d, i, logtimediff[i] ,3, unname(start_stats[i,d,]),rep(0,P_end)))
                }else{                    
                        # if sender or receiver of d is currently involved in an event then stat_durem = c(n,0) otherwise c(0,0)
                        if(end_undirected & !start_undirected){
                            sender = as.numeric(dur.rs[d,1])
                            receiver = as.numeric(dur.rs[d,2])
                        }else{
                            sender = as.numeric(rs[d,1])
                            receiver = as.numeric(rs[d,2])
                        }
                        if(engaged_stat){
                            # Count engaged start overlaps
                            if(engaged_directed){
                                engaged_send = engaged_actors[sender] - 1 
                                engaged_recv = engaged_actors[receiver] -1
                                stat_durem = c(engaged_send, engaged_recv, 0, 0)
                            }else{
                                engaged_start = engaged_actors[sender] + engaged_actors[receiver] - 2                               
                            stat_durem = c(engaged_start, 0)
                            }
                        }                  
                    return(c(1, d, i, logtimediff[i] ,3, unname(start_stats[i,d,]),rep(0,P_end),stat_durem))
                }            
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
                if(P_durem==0){
                    return(c(0, d, i, logtimediff[i] ,4,unname(start_stats[i,d,]),rep(0,P_end)))
                }else{
                        # if sender or receiver of d is currently involved in an event then stat_durem = c(1,0) otherwise c(0,0)
                        # Count engaged start overlaps
                        if(end_undirected & !start_undirected){
                            sender = as.numeric(dur.rs[d,1])
                            receiver = as.numeric(dur.rs[d,2])
                        }else{
                            sender = as.numeric(rs[d,1])
                            receiver = as.numeric(rs[d,2])
                        }
                        if(engaged_stat){
                            # Count engaged start overlaps
                            if(engaged_directed){
                                engaged_send = engaged_actors[sender]
                                engaged_recv = engaged_actors[receiver]
                                stat_durem = c(engaged_send, engaged_recv, 0, 0)
                            }else{
                                engaged_start = engaged_actors[sender] + engaged_actors[receiver]                               
                                stat_durem = c(engaged_start, 0)
                            }
                        }                                    
                    return(c(0, d, i, logtimediff[i] ,4,unname(start_stats[i,d,]),rep(0,P_end),stat_durem))
                }
            }))))            
        }
        return(stack_row)
    }))
    
    if(P_durem==0){
        colnames(stat_stack) = c(c("obs","tie","i","logtimediff","type"),dimnames(start_stats)[[3]],dimnames(end_stats)[[3]])
    }else if(engaged_stat){
        if(engaged_directed){
            colnames(stat_stack) = c(c("obs","tie","i","logtimediff","type"),dimnames(start_stats)[[3]],dimnames(end_stats)[[3]],"engaged_sender.start","engaged_receiver.start","engaged_sender.end","engaged_receiver.end")
        }else{
            colnames(stat_stack) = c(c("obs","tie","i","logtimediff","type"),dimnames(start_stats)[[3]],dimnames(end_stats)[[3]],"engaged.start","engaged.end")
        }
        
    }
    
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

#' Estimate hyperparameters for Duration Relational Event Model (DuREM)
#' 
#' @description Function to estimate the hyperparameters for duration relational event model (DuREM) using grid search
#' 
#'  
#' @param start_effects formula object for remstats, used to compute the start statistics
#' @param end_effects formula object for remstats, used to compute the end statistics
#' @param edgelist data.frame with columns (start_time, sender, receiver, end_time)
#' @param psi_start_candidates numeric vector of candidate values for psi_start. default value is 1
#' @param psi_end_candidates numeric vector of candidate values for psi_end. default value is 1
#' @param memory if "full" then no memory effects are incorporated. If "decay" then decay memory effects with specified half life
#' @param half_life_candidates numeric vector of candidate values for half life parameters
#' @param end_undirected TRUE if riskset for the end DuREM model needs to be undirected. See details
#' @param start_undirected TRUE if riskset for both start and end DuREM models needs to be undirected
#' @param strip_return logical, if TRUE then strip the heavy elements from a glm output object
#' @param save_dir character, local directory where to save fitted candidate model files
#' 
#' @return list with element \code{loglik}, a matrix or array of loglikelihood of fitted DuREM candidate models and \code{mle}, a vector of candidate values (psi_start,psi_end(,half_life)) with maximum likelihood
#' @details
#' \code{end_undirected} is set to \code{TRUE} if riskset for the duration model needs to be undirected. i.e if dyad A->B is in an event, the undirected dyad (AB==BA) is at risk to end the event. This argument can be used when it is not directly observable whether the sender or receiver ended the event
#' 
#' A list of available effects for the start and end models of DuREM can be obtained with \code{\link[remstats:tie_effects]{remstats::tie_effects()}} and
#' for a list of undirected effects \code{\link[remstats:tie_effects]{remstats::tie_effects(directed = FALSE)}}
#' @examples
#' 
#' # Define effects for the start and end model of DuREM
#' start_effects <- ~ 1 + remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
#' end_effects <- ~ 1 + remstats::outdegreeSender(scaling = "std")
#'
#' # Fit a DuREM model
#' durem::duremstimate.grid(start_effects, end_effects, dat$edgelist)
#' 
#' @export
duremstimate.grid.parallel <- function(cl,
                                      start_effects,
                                      end_effects,
                                      edgelist,    
                                      psi_start_candidates = 1,
                                      psi_end_candidates = 1,
                                      memory = c("full", "decay"),
                                      half_life_candidates = NA,
                                      end_undirected = FALSE,
                                      start_undirected = FALSE,
                                      strip_return = TRUE,
                                      save_dir = NULL) {
  
  if (any(edgelist$end_time < edgelist$start_time, na.rm = TRUE)) {
      stop("End time of an event cannot be before its start time.")
  }
  if (anyNA(edgelist$end_time)) {
      stop("Missing event end time")
  }
  
  memory <- match.arg(memory)
 
  if (memory == "decay" && any(is.na(half_life_candidates))) {
      stop("Incorrect half life supplied for decay memory.")
  }
  
  # Prepare the parameter grid
  param_grid <- expand.grid(psi_start = psi_start_candidates,
                            psi_end = psi_end_candidates,
                            half_life = ifelse(memory == "decay", half_life_candidates, NA),
                            stringsAsFactors = FALSE)
  
  # Check and convert param_grid to ensure proper data type handling
  if (!is.data.frame(param_grid)) {
      param_grid <- as.data.frame(param_grid)
  }

  # Execute the function in parallel
  results <- parLapply(cl, param_grid, function(params) {
    if(memory == "decay"){
        fit <- duremstimate(start_effects, end_effects, edgelist,
            psi_start = params[['psi_start']],
            psi_end = params[['psi_end']],
            memory = "decay",
            half_life = params[['half_life']],
            end_undirected = end_undirected,
            start_undirected = start_undirected,
            strip_return = strip_return)
        if(!is.null(save_dir)){
            file_path = paste0(save_dir,"durem_fit_psi_start=",psi_start_candidates[j],"_psi_end=",psi_end_candidates[k],"_half_life=",half_life_candidates[l],".rdata")

            save(fit,file = file_path)
        }                
    }else{
        fit <- duremstimate(start_effects, end_effects, edgelist,
            psi_start = params[['psi_start']],
            psi_end = params[['psi_end']],
            memory = "full",
            end_undirected = end_undirected,
            start_undirected = start_undirected,
            strip_return = strip_return)

        if(!is.null(save_dir)){
            file_path = paste0(save_dir,"durem_fit_psi_start=",psi_start_candidates[j],"_psi_end=",psi_end_candidates[k],".rdata")

            save(fit,file = file_path)
        }
    }    

    loglik <- as.numeric(logLik(fit))

    return(loglik)
  })

  if(memory == "decay") {
    loglik <- array(data = unlist(results), dim = c(length(psi_start_candidates), length(psi_end_candidates), length(half_life_candidates)))
    dimnames(loglik)[[3]] = half_life_candidates
  } else {
    loglik <- matrix(data = unlist(results), nrow = length(psi_start_candidates), ncol = length(psi_end_candidates))
  }

  dimnames(loglik)[[1]] = psi_start_candidates
  dimnames(loglik)[[2]] = psi_end_candidates

  # Determine the maximum likelihood estimates
  indices <- which(loglik == max(loglik), arr.ind = TRUE)
  mle = c(dimnames(loglik)[[1]][indices[1]], dimnames(loglik)[[2]][indices[2]])
  if(memory == "decay") {
    mle = c(mle, dimnames(loglik)[[3]][indices[3]])
  }

  return(list(loglik = loglik, mle = mle))
}
