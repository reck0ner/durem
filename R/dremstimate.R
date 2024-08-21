#' Function to estimate duration relational event model (DREM)
#' 
#'  
#' @param start_effects formula object for remstats, used to compute the start statistics
#' @param end_effects formula object for remstats, used to compute the end statistics
#' @param edgelist data.frame with columns (start_time, sender, receiver, end_time)
#' @param actors Vector, character vector of actor names
#' @param dur_undirected Logical, if TRUE, riskset for the duration model needs to be undirected. See details
#' @param strip_return Logical, if TRUE then strip the heavy elements from a glm output object
#' 
#' @details
#' dur_undirected is set to TRUE if riskset for the duration model needs to be undirected. i.e if dyad A->B is in an event, the undirected dyad (AB==BA) is at risk to end the event. This argument can be used when it is not directly observable whether the sender or receiver ended the event
#' 
#' @examples
#' 
#' # Define effects for the start and end model of DREM
#' start_effects <- ~ 1 + remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
#' end_effects <- ~ 1 + remstats::outdegreeSender(scaling = "std")
#'
#' Fit a DREM model
#' drem::dremstimate(start_effects, end_effects, dat$edgelist)
#' 
#' @export
dremstimate <- function(
    start_effects,
    end_effects,
    edgelist,
    dur_undirected = FALSE,
    strip_return = TRUE){

    stats_list = dur.remstats(start_effects, end_effects, edgelist,dur_undirected = dur_undirected)
    
    stat_stack <- dur.statStack(edgelist = edgelist, stats_list$start_stats, stats_list$end_stats, dur_undirected = dur_undirected)
    
    stat_names = colnames(stat_stack)[6:length(colnames(stat_stack))]
    
    est <- glm(create_glm_formula(stat_names),
               family = poisson(),
               data = stat_stack)
    
    
    #Make the output smaller
    if(strip_return){
        est <- strip_glm(est)
    }
    
    return(est)
}



#' function to compute preliminary stats for the duration rem using remstats 
#' The stats generated from this function still require processing to be used to estimate the model
#' 
#' @param start_effects formula object for remstats, used to compute the start statistics
#' @param end_effects formula object for remstats, used to compute the end statistics
#' @param edgelist data.frame with columns (start_time, sender, receiver, end_time)
#' @param actors character vector of actor names
#' @param memory_value value of half life parameter for exponential decay of memory
#' @param dur_undirected TRUE if riskset for the duration model needs to be undirected. See details
#' 
#' @details
#' dur_undirected is set to TRUE if riskset for the duration model needs to be undirected. i.e if dyad A->B is in an event, the undirected dyad (AB==BA) is at risk to end the event. This argument can be used when it is not directly observable whether the sender or receiver ended the event
#' 
#' @keywords internal
dur.remstats <- function(start_effects, end_effects, edgelist, actors = NULL, psi = 0,memory = "full",memory_value = Inf, dur_undirected = FALSE){
    colnames(edgelist) = c("start_time","sender","receiver","end_time")
    #edgelist$time_end <- edgelist$time + edgelist$duration
    edgelist$weight <- (edgelist$end_time - edgelist$start_time)^psi
    
    dur.edgelist = edgelist[rep(seq_len(nrow(edgelist)), each = 2), ]
    dur.edgelist$type = rep(c("start","end"),nrow(edgelist))
    
    dur.edgelist$time = ifelse(dur.edgelist$type=="start",dur.edgelist$start_time,dur.edgelist$end_time)
    
    #memory
    memory <- ifelse(memory_value == Inf, "full", "decay")
    
    #remify
    suppressWarnings({
        reh_dir <- remify::remify(dur.edgelist[,c("time","sender","receiver","weight","type")], actors=actors, types = c("start","end"))
    })
    
    start_stats <- remstats::tomstats(effects = start_effects, reh = reh_dir, memory, memory_value)
    dimnames(start_stats)[[3]] = paste0(dimnames(start_stats)[[3]],".start")
    
    if(dur_undirected){
        suppressWarnings({
            reh_undir <- remify::remify(dur.edgelist[,c("time","sender","receiver","weight","type")],actors=actors, types = c("start","end"), directed = FALSE)
        })
        end_stats <- remstats::tomstats(effects = end_effects, reh = reh_undir, memory, memory_value)    
    }else{
        end_stats <- remstats::tomstats(effects = end_effects, reh = reh_dir, memory, memory_value)    
    }
    
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
#' @param start_stats remstats object for start stats
#' @param end_stats remstats object for end stats
#' @param dur_undirected (default=FALSE) boolean for if the end model should be undirected
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
dur.statStack <- function(edgelist, start_stats, end_stats, dur_undirected = FALSE){
    
    colnames(edgelist) = c("start_time","sender","receiver","end_time")
    
    rs = attr(start_stats,"riskset")
    #start and end event type dyad id
    edgelist$start_dyad <- apply(edgelist,1,function(x){
        return(which(rs[,1]== x["sender"] & rs[,2]==x["receiver"] & rs[,3] == "start"))
    })
    edgelist$end_dyad <- apply(edgelist,1,function(x){
        return(which(rs[,1]== x["sender"] & rs[,2]==x["receiver"] & rs[,3] == "end"))
    })    
    
    if(dur_undirected){
        dur.rs = attr(end_stats,"riskset")
        edgelist$end_dyad_und <- apply(edgelist,1,function(x){
            ind1 = which(dur.rs[,2]== as.numeric(x["sender"]) & dur.rs[,1]==as.numeric(x["receiver"]) & dur.rs[,3] == "end")
            ind2 = which(dur.rs[,1]== as.numeric(x["sender"]) & dur.rs[,2]==as.numeric(x["receiver"]) & dur.rs[,3] == "end")
            return(max(ind1,ind2))
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
            if(dur_undirected){
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
        if(dur_undirected){
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