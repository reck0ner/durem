#' Simulate Temporal Events with Duration Relational Event Model (DREM)
#' @description 
#'  A function to simulate relational event data by sampling from a
#' tie based duration relational event model.
#'
#'
#' @param start_effects formula object:  a symbolic description of the effects to
#' simulate the start model of DREM. see 'Details' for
#' the available effects and their corresponding statistics 
#' @param end_effects formula object:  a symbolic description of the effects to
#' simulate the end model of DREM. see 'Details' for
#' the available effects and their corresponding statistics 
#' @param num_actors Integer, number of actors in the network.
#' @param num_events Integer, maximum number of events to simulate.
#' @param psi_start Numeric, value of psi parameter for start rate model
#' @param psi_end Numeric, value of psi parameter for end rate model
#' @param dur_undirected [Optional] Logical, if TRUE, riskset for the duration model needs to be undirected. See details
#' @param event_threshold [Optional] Integer, maximum number of  incomplete start or end 
#' events to simulate before stopping the simulation
#'
#' @return \describe{
#' \item{edgelist}{data.frame object with columns (start_time,sender,receiver,end_time)}
#' \item{evls}{matrix containing the event list  with columns (dyad,start_time,end_time) where dyad represents the index of the dyad or the (sender,receiver) pair in the riskset}
#' \item{riskset}{pair of sender,receiver ids corresponding to dyad index in the riskset}
#' }
#' @examples 
#' # Define effects for the start and end model of DREM
#' start_effects <- ~ 1 + remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
#' end_effects <- ~ 1 + remstats::outdegreeSender(scaling = "std")
#'
#' # Set model parameters
#' start_params <- c(-7, 0.2, 0.1)
#' end_params <- c(-4, -0.2)
#'
#' # Run the simulation with specified parameters
#' drem::dremulate(start_effects, end_effects, start_params, end_params, num_actors = 10, num_events = 1000, event_threshold = 1500)
#' @details 
#' A list of available effects for the start and end models of DREM can be obtained with \code{\link[remstats:tie_effects]{remstats::tie_effects()}} and
#' for a list of undirected effects \code{\link[remstats:tie_effects]{remstats::tie_effects(directed = FALSE)}}
#' 
#' @export
dremulate <- function(
  start_effects,
  end_effects,
  start_params,
  end_params,
  num_actors,
  num_events,
  psi_start = 1,
  psi_end = 1,
  dur_undirected = FALSE,
  event_threshold = NULL
  ){
	if(is.null(event_threshold)){
		event_threshold = 3 * num_events
	}
	waiting_time = "exp"

	dummy <- data.frame(time = c(1:2), actor1 = c(1,2), actor2 = c(2,3))

	start_reh <- remify::remify(dummy, actors = 1:num_actors,model="tie",directed = TRUE)

	start_stats <- remstats::tomstats(effects = ~ 1 , reh = start_reh)

	rs <- attr(start_stats, 'riskset')

	if(dur_undirected){
		end_reh <- remify::remify(dummy, actors = 1:num_actors,model="tie",directed = FALSE)
		end_stats <- remstats::tomstats(effects = ~ 1 , reh = end_reh)
		undir_rs <- attr(end_stats, 'riskset')
	}else{
		end_stats <- remstats::tomstats(effects = ~ 1 , reh = start_reh)		
	}
	edgelist <- data.frame() #start time, sender, receiver, end_time
	evls <- data.frame()

	start_weighted_edgelist = data.frame() #start_time,sender,recv,start_weight
	end_weighted_edgelist = data.frame() # end_time,sender,recv,end_weight

	dyads_is.active = rep(FALSE,nrow(rs))

	#counters for sampled events
	end_count = 0
	start_count = 0
	t = 0
	while((end_count < num_events & start_count < event_threshold) | (start_count < num_events & end_count < event_threshold)){
		#updating event rate / lambda
		lambda <- sapply(1:nrow(rs),function(x){
			if(dyads_is.active[x]){
				if(end_count == 0){
					return(exp(end_params[1])) #just baseline effect for first event
				}				
				if(dim(end_stats)[3]==1){
					return(exp(end_stats[1,x,] * end_params))
				}
				return(exp(end_stats[1,x,] %*% end_params)) 
			}else{
				if(start_count == 0){
					return(exp(start_params[1])) #just baseline effect for first event
				}
				if(dim(start_stats)[3]==1){
					return(exp(start_stats[1,x,] * start_params))
				}
				return(exp(start_stats[1,x,] %*% start_params))
			}
		})	

		#sampling waiting time dt
		if (waiting_time == "exp") {
			dt <- rexp(1, rate = sum(lambda))
			t = t + dt
		}

		#sampling dyad for next event
		dyad <- sample(1:nrow(rs), 1, prob = lambda / sum(lambda))

		if(dyads_is.active[dyad]){

			#end event
			x = max(which(evls[,1]==dyad))
			#update end times
			edgelist[x,4] = t
			#update event weight
			end_weighted_edgelist[x,1] = t
			start_weighted_edgelist[x,4] = (edgelist[x,4] - edgelist[x,1])^psi_start
			end_weighted_edgelist[x,4] = (edgelist[x,4] - edgelist[x,1])^psi_end
			evls[x,3] = t
			end_count = end_count + 1
			dyads_is.active[dyad] = FALSE				
			end_indx = which(edgelist[,4]>0)
			suppressWarnings({
			#weighted edgelist
			#end_reh <- remify::remify(edgelist[end_indx,c(4,2,3,6)],actors = 1:num_actors,model="tie")
			end_reh <- remify::remify(end_weighted_edgelist, actors = 1:num_actors,model="tie")
			})
			end_stats <- remstats::tomstats(effects = end_effects, reh = end_reh, start = length(end_indx), stop = length(end_indx))
			
		}else{

			#start an event
			if(start_count < num_events){
				#update dyad/sender/receiver
				evls = rbind(evls,data.frame(dyad = dyad, time = t, end_time = NA))
				edgelist = rbind(edgelist,data.frame(time = t, actor1 = rs[dyad,1], actor2 = rs[dyad,2], end_time = NA))
				
				start_weighted_edgelist = rbind(start_weighted_edgelist,data.frame(time = t, actor1 = rs[dyad,1], actor2 = rs[dyad,2], weight = 1))
				end_weighted_edgelist = rbind(end_weighted_edgelist,data.frame(time = NA, actor1 = rs[dyad,1], actor2 = rs[dyad,2], weight = 0))
				dyads_is.active[dyad] = TRUE				
			}     
			#weighted edgelist
			#start_reh <- remify::remify(edgelist[,c(1,2,3,5), actors = 1:num_actors,model="tie")			
			start_reh <- remify::remify(start_weighted_edgelist, actors = 1:num_actors,model="tie")			

			start_stats <- remstats::tomstats(effects = start_effects, reh = start_reh, start = nrow(edgelist), stop = nrow(edgelist))
			start_count = start_count + 1
		}
		#end if max number of events reached
		if(start_count>= event_threshold){
			cat(paste0("\n endping: specified model generates insufficient end events \n number of end events sampled = ",end_count," \n"))
			return(list(
				edgelist = edgelist,
				evls = evls,
				start_params = start_params,
				end_params = end_params
    			))
		}
		if(end_count>= event_threshold){
			cat(paste0("\n endping: specified model generates insufficient start events \n number of start events sampled = ",start_count," \n"))
			return(list(
				edgelist = edgelist,
				evls = evls,
				start_params = start_params,
				end_params = end_params
			))
		}
		# Calculate progress percentage
  		progress <- min(99,end_count / num_events * 100)

  		# Update the progress bar in the console
  		cat(sprintf("\rProgress: [%-20s] %.1f%%", paste0(rep("=", progress/5), collapse = ""), progress))
  		
	}

	cat(sprintf("\rProgress: [%-20s] %.1f%%", paste0(rep("=", 20), collapse = ""), 100))
	flush.console()
  return(
    list(
        edgelist = edgelist,
		evls = evls,
		start_params = start_params,
		end_params = end_params,
		psi_start = psi_start,
		psi_end = psi_end,
		riskset = rs
    )
  )
}
