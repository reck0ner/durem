#' Simulate Temporal Events with Duration Relational Event Model (DREM)
#' @description 
#'  A function to simulate relational event data by sampling from a
#' tie based duration relational event model.
#'
#'
#' @param start_effects an object of class \code{"\link[stats]{formula}"}:  a symbolic description of the effects to
#' simulate the start model of DREM. see 'Details' for
#' the available effects and their corresponding statistics 
#' @param end_effects an object of class \code{"\link[stats]{formula}"}:  a symbolic description of the effects to
#' simulate the end model of DREM. see 'Details' for
#' the available effects and their corresponding statistics 
#' @param num_actors Integer, number of actors in the network.
#' @param num_events Integer, maximum number of events to simulate.
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
  dur_undirected = FALSE,
  event_threshold = NULL
  ){
	if(is.null(event_threshold)){
		event_threshold = 3 * num_events
	}
	waiting_time = "exp"

	dummy <- data.frame(time = c(1:2), actor1 = c(1,2), actor2 = c(2,3))

	rate_reh <- remify::remify(dummy, actors = 1:num_actors,model="tie",directed = TRUE)

	rate_stats <- remstats::tomstats(effects = ~ 1 , reh = rate_reh)

	rs <- attr(rate_stats, 'riskset')

	if(dur_undirected){
		dur_reh <- remify::remify(dummy, actors = 1:num_actors,model="tie",directed = FALSE)
		dur_stats <- remstats::tomstats(effects = ~ 1 , reh = dur_reh)
		undir_rs <- attr(dur_stats, 'riskset')
	}else{
		dur_stats <- remstats::tomstats(effects = ~ 1 , reh = rate_reh)		
	}
	edgelist <- data.frame() #start time, sender, receiver, end_time
	evls <- data.frame()

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
				if(dim(dur_stats)[3]==1){
					return(exp(dur_stats[1,x,] * end_params))
				}
				return(exp(dur_stats[1,x,] %*% end_params)) 
			}else{
				if(start_count == 0){
					return(exp(start_params[1])) #just baseline effect for first event
				}
				if(dim(rate_stats)[3]==1){
					return(exp(rate_stats[1,x,] * start_params))
				}
				return(exp(rate_stats[1,x,] %*% start_params))
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
			evls[x,3] = t
			end_count = end_count + 1
			dyads_is.active[dyad] = FALSE				
			end_indx = which(edgelist[,4]>0)
			suppressWarnings({
			dur_reh <- remify::remify(edgelist[end_indx,c(4,2,3)],actors = 1:num_actors,model="tie")
			})
			dur_stats <- remstats::tomstats(effects = end_effects, reh = dur_reh, start = length(end_indx), stop = length(end_indx))
			
		}else{

			#start an event
			if(start_count < num_events){
				#update dyad/sender/receiver
				evls = rbind(evls,data.frame(dyad = dyad, time = t, end_time = NA))
				edgelist = rbind(edgelist,data.frame(time = t, actor1 = rs[dyad,1], actor2= rs[dyad,2], end_time = NA))

				dyads_is.active[dyad] = TRUE				
			}     
			
			rate_reh <- remify::remify(edgelist[,c(1,2,3)],actors = 1:num_actors,model="tie")			
			
			rate_stats <- remstats::tomstats(effects = start_effects, reh = rate_reh, start = nrow(edgelist), stop = nrow(edgelist))
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
		riskset = rs
    )
  )
}
