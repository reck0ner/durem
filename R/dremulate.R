#' Simulate Temporal Events with Duration - Tie based duration model
#' @description 
#'  A function to simulate relational event data by sampling from a
#' tie based duration relational event model.
#'
#' @details . 
#'
#' @param start_effects an object of type \code{formula} for specification of statistics used to simulate the network. 
#' @param stop_effects an object of type \code{formula} for specification of statistics used to simulate the network. 
#' @param num_actors Integer, number of actors in the network.
#' @param num_events [Optional] Integer, maximum number of events to simulate.
#' @return \describe{
#' \item{edgelist}{data.frame object with columns (start_time,sender,receiver,end_time)}
#' \item{evls}{matrix containing the event list  with columns (dyad,start_time,end_time) where dyad represents the index of the dyad or the (sender,receiver) pair in the riskset}
#' \item{riskset}{pair of sender,receiver ids corresponding to dyad index in the riskset}
#' }
#' @examples 
#'  
#' @export
dremulateTie <- function(
  start_effects,
  stop_effects,
  start_params,
  stop_params,
  num_actors,
  num_events,
  dur_undirected = FALSE,
  stop_threshold = NULL
  ){
	if(is.null(stop_threshold)){
		stop_threshold = 3 * num_events
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
	#while(t <= time){
	while((end_count < num_events & start_count < stop_threshold) | (start_count < num_events & end_count < stop_threshold)){
		#updating event rate / lambda
		lambda <- sapply(1:nrow(rs),function(x){
			if(dyads_is.active[x]){
				if(end_count == 0){
					return(exp(stop_params[1])) #just baseline effect for first event
				}				
				if(dim(dur_stats)[3]==1){
					return(exp(dur_stats[1,x,] * stop_params))
				}
				return(exp(dur_stats[1,x,] %*% stop_params)) 
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
#print(2)
		#sampling waiting time dt
		if (waiting_time == "exp") {
			dt <- rexp(1, rate = sum(lambda))
			t = t + dt
		}
#print(3)
		#sampling dyad for next event
		dyad <- sample(1:nrow(rs), 1, prob = lambda / sum(lambda))
#print(4)
		if(dyads_is.active[dyad]){
			# cat("stop time: ",t,"\n")
			#stop event
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
			dur_stats <- remstats::tomstats(effects = stop_effects, reh = dur_reh,start = length(end_indx), stop = length(end_indx))
			
		}else{
			#  cat("start time: ",t,"\n")
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
		# cat("start count: ",start_count," end count: ",end_count)
		#stop if max number of events reached
		if(start_count>= stop_threshold){
			cat(paste0("\n Stopping: specified model generates insufficient stop events \n number of stop events sampled = ",end_count," \n"))
			return(list(
				edgelist = edgelist,
				evls = evls,
				start_params = start_params,
				stop_params = stop_params
    			))
		}
		if(end_count>= stop_threshold){
			cat(paste0("\n Stopping: specified model generates insufficient start events \n number of start events sampled = ",start_count," \n"))
			return(list(
				edgelist = edgelist,
				evls = evls,
				start_params = start_params,
				stop_params = stop_params
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
		stop_params = stop_params,
		riskset = rs
    )
  )
}
