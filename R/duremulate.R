#' Simulate Temporal Events with Duration Relational Event Model (DuREM)
#'
#' @description 
#' Simulates relational event data by sampling from a tie-based Duration Relational Event Model (DuREM).
#'
#' @param start_effects Formula object. A symbolic description of the effects used to simulate the start model of DuREM. 
#' See 'Details' for the available effects and their corresponding statistics.
#' @param end_effects Formula object. A symbolic description of the effects used to simulate the end model of DuREM. 
#' See 'Details' for the available effects and their corresponding statistics.
#' @param start_params Numeric vector. Parameters corresponding to the `start_effects` in the DuREM start model.
#' @param end_params Numeric vector. Parameters corresponding to the `end_effects` in the DuREM end model.
#' @param num_actors Integer. The number of actors in the network.
#' @param num_events Integer. The maximum number of events to simulate.
#' @param psi_start Numeric. The value of the psi parameter for the start rate model. Default is `1`.
#' @param psi_end Numeric. The value of the psi parameter for the end rate model. Default is `1`.
#' @param end_undirected Logical (optional). If `TRUE`, the risk set for the duration model is undirected (see Details). Default is `FALSE`.
#' @param event_threshold Integer (optional). The maximum number of incomplete start or end events to simulate before stopping the simulation. Default is `NULL` (no threshold).
#' @param engaged_stat Logical (optional). If `TRUE`, includes statistics to account for engagement effects. Default is `FALSE`.
#' @param start_engaged_params Numeric vector (optional). Parameters for engagement-related effects in the start model. Default is `NULL`.
#' @param end_engaged_params Numeric vector (optional). Parameters for engagement-related effects in the end model. Default is `NULL`.
#'
#' @return 
#' A list containing:
#' \describe{
#'   \item{edgelist}{A `data.frame` with columns `start_time`, `sender`, `receiver`, and `end_time`, representing the sequence of observed relational events.}
#'   \item{evls}{A matrix containing the event list with columns `dyad`, `start_time`, and `end_time`, 
#'   where `dyad` represents the index of the dyad or the (`sender`, `receiver`) pair in the risk set.}
#'   \item{start_params}{A numeric vector containing the parameters used for the start rate model.}
#'   \item{end_params}{A numeric vector containing the parameters used for the end rate model.}
#'   \item{psi_start}{A numeric value representing the psi parameter used in the start rate model.}
#'   \item{psi_end}{A numeric value representing the psi parameter used in the end rate model.}
#'   \item{riskset}{A matrix with columns `sender` and `receiver`, representing the possible dyads included in the risk set.}
#' }

#' @details 
#' - A list of available effects for the start and end models of DuREM can be obtained with 
#'   \code{\link[remstats:tie_effects]{remstats::tie_effects()}}. For a list of undirected effects, use 
#'   \code{\link[remstats:tie_effects]{remstats::tie_effects(directed = FALSE)}}.
#'
#' - The `end_undirected` parameter should be set to `TRUE` if the risk set for the duration model needs to be undirected. 
#'   For example, if a dyad A → B is part of an event, the undirected dyad (AB = BA) is considered at risk to end the event.
#'
#' - If `engaged_stat` is `TRUE`, additional engagement-related statistics are included in the simulation. 
#'   The parameters for these effects can be specified using `start_engaged_params` for the start model and `end_engaged_params` for the end model.
#'
#' @examples 
#' # Define effects for the start and end models of DuREM
#' start_effects <- ~ 1 + remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
#' end_effects <- ~ 1 + remstats::outdegreeSender(scaling = "std")
#'
#' # Set model parameters
#' start_params <- c(-7, 0.2, 0.1)
#' end_params <- c(-4, -0.2)
#'
#' # Run the simulation with specified parameters
#' simulated_data <- durem::duremulate(
#'   start_effects, 
#'   end_effects, 
#'   start_params, 
#'   end_params, 
#'   num_actors = 10, 
#'   num_events = 1000, 
#'   event_threshold = 1500, 
#'   engaged_stat = TRUE, 
#'   start_engaged_params = c(0.2, 0.2), 
#'   end_engaged_params = c(0.2, 0.2)
#' )
#'
#' @export
duremulate <- function(
  start_effects,
  end_effects,
  start_params,
  end_params,
  num_actors,
  num_events,
  psi_start = 1,
  psi_end = 1,
  start_undirected = FALSE,
  end_undirected = FALSE,
  event_threshold = NULL,
  engaged_stat = FALSE,
  start_engaged_params = NULL,
  end_engaged_params = NULL
  ){
	if(is.null(event_threshold)){
		event_threshold = 3 * num_events
	}
	waiting_time = "exp"

	dummy <- data.frame(time = c(1:2), actor1 = c(1,2), actor2 = c(2,3))

	end_undirected = ifelse(start_undirected,TRUE,FALSE) #if reh is undirected dur mst be undirected

	if(!start_undirected){
		start_reh <- remify::remify(dummy, actors = 1:num_actors,model="tie",directed = TRUE)
	}else{
		start_reh <- remify::remify(dummy, actors = 1:num_actors,model="tie",directed = FALSE)
	}

	start_stats <- remstats::tomstats(effects = ~ 1 , reh = start_reh)

	rs <- attr(start_stats, 'riskset')
	rs[,1] <- as.numeric(rs[,1])
	rs[,2] <- as.numeric(rs[,2])

	if(end_undirected){
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
		
		if (nrow(edgelist) > 0) {
		 	engaged_stats <- tabulate(as.numeric(unlist(edgelist[edgelist[,1] <= t & edgelist[,4] > t, 2:3], use.names = FALSE)), nbins = num_actors) + tabulate(as.numeric(unlist(edgelist[is.na(edgelist[,4]), 2:3], use.names = FALSE)), nbins = num_actors)		
		} else {
			engaged_stats <- rep(0, num_actors) # No engaged events
		}

		#updating event rate / lambda
		lambda <- sapply(1:nrow(rs),function(x){
			if(engaged_stat){
				if(dyads_is.active[x]){
					if(end_count == 0){
						return(exp(end_params[1])) #just baseline effect for first event
					}				
				if(dim(end_stats)[3]==1){
					if(length(end_engaged_params)==1){ #engaged undirected
						return(exp(end_stats[1,x,] * end_params +(engaged_stats[rs[x,1]]+engaged_stats[rs[x,2]])*end_engaged_params[1] ))
					}
					return(exp(end_stats[1,x,] * end_params +(engaged_stats[rs[x,1]]*end_engaged_params[1] + engaged_stats[rs[x,2]]*end_engaged_params[2])))
				}
				if(length(end_engaged_params)==1){ #engaged undirected
						return(exp(end_stats[1,x,] %*% end_params +(engaged_stats[rs[x,1]]+engaged_stats[rs[x,2]])*end_engaged_params[1] ))
					}
				return(exp(end_stats[1,x,] %*% end_params +(engaged_stats[rs[x,1]]*end_engaged_params[1] + engaged_stats[rs[x,2]]*end_engaged_params[2]))) 
			}else{
				if(start_count == 0){
					return(exp(start_params[1])) #just baseline effect for first event
				}
				if(dim(start_stats)[3]==1){
					if(length(start_engaged_params)==1){ #engaged undirected
						return(exp(start_stats[1,x,] * start_params +(engaged_stats[rs[x,1]]+engaged_stats[rs[x,2]])*start_engaged_params[1]))
					}
					return(exp( (start_stats[1,x,] * start_params) +(engaged_stats[rs[x,1]]*start_engaged_params[1] + engaged_stats[rs[x,2]]*start_engaged_params[2])))
				}			
				if(length(start_engaged_params)==1){ #engaged undirected
						return(exp(start_stats[1,x,] %*% start_params +(engaged_stats[rs[x,1]]+engaged_stats[rs[x,2]])*start_engaged_params[1]))
					}
				return(exp( (start_stats[1,x,] %*% start_params)  +(engaged_stats[rs[x,1]]*start_engaged_params[1] + engaged_stats[rs[x,2]]*start_engaged_params[2])))
			}
			}else{
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
			if(!end_undirected){
				end_reh <- remify::remify(end_weighted_edgelist, actors = 1:num_actors,model="tie")
			}else{
				end_reh <- remify::remify(end_weighted_edgelist, actors = 1:num_actors,model="tie",directed=FALSE)
			}
			
			})
			end_stats <- remstats::tomstats(effects = end_effects, reh = end_reh, start = length(end_indx), stop = length(end_indx))
			
		}else{

			#start an event
			if(start_count < num_events){
				#update dyad/sender/receiver
				evls = rbind(evls,data.frame(dyad = dyad, start_time = t, end_time = NA))
				edgelist = rbind(edgelist,data.frame(start_time = t, actor1 = rs[dyad,1], actor2 = rs[dyad,2], end_time = NA))
				
				start_weighted_edgelist = rbind(start_weighted_edgelist,data.frame(time = t, actor1 = rs[dyad,1], actor2 = rs[dyad,2], weight = 1))
				end_weighted_edgelist = rbind(end_weighted_edgelist,data.frame(time = NA, actor1 = rs[dyad,1], actor2 = rs[dyad,2], weight = 0))
				dyads_is.active[dyad] = TRUE				
			}     
			#weighted edgelist
			#start_reh <- remify::remify(edgelist[,c(1,2,3,5), actors = 1:num_actors,model="tie")			
			if(!start_undirected){
				start_reh <- remify::remify(start_weighted_edgelist, actors = 1:num_actors,model="tie")			
			}else{
				start_reh <- remify::remify(start_weighted_edgelist, actors = 1:num_actors,model="tie",directed = FALSE)			
			}	

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