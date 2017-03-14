## Some basic step-growth polymerization results:
SGP_Plot_Number_Fractions <- function(p,N_max)
{
  Ns <- seq(from=1,to=N_max,by=1)
  ps <- array(0,dim=c(length(Ns)))
  for(i in 1:N_max)
  {
    ps[i] <- (1.0-p)*p^(i-1)
  }
  df_nfs <- data.frame(Ns,ps)
  return(df_nfs)
}

SGP_Plot_Number_Fractions_For_Time <- function(reac_time,kp,c0,N_max)
{
  p <- (kp*c0*reac_time)/(1.0+kp*c0*reac_time)
  return(SGP_Plot_Number_Fractions(p,N_max));
}

## Simulate_Step_Growth_Polym(init_numbers,num_iters,N_max,T_max)
## simulates the time evolution of a system 
Simulate_Step_Growth_Polym <- function(init_numbers,num_iters,N_max,T_max)
{
  ## calculate the number of reactions:
  num_reactions <- N_max*(N_max+1.0)/2.0;
  print("number of reactions: ")
  print(num_reactions)
  ## calculate the state-change vectors:
  scv <- matrix()
}


Simulate_System_With_Backreactions <- function(init_particle_numbers,N_max,num_iters,thinning_factor)
{
  ## number of on-reactions
  num_on_reactions <- N_max*(N_max+1)/2
  ## number of off-reactions
  num_off_reactions <- N_max*(N_max-1)/2
  num_total_reactions <- num_on_reactions+num_off_reactions
  
  ## initialize the number of particles:
  Ns <- array(0,dim=c(N_max))
  Ns <- init_particle_numbers
  
  Ns_record <- array(0,dim=c(num_iters/thinning_factor,N_max))
  time_record <- array(0,dim=c(num_iters/thinning_factor))
##  Ns_record[1,] <- Ns
  
  cs_on_0 <- 2.66e-8
  cs_off_0 <- 3e-4
  
  ## initialize the on-rates:
  cs_on <- array(cs_on_0,dim=c(num_on_reactions))
  cs_off <- array(cs_off_0,dim=c(num_off_reactions))
  
  ## create state change vector (one for each reaction)
  scv <- array(0,dim=c(num_total_reactions,N_max))
  counter <- 0
  for(i in 1:N_max)
  {
    
    for(j in i:N_max)
    {
      curr_scv <- array(0,dim=c(N_max))
      curr_scv[i] <- curr_scv[i]-1
      curr_scv[j] <- curr_scv[j]-1
      if((i+j)<=N_max)
        curr_scv[i+j] <- 1
      counter <- counter+1
      scv[counter,] <- curr_scv
    }
    
  }
  ## state change vectors for off-reactions
  counter <- num_on_reactions
  for(i in 2:N_max)
  {
    for(j in 1:(i-1))
    {
      curr_scv <- array(0,dim=c(N_max))
      ## an i-mer, so indices will be j and i-j
      curr_scv[i] <- curr_scv[i] - 1
      curr_scv[j] <- curr_scv[j] + 1
      curr_scv[i-j] <- curr_scv[i-j] + 1
      counter <- counter+1
      scv[counter,] <- curr_scv
    }
  }
  
  ## create the propensity functions:
  ## propensity functions for on-reactions
  ## propensity functions for off-reactions (first-order reactions...)
 
  print(scv)
  print(as)
  
  ## start the simulation (iterate)
  test <- array(0,dim=c(num_iters))
  curr_time <- 0
  curr_index <- 1
  for(iter in 1:num_iters)
  {
    counter <- 0
    as <- array(0,dim=c(num_total_reactions))
    for(i in 1:N_max)
    {
      
      for(j in i:N_max)
      {
        counter <- counter+1
        if(i==j)
        {
          as[counter] <- 0.5*cs_on[counter]*Ns[i]*(Ns[i]-1)
        }
        else
        {
          as[counter] <- cs_on[counter]*Ns[i]*Ns[j]
        }
      }
    }
    counter <- num_on_reactions
    for(i in 2:N_max)
    {
      for(j in 1:(i-1))
      {
        counter <- counter+1
        as[counter] <- cs_off[counter-num_on_reactions]*Ns[i] 
      }
    }
    a0 <- sum(as)
    tau <- rexp(n=1,rate=a0)
    j <- sample.int(n=num_total_reactions,size=1,replace=TRUE,prob=as)
    curr_time <- curr_time + tau
    Ns <- Ns  + scv[j,]
    if(j>num_on_reactions)
      print("YUHU, off reaction took place")
    if(iter%%thinning_factor==1)
    {
      print(curr_index)
      Ns_record[curr_index,] <- Ns
      time_record[curr_index] <- curr_time
      print(Ns_record[curr_index,])
      curr_index <- curr_index + 1
    }
  }
  result_list <- list(times=time_record,Ns=Ns_record)
  return(result_list)
}

Test_Simulation <- function()
{
  N_max <- 50
  init_numbers <- array(0,dim=c(N_max))
  init_numbers[1] <- 241000
  num_iters <- 200000
  thinning_factor <- 100
  res_list <- Simulate_System_With_Backreactions_And_Length_Dependence(init_numbers,N_max,num_iters,thinning_factor)
  return(res_list)
}

Get_Initial_Numbers_And_Rates <- function()
{
  volume_in_m3 <- 1e-13 ## one pico-cubic-metre
  volume_in_litre <- volume_in_m3*1000
  
  des_conc <- 4e-9 ## 4 nM
  
  ## how many particles are in this volume ?
  particles_per_mole <- 6.022e23
  des_particles_per_litre <- particles_per_mole*des_conc
  des_particles_in_volume <- des_particles_per_litre*volume_in_litre
  
  print("particles per litre: ")
  print(des_particles_per_litre)
  
  print("particles in the set volume: ")
  print(des_particles_in_volume)
  
  ## get an estimate for the c_j,.....
  k_on <- 8e5/particles_per_mole ## 1/(s*(mol/l))
  ##k_on <- 8e8 ## 1/(s*(mol/m^3))
  print("desired c_on:")
  c_on <- 2*k_on / volume_in_litre
  print(c_on)
}

Plot_Degree_Of_Polymerization <- function(result_list)
{
  num_rows <- nrow(result_list$times)
  num_of_polym_sites <- array(0,dim=c(num_rows))
  for(i in 1:num_rows)
  {
    curr_dop <-0
    for(j in 2:ncol(result_list$Ns))
    {
      curr_dop <- curr_dop + result_list$Ns[i,j]*(j-1)
    }
    num_of_polym_sites[i] <- curr_dop
  }
  df_res <- data.frame(ts=result_list$times,dop=num_of_polym_sites)
  ##plot(result_list$times,num_of_polym_sites)
  return(df_res)
}

Plot_Degree_OF_Polymerization_theor <- function()
{
  ## degree of polymerization is given by:
  ## X_mean =  1/(1-p) 
  ## with p(t) = (k_p c0 t)/(1+k_p c0 t)
  
  times <- seq(from=0,to=10000,by=2)
  k_p <- 8e5;
  c0 <-  4e-9;
  num_points <- length(times)
  ps <- array(0,dim=c(num_points))
  dops <- array(0,dim=c(num_points))
  for(i in 1:num_points)
  {
    ps[i] <-  k_p*c0*times[i]/(1.0+k_p*c0*times[i])
    dops[i] <- 1.0/(1.0 - ps[i]);
  }
  plot(times,dops,type="l")
}

Plot_Mean_Length_Over_Time <- function(result_list)
{
  num_rows <- nrow(result_list$times)
  mean_length <- array(0,dim=c(num_rows))
  for(i in 1:num_rows)
  {
    curr_ml <-0
    
    for(j in 1:ncol(result_list$Ns))
    {
      norm <- norm+
      curr_ml <- curr_ml + result_list$Ns[i,j]*j^2
    }
    mean_length[i] <- curr_ml
  }
  df_res <- data.frame(ts=result_list$times,ml=mean_length)
  ##plot(result_list$times,num_of_polym_sites)
  return(df_res)
}

## make plots of number fraction distributions:
Make_Number_fraction_distribution_plots <- function(res_list)
{
  indices <- c(391,1049,1789)
  norm_res_list <- res_list$Ns
  times <- res_list$times
  colors=c("black","red","green")
  for(i in 1:length(indices))
  {
    print("time: ")
    print(times[indices[i]])
    norm_res_list[indices[i],] <- norm_res_list[indices[i],]/sum(norm_res_list[indices[i],])
    if(i==1)
    {
      plot(norm_res_list[indices[i],],xlab="N-mer",ylab="relative Häufigkeit",xlim=c(0,12),col=colors[i],pch=19,ylim=c(0,1),cex=2)
    }
    else
    {
      points(norm_res_list[indices[i],],xlab="N-mer",ylab="relative Häufigkeit",xlim=c(0,12),col=colors[i],pch=19,cex=2)
    }
    ## do the simple minded simulation:
    df_theor <- SGP_Plot_Number_Fractions_For_Time(times[indices[i]],kp=8e5,c0=4e-9,N_max=20)
    points(df_theor$Ns,df_theor$ps,col=colors[i],type="l",lwd=2)
    ##legend("top",c("1 Minute","4 Minuten", "15 Minuten"),text.col=c("black","red","green"),box.lty=0)
    ##legend("topright",inset=.02,title="Polymerisierungszeit in Minuten", c("1","4","15"),text.col=c("black","red","green"),horiz=TRUE,box.lty=0)
    ##legend("right",c("SGP-Modell","Simulation"),box.lty=0,pch=c(19,NA),lty=c(NA,1))
    legend("topright",c("1 Minute (SGP)","4 Minuten (SGP)","15 Minuten (SGP)","1 Minute (Sim.)","4 Minuten (Sim.)","15 Minuten (Sim.)"),box.lty=0,pch=c(19,19,19,NA,NA,NA),lty=c(NA,NA,NA,1,1,1),
           col=c("black","red","green","black","red","green"))
  }
}

## Make_Mean_Polymerlength_Plot
Make_Mean_Polymerlength_Plot <- function(res_list1,res_list2,res_list3)
{
  num_rows <- nrow(res_list1$Ns)
  print("number of times: ")
  print(num_rows)
  num_ns <- ncol(res_list$Ns)
  
  df_pl_dist1 <- Get_Mean_Polymerlength(res_list1)
  df_pl_dist2 <- Get_Mean_Polymerlength(res_list2)
  df_pl_dist3 <- Get_Mean_Polymerlength(res_list3)
  
  plot(df_pl_dist1$times,df_pl_dist1$mpls,type="l",col="black",xlab="Zeit [s]",ylab="mittlere Polymerlänge",xlim=c(0,1100),ylim=c(0,5),lwd=3)
  points(df_pl_dist2$times,df_pl_dist2$mpls,type="l",col="red",lwd=3)
  points(df_pl_dist3$times,df_pl_dist3$mpls,type="l",col="green",lwd=3)
  
  legend("topleft",c("c_off = 0 1/s","c_off = 3e-4 1/s","c_off = 1e-3 1/s"),col=c("black","red","green"),box.lty=0,lty=c(1,1,1),lwd=3)
  points(c(60,60),c(0,1.2),type="l",lty=2)
  points(c(300,300),c(0,1.95),type="l",lty=2)
  points(c(900,900),c(0,3.85),type="l",lty=2)
  
}

## Get_Mean_Polymerlength <- function(rl)
Get_Mean_Polymerlength <- function(rl)
{
  num_times <- nrow(rl$Ns)
  num_Ns <- ncol(rl$Ns);
  print("number of times: ")
  print(num_times)
  print("number of NS:");
  print(num_Ns)
  mpl <- array(0,dim=c(num_times))
  for(i in 1:num_times)
  {
    ml <- 0
     
    for(j in 1:ncol(rl$Ns))
    {
      ml <- ml+rl$Ns[i,j]*j
    }
    ml <- ml/sum(rl$Ns[i,])
    mpl[i] <- ml
  }
  df_pl_dist <- data.frame(times=rl$times,mpls=mpl)
  return(df_pl_dist)
}

## Make_Comparative_Plot_Of_Number_fraction_distribution(res_list1,res_list2,res_list3)
Make_Comparative_Plot_Of_Number_fraction_distribution <- function (time_value,res_list1,res_list2,res_list3)
{
  ## search for indices of time:
  indices1 <- findInterval(x=time_value,vec=unlist(res_list1$times))
  indices2 <- findInterval(x=time_value,vec=unlist(res_list2$times))
  indices3 <- findInterval(x=time_value,vec=unlist(res_list3$times))
  
  number_fractions_dist1 <- res_list1$Ns[indices1,]/sum(res_list1$Ns[indices1,])
  number_fractions_dist2 <- res_list2$Ns[indices2,]/sum(res_list2$Ns[indices2,])
  number_fractions_dist3 <- res_list3$Ns[indices3,]/sum(res_list3$Ns[indices3,])
  
  colors = c("black","red","green")
  
  plot(number_fractions_dist1,xlab="N-mer",ylab="relative Häufigkeit",xlim=c(0,12),col=colors[1],type="o",lty=1,lwd=2,pch=19,ylim=c(0,1))
  points(number_fractions_dist2,col=colors[2],type="o",lty=1,lwd=2,pch=19)
  points(number_fractions_dist3,col=colors[3],type="o",lty=1,lwd=2,pch=19)
  legend("topright",c("c_off = 0 1/s","c_off = 3e-4 1/s","c_off = 1e-3 1/s"),box.lty=0,lty=1,pch=19,
         col=c("black","red","green"))
}


## Simulate_System_With_Backreactions_And_Length dependence
Simulate_System_With_Backreactions_And_Length_Dependence <- function(init_particle_numbers,N_max,num_iters,thinning_factor)
{
  ## number of on-reactions
  num_on_reactions <- N_max*(N_max+1)/2
  ## number of off-reactions
  num_off_reactions <- N_max*(N_max-1)/2
  num_total_reactions <- num_on_reactions+num_off_reactions
  
  ## initialize the number of particles:
  Ns <- array(0,dim=c(N_max))
  Ns <- init_particle_numbers
  
  Ns_record <- array(0,dim=c(num_iters/thinning_factor,N_max))
  time_record <- array(0,dim=c(num_iters/thinning_factor))
  ##  Ns_record[1,] <- Ns
  
  cs_on_0 <- 2.66e-8
  cs_off_0 <- 0
  delta <- 0.5 ## parameter for power law dependence of on-rate on polymer length...
  
  ## initialize the on-rates:
  cs_on <- array(cs_on_0,dim=c(num_on_reactions))
  cs_off <- array(cs_off_0,dim=c(num_off_reactions))
  
  ## initialize the on-rates (case of length-dependence)
  counter <- 0
  for(i in 1:N_max)
  {
    for(j in i:N_max)
    {
      counter <- counter+1
      cs_on[counter] <- cs_on_0/(1+i^(delta)) + cs_on_0/(1+j^(delta));
      print(cs_on[counter])
    }
  }
  
  ## create state change vector (one for each reaction)
  scv <- array(0,dim=c(num_total_reactions,N_max))
  counter <- 0
  for(i in 1:N_max)
  {
    for(j in i:N_max)
    {
      curr_scv <- array(0,dim=c(N_max))
      curr_scv[i] <- curr_scv[i]-1
      curr_scv[j] <- curr_scv[j]-1
      if((i+j)<=N_max)
        curr_scv[i+j] <- 1
      counter <- counter+1
      scv[counter,] <- curr_scv
    }
  }
  
  ## state change vectors for off-reactions
  counter <- num_on_reactions
  for(i in 2:N_max)
  {
    for(j in 1:(i-1))
    {
      curr_scv <- array(0,dim=c(N_max))
      ## an i-mer, so indices will be j and i-j
      curr_scv[i] <- curr_scv[i] - 1
      curr_scv[j] <- curr_scv[j] + 1
      curr_scv[i-j] <- curr_scv[i-j] + 1
      counter <- counter+1
      scv[counter,] <- curr_scv
    }
  }
  
  ## create the propensity functions:
  ## propensity functions for on-reactions
  ## propensity functions for off-reactions (first-order reactions...)
  print(scv)
  print(as)
  
  ## start the simulation (iterate)
  test <- array(0,dim=c(num_iters))
  curr_time <- 0
  curr_index <- 1
  for(iter in 1:num_iters)
  {
    counter <- 0
    as <- array(0,dim=c(num_total_reactions))
    for(i in 1:N_max)
    {
      
      for(j in i:N_max)
      {
        counter <- counter+1
        if(i==j)
        {
          as[counter] <- 0.5*cs_on[counter]*Ns[i]*(Ns[i]-1)
        }
        else
        {
          as[counter] <- cs_on[counter]*Ns[i]*Ns[j]
        }
      }
    }
    counter <- num_on_reactions
    for(i in 2:N_max)
    {
      for(j in 1:(i-1))
      {
        counter <- counter+1
        as[counter] <- cs_off[counter-num_on_reactions]*Ns[i] 
      }
    }
    a0 <- sum(as)
    tau <- rexp(n=1,rate=a0)
    j <- sample.int(n=num_total_reactions,size=1,replace=TRUE,prob=as)
    curr_time <- curr_time + tau
    Ns <- Ns  + scv[j,]
    if(j>num_on_reactions)
      print("YUHU, off reaction took place")
    if(iter%%thinning_factor==1)
    {
      print(curr_index)
      Ns_record[curr_index,] <- Ns
      time_record[curr_index] <- curr_time
      print(Ns_record[curr_index,])
      curr_index <- curr_index + 1
    }
  }
  result_list <- list(times=time_record,Ns=Ns_record)
  return(result_list)
}

