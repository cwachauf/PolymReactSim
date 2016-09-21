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