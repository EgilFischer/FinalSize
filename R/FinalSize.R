######################################################################################
#                                                                                    #
#                   Final size estimation
#                   with a fast implementation                 #
#                                                                                    #
#                   Copyright R-code: Egil A.J.Fischer                             #
#                   e.a.j.fischer@uu.nl/egil@egilfischer.nl
######################################################################################
require(ggplot2)
###
#' @title Final Size distribution for given R and intitial values
#'
#' @description
#' this function will first create a generation based table with states and then gives the final size distribution as output
#' @param R Value of R0
#' @param s0in Number of initial susceptibles
#' @param i0in Number of initial infectious
#' @param r0in Number of initial recovered
#'
#' @return distribution of number of individuals that were infected during the outbreak
#'
#' @examples
#' distFS(1, 3, 3, 0)
#' 
#'          [,1]      [,2]      [,3]      [,4]
#' [1,] 0.2962963 0.2825521 0.2505569 0.1705947
#' 
#' distFS(2, 3, 3, 0)  
#'       [,1]   [,2]      [,3]      [,4]
#' [1,] 0.125 0.1638 0.2552766 0.4559234
#' 
#' distFS(2, 3, 3, 1)  
#'         [,1]      [,2]      [,3]      [,4]
#'[1,] 0.156122 0.1939235 0.2688927 0.3810617
#' 
distFS <- function(R,s0in,i0in ,r0in = NULL)
  {

    #set size of the population
    if(is.null(r0in)){r0in <- 0*s0in}
    nin = s0in+i0in + r0in;


    #infection probability
    pinf =function(s,n){s*R/(s*R + n)};

     #create final size distribution for each s0 and i0 combination
     max.cases <- max(s0in);
     gen.mat.out = matrix(0,ncol = (max.cases+1),nrow = length(s0in));
     for(m in 1:length(s0in))
       {
       s0= s0in[m];
       i0= i0in[m];
       ifelse(is.null(r0in),r0 <- 0,r0 <- r0in[m] )
       n = nin[m];
       #set the state matrix
       gen.mat <- matrix(rep(0, (s0+i0+1)*(s0+1)),nrow= (s0+1));
       gen.mat[s0+1,i0+1]<- 1;
       #fill first row
       for(j in c((i0+1):1)){
         gen.mat[s0+1,j-1]<- gen.mat[s0+1,j]*(1-pinf(s0,n))
         }

        #next rows added by infection
        for(i in c(1:(s0+1))){
          for(j in c((s0+i0+1):1)){
            #enter state by infection (only for j > 2, j = 1 no infections, j = 2 is 1 infection)
            if(j>2){
              gen.mat[s0+1-i,j]<-gen.mat[s0+1-i,j] + gen.mat[s0+2-i,j-1]*pinf(s0+1-i,n)
             };
          if(j < s0+i0+1){
            #enter state by recovery
            gen.mat[s0+1-i,j]<- gen.mat[s0+1-i,j] + gen.mat[s0+1-i,j+1]*(1-pinf(s0-i,n));
            }

          }
        }
    #store the final size distribution -> reversed so that it is 0 cases first
    gen.mat.out[m,1:(s0+1)] = rev(gen.mat[,1]);
     }
     #report it reversed
     return(gen.mat.out)
}

#plot final size distribution given R, initial susceptibles, initial infectious and initial recovered
plotFSdist <- function(R,s0in,i0in ,r0in = NULL){
  n <- s0in + i0in + ifelse(is.null(r0in),0,r0in);
  input <- data.frame(final.size = c(1:(n-r0in)),
                   prob = c(distFS(R,s0in, i0in,r0in)));
  return(ggplot(data = input)+
           geom_col(aes(x = final.size, y = prob)))


}

#, initial susceptibles, initial infectious and initial recovered
#' @title probability of final size given R
#'
#' @param R Value of R0
#' @param x Number of cases at the end of the outbreak
#' @param s0in Number of initial susceptibles
#' @param i0in Number of initial infectious
#' @param r0in Number of initial recovered
#'
#' @return probability to observe x given other values
#'
#' @examples
#' pFS(2,0,1,1)
#' pFS(2,1,1,1)
#' 
pFS <- function(R,x,s0,i0){
  
  #produce final size distribution for value R
  final.size.dist <- distFS(R,s0,i0);
  #return the outcome for x
  return(prod(mapply(function(i,j)final.size.dist[i,j+1],c(1:length(s0)),x)))
}

###Function to determine the probability of more extreme values given R####
#' @title Probability of extreme values given R
#'
#' @param R Value of R0
#' @param x Number of cases at the end of the outbreak
#' @param s0in Number of initial susceptibles
#' @param i0in Number of initial infectious
#' @param r0in Number of initial recovered
#' @param comp Direction of comparison i.e.(`>`or`<`)
#'
#' @return Exact probability that x cases or more extreme are found for this value of R
#'
#' @examples
#' pExtremes(1, 1, 5,5, `<`)
#' pExtremes(1, 1, 5,5, `>`)
#' 
pExtremes<-  function(R,x,s0in,i0in,comp = `<`){
  
  #create all possible outcomes of these transmission experiments
  #this means all possibilities between 0 and s0 contact infection (hence s0 + 1 options per trial)
  out <- matrix(ncol = length(s0in),nrow = prod(s0in+1))
  #repeat the outcome as many times as the previous trial possibilities
  for(k in c(1:length(s0in)))
  {
    #check how often to repeat the same number given previous trials
    repetition <- ifelse(k > 1,prod(s0in[1:k-1]+1),1)
    #put it in the matrix
    out[,k]<- matrix(sapply(c(0:s0in[k]),
                            FUN = function(x){rep(x,repetition)}),
                     ncol = 1,
                     nrow =prod(s0in+1))[,1]
  }
  #produce final size distribution for value R
  final.size.dist <- distFS(R,s0in,i0in);
  #define function for this distribution for the probability of a certain number of cases x in each of the trials
  pFSloc<- function(v){
    return(prod(mapply(function(i,j) {final.size.dist[i,j+1]},
                       c(1:length(s0in)),
                       v)))
    }
  #select the extremes by selecting those for which the total number of cases
  #and the sum of x are given by the comparison "comp"  thus either <,>,<= or >=
  #calculate for each extreme the probability of the Final Size under the hypothesis R = r
  #and sum all probabilities of these extreme outcomes
  #if only one experiment
  if(length(s0in)==1)
  {
    return(sum(final.size.dist[out[comp(apply(out,1,sum),sum(x)),]+1]))
  }
  #If only one extreme outcome possible:
   if(is.null(dim(out[comp(apply(out,1,sum),sum(x)),]))){
     #only requires the final size distribution
     return(pFSloc(out[comp(apply(out,1,sum),sum(x)),]))
     }
    #else
       return(c(sum(apply(out[comp(apply(out,1,sum),sum(x)),],
                          1,
                          function(ext){pFSloc(ext)}))))
}

##################################################################################################
#                                                                                                #
#                  Function to estimate R using the final size method                            #
#                  Includes confidence interval size 1- alpha and R >= 1 test                    #
#                                                                                                #
##################################################################################################
# Final Size function for x cases, s0 initially susceptible and i0 initially infectious animals.
# Optional set sigficance level alpha, onesided testing
# max.val is the maximum value used for optimization.
#' @title Final size analysis
#'
#' @param x Number of cases at the end of the outbreak (single value or vector) 
#' @param s0 Number of initial susceptibles (single value or vector)
#' @param i0 Number of initial infectious (single value or vector)
#' @param comp Direction of comparison i.e.(`>`or`<`)#' @param i0 
#' @param alpha Significance level default = 0.05
#' @param onesided Test onesided default = FALSE
#' @param max.val Maximum value of R for numeric optimization
#' @param decimals Number of decimals
#'
#' @return Table with estimator, confidence interval and test R<1, R>1 and R = 1
#'  point.est     ci.ll     ci.ul pval.above1 pval.below1 pval.equal1
#' @examples
#' #tests for same code as Mathematica book of Mart de Jong
#' FinalSize(c(0,0,0,0),c(2,2,2,2),c(2,2,2,2),onesided = T)
#' FinalSize(c(2,2,2,2),c(2,2,2,2),c(2,2,2,2),onesided = T)
#' FinalSize(c(2,1,1,0),c(2,2,2,2),c(2,2,2,2))
#' FinalSize(c(3),c(20),c(20),max.val = 50)
#' 
#' FinalSize(c(7,4),c(40,19)-1, c(1,1))
#' FinalSize(c(5,5),c(40,40)-1, c(1,1))
#' FinalSize(c(7),c(40)-1, c(1))
#' FinalSize(c(6,2),c(24,12)-1, c(1,1))
#' 
FinalSize.Estimator<- function(x,s0,i0, alpha = 0.05, onesided = FALSE, max.val = 250, decimals =2){
  #check data consistency
  if(min(s0-x)<0)stop("more cases x than susceptibles s0")
  #
  res <- data.frame(point.est = -999)
  #determine the point estimate by optimization of the log-likelihood function
  res$point.est <- ifelse(sum(x)==0,0,
                          ifelse(sum(x)==sum(s0),Inf,
                    optimize(interval = c(0.,max.val),
                            f = function(R.optimize){-log(pFS(R.optimize,x,s0,i0))})$minimum))
  #determine the confidence intervals
  #if one-sided is FALSE both sides, either only lower or upper limit of CI
  #lowerlimit is found for values of R for which the probability of extremes below the observations
  res$ci.ll <- ifelse(sum(x)==0,0,
                      uniroot(interval = c(10^-10,max.val),extendInt = "yes",
                              f = function(R.optimize){(pExtremes(R.optimize,x,s0,i0,comp = `>=`) - alpha / (2 - onesided))})$root)
  #upperlimit is found for values of R for which the probability of extremes above the observations
  res$ci.ul <- ifelse(sum(x)==sum(s0),Inf,
                      uniroot(interval = c(0,max.val),extendInt = "yes",
                              f = function(R.optimize){( pExtremes(R.optimize,x,s0,i0,comp = `<=`) - alpha / (2 - onesided))})$root)

  #probability of R >= 1 is found be calculating the probability to find an equal or less positive under the assumption R0 = 1
  res$pval.above1 = pExtremes(1,x,s0,i0,comp = `<=`)
  res$pval.below1 = pExtremes(1,x,s0,i0,comp = `>=`)
  res$pval.equal1 = pExtremes(1,x,s0,i0,comp = `>=`) * pExtremes(1,x,s0,i0,comp = `<=`)
  return(round(res, digits = decimals))
}

#
FinalSize<- function(x,s0,i0,group = NULL ,alpha = 0.05, onesided = FALSE, max.val = 250, decimals =2){
  if(is.null(group))  {
    return(FinalSize.Estimator(x,s0,i0, alpha, onesided, max.val, decimals))
  }else{res <-NULL;
    for(j in unique(group)){
      #iterate over groups
      res <-rbind(res,cbind(Group = j,FinalSize.Estimator(x[group == j],s0[group == j],i0[group == j], alpha, onesided, max.val, decimals) ))
      }
    }
    return(res);
}


