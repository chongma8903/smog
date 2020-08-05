# Copyright (c) 2018 - 2020 Chong Ma
 
# This file contains the kernel function for the R package smog. 
# The function smog is written for the generalized linear model constraint 
# on specified hierarchical structures by using overlapped group penalty. 
# It is implemented by combining the ISTA and ADMM algorithms, and works 
# for continuous, multimonial and survival data. 


#' Generalized linear model constraint on hierarchical structure
#' by using overlapped group penalty
#' 
#' \code{smog} fits a linear non-penalized phynotype (demographic) variables such as 
#' age, gender, treatment, etc, and penalized groups of prognostic effect (main effect)
#' and predictive effect (interaction effect), by satisfying the hierarchy structure:
#' if a predictive effect exists, its prognostic effect must be in the model. It can deal
#' with continuous, binomial or multinomial, and survival response variables, underlying 
#' the assumption of Gaussian, binomial (multinomial), and Cox proportional hazard models,
#' respectively. It can accept \code{\link[stats]{formula}}, and output coefficients table,
#' fitted.values, and convergence information produced in the algorithm iterations.   
#' 
#' @param x a model matrix, or a data frame of dimensions n by p, 
#'          in which the columns represents the predictor variables. 
#' @param y response variable, corresponds to the family description. 
#'          When family is ``gaussian'' or ``binomial'', \code{y} ought to
#'          be a numeric vector of observations of length n; when family 
#'          is ``coxph'', \code{y} represents the survival objects, containing the 
#'          survival time and the censoring status. See \code{\link[survival]{Surv}}.
#' @param g a vector of group labels for the predictor variables.
#' @param v a vector of binary values, represents whether or not the 
#'          predictor variables are penalized. Note that 1 indicates 
#'          penalization and 0 for not penalization.
#' @param label a character vector, represents the type of predictors in terms of treatment,
#'              prognostic, and predictive effects by using ``t'', ``prog'', and ``pred'',
#'              respectively.
#' @param lambda1 penalty parameter for the L2 norm of each group of prognostic and predictive effects.
#' @param lambda2 ridge penalty parameter for the squared L2 norm of each group of prognostic and predictive effects.
#' @param lambda3 penalty parameter for the L1 norm of predictive effects.              
#' @param family a description of the distribution family for the response 
#'               variable variable. For continuous response variable,
#'               family is ``gaussian''; for multinomial or binary response
#'               variable, family is ``binomial''; for survival response
#'               variable, family is ``coxph'', respectively.
#' @param subset an optional vector specifying a subset of observations to be 
#'               used in the model fitting. Default is \code{NULL}.
#' @param rho   the penalty parameter used in the alternating direction method 
#'              of multipliers (ADMM) algorithm. Default is 10.
#' @param scale whether or not scale the design matrix. Default is \code{TRUE}.
#' @param eabs  the absolute tolerance used in the ADMM algorithm. Default is 1e-3.
#' @param erel  the reletive tolerance used in the ADMM algorithm. Default is 1e-3.
#' @param LL    initial value for the Lipschitz continuous constant for 
#'              approximation to the objective function in the Majorization-
#'              Minimization (MM) (or iterative shrinkage-thresholding algorithm 
#'              (ISTA)). Default is 1.
#' @param eta   gradient stepsize for the backtrack line search for the Lipschitz
#'              continuous constant. Default is 1.25. 
#' @param maxitr the maximum iterations for convergence in the ADMM algorithm. 
#'               Default is 1000.
#' @param formula an object of class ``formula'': a symbolic description of the
#'                model to be fitted. Should not include the intercept. 
#' @param data    an optional data frame, containing the variables in the model. 
#' @param ...   other relevant arguments that can be supplied to smog.
#' 
#' @return \code{smog} returns an object of class inhering from ``smog''. The 
#'         generic accessor functions \code{coef}, \code{coefficients}, 
#'         \code{fitted.value}, and \code{predict} can be used to extract
#'         various useful features of the value returned by \code{smog}.
#'         
#'         An object of ``smog'' is a list containing at least the following 
#'         components: 
#'         
#'         \item{coefficients}{a data frame containing the nonzero predictor
#'                             variables' indexes, names, and estimates. When
#'                             family is ``binomial'', the estimates have K-1 
#'                             columns, each column representing the weights for the 
#'                             corresponding group. The last group behaves the
#'                             ``pivot''.}
#'         \item{fitted.values}{the fitted mean values for the response variable,
#'                              for family is ``gaussian''. When family is 
#'                              ``binomial", the fitted.values are the probabilies
#'                              for each class; when family is ``coxph'', 
#'                              the fitted.values are survival probabilities.}
#'         \item{model}{a list of estimates for the intercept, treatment effect, 
#'                      and prognostic and predictive effects for the selectd
#'                      biomarkers.}
#'         \item{weight}{the weight of predictors resulted from the penalty funciton,
#'                       is used to calculate the degrees of freedom.}
#'         \item{DF}{the degrees of freedom. When family = ``gaussian'', 
#'                   \eqn{DF = tr(x_{\lambda}'(x_{\lambda}'x_{\lambda}+W)x_{\lambda})}. 
#'                   For other families, DF is approximated by \eqn{diag(1/(1+W))}.}
#'         \item{criteria}{model selection criteria, including the correction Akaike's Information 
#'                         Criterion (AIC), AIC, Bayesian Information Criterion (BIC), and the generalized 
#'                         cross-validation score (GCV), respectively. See also \code{\link{cv.smog}}.
#'                         \describe{
#'                         \item{cAIC}{\eqn{\frac{n}{2}}log(|2*log-likelihood|) + 
#'                                     \eqn{\frac{n}{2} \left( \frac{1+k/n}{1-k+2/n} \right)}.}
#'                         \item{AIC}{log(|2*log-likelihood|/n) + \eqn{2\frac{k}{n}}.}
#'                         \item{BIC}{log(|2*log-likelihood|/n) + \eqn{2\frac{k}{n}}log(n).}
#'                         \item{GCV}{|2*log-likelihood| / \eqn{(n(1-k/n)^2)}.}
#'                         }
#'                         Where k is the degrees of freedom \code{DF}, which is related to the 
#'                         penalty parameters \eqn{\lambda}'s. 
#'                         }
#'         \item{llikelihood}{the log-likelihood value for the converged model.}
#'         \item{loglike}{the penalized log-likelihood values for each 
#'                        iteration in the algorithm.}
#'         \item{PrimalError}{the averged norms \eqn{||\beta-Z||/\sqrt{p}} for each iteration,
#'                            in the ADMM algorithm.}
#'         \item{DualError}{the averaged norms \eqn{||Z^{t+1}-Z^{t}||/\sqrt{p}} for 
#'                          each iteration, in the ADMM algorithm.}
#'         \item{converge}{the number of iterations processed in the ADMM algorithm.}
#'         \item{call}{the matched call.}
#'         \item{formula}{the formula supplied.}
#'    
#' 
#' @details The formula has the form \code{response ~ 0 + terms} where \code{terms} is
#'          a series of predictor variables to be fitted for \code{response}. For \code{gaussian} 
#'          family, the response is a continuous vector. For \code{binomial} family, 
#'          the response is a factor vector, in which the last level denotes the ``pivot''.
#'          For \code{coxph} family, the response is a \code{\link[survival]{Surv}} 
#'          object, containing the survival time and censoring status.
#'          
#'          The terms contains the non-penalized predictor variables, and many groups 
#'          of prognostic and predictive terms, where in each group the prognostic 
#'          term comes first, followed by the predictive term.
#'          
#'          Different hierachical structures within groups can result from adjusting 
#'          the penalty parameters in the penalty function:          
#'          \deqn{\Omega(\bm{\beta}) = \lambda_1||\bm{\beta}|| +
#'          \lambda_2||\bm{\beta}||^2+\lambda_3|\beta_2|}
#'          Where \eqn{\bm{\beta}=(\beta_1,\beta_2)}. Note that \eqn{\beta_1} denotes
#'          the prognostic effect (main effect), and \eqn{\beta_2} for the predictive effect 
#'          (interactive effect), respectively. When \eqn{\lambda_2 = 0}
#'          and \eqn{\lambda_3 = 0}, it indicates no structure within groups.
#'          When \eqn{\lambda_2 \ne 0}, the penalty function honors the structure within groups
#'          such that: predictive effect \eqn{\ne 0 \Longrightarrow} prognostic effect \eqn{\ne 0}. 
#'          
#'          \code{rho,eabs,erel,LL,eta} are the corresponding parameters used in the itervative
#'          shrinkage-thresholding algorithm (ISTA) and the alternating direction method of 
#'          multipliers algorithm (ADMM). 
#'          
#'          Note that the missing values in the data are supposed to be dealt with in the 
#'          data preprocessing, before applying the method. 
#' @references \insertRef{ma2019structural}{smog}
#' 
#' @examples  
#' require(coxed)
#' 
#' n=50;p=100
#' set.seed(2018)
#' # generate design matrix x
#' s=10
#' x=matrix(0,n,1+2*p)
#' x[,1]=sample(c(0,1),n,replace = TRUE)
#' x[,seq(2,1+2*p,2)]=matrix(rnorm(n*p),n,p)
#' x[,seq(3,1+2*p,2)]=x[,seq(2,1+2*p,2)]*x[,1]
#' 
#' g=c(p+1,rep(1:p,rep(2,p)))  # groups 
#' v=c(0,rep(1,2*p))           # penalization status
#' label=c("t",rep(c("prog","pred"),p))  # type of predictor variables
#' 
#' # generate beta
#' beta=c(rnorm(13,0,2),rep(0,ncol(x)-13))
#' beta[c(2,4,7,9)]=0
#' 
#' # generate y
#' data1=x%*%beta
#' noise1=rnorm(n)
#' snr1=as.numeric(sqrt(var(data1)/(s*var(noise1))))
#' y1=data1+snr1*noise1
#' lfit1=smog(x,y1,g,v,label,lambda1=8,lambda2=0,lambda3=8,family = "gaussian")
#' 
#' ## generate binomial data
#' prob=exp(as.matrix(x)%*%as.matrix(beta))/(1+exp(as.matrix(x)%*%as.matrix(beta)))
#' y2=ifelse(prob<0.5,0,1)
#' lfit2=smog(x,y2,g,v,label,lambda1=0.03,lambda2=0,lambda3=0.03,family = "binomial")
#' 
#' ## generate survival data
#' data3=sim.survdata(N=n,T=100,X=x,beta=beta)
#' y3=data3$data[,c("y","failed")]
#' y3$failed=ifelse(y3$failed,1,0)
#' colnames(y3)=c("time","status")
#' lfit3=smog(x,y3,g,v,label,lambda1=0.2,lambda2=0,lambda3=0.2,family = "coxph")
#' 
#' @export
smog.default <- function(x, y, g, v, label, lambda1, lambda2, lambda3, family = "gaussian", subset = NULL, rho = 10, 
                         scale = TRUE, eabs = 1e-3, erel = 1e-3, LL = 1, eta = 1.25, maxitr = 1000, ...){
  
  lambda=c(lambda1,lambda2,lambda3)
  hierarchy = ifelse(lambda2 == 0, 
                     ifelse(lambda3 == 0, 0, 1),2)
  
  if(!is.null(subset)){
    x <- as.matrix(as.matrix(x)[subset,])
    y <- as.matrix(as.matrix(y)[subset,])
  }else{
    x <- as.matrix(x)
    y <- as.matrix(y)
  }

  g <- as.numeric(as.factor(g))
  est <- glog(y,x,g,v,lambda,hierarchy,family,rho,
              scale,eabs,erel,LL,eta,maxitr)
  
  if(nrow(est$coefficients)){  # continue for some variables are selected  
    if(family == "gaussian"){
      wx <- cbind(rep(1,nrow(x)),x)
      est$fitted.value = as.vector(as.matrix(wx[,est$coefficients$Id+1])%*%
                                     as.matrix(est$coefficients$Estimate))
      est$residuals = as.vector(y - est$fitted.value)
      
      if(!is.null(colnames(x))){
        est$coefficients$Beta = c("Intercept",colnames(x))[est$coefficients$Id+1]
        est$coefficients = est$coefficients[,c("Id","Beta","Estimate")]
      }
    }
    
    if(family == "binomial"){
      probTAB = exp(as.matrix(x[,est$coefficients$Id])%*%as.matrix(est$coefficients$Estimate))/
        (1+exp(as.matrix(x[,est$coefficients$Id])%*%as.matrix(est$coefficients$Estimate)))
      probTAB = cbind(1-rowSums(as.matrix(probTAB)),probTAB)
      predClass = apply(probTAB,1,which.max)
      predProb = apply(probTAB,1,max)
      
      est$levels = sort(unique(y))
      est$fitted.value = data.frame(Class = est$levels[predClass],
                                    Prob = predProb)
      
      if(!is.null(colnames(x))){
        est$coefficients$Beta = colnames(x)[est$coefficients$Id]
        est$coefficients = est$coefficients[,c(1,ncol(probTAB)+1,2:ncol(probTAB))]
      }
    }
    
    if(family == "coxph"){
      est$fitted.value = as.vector(exp(-exp(as.matrix(x[,est$coefficients$Id])%*%
                                            as.matrix(est$coefficients$Estimate))))
      
      if(!is.null(colnames(x))){
        est$coefficients$Beta = colnames(x)[est$coefficients$Id]
        est$coefficients = est$coefficients[,c("Id","Beta","Estimate")]
      }
    }
  }
  
  # calculate the degrees of freedom 
  idx=est$coefficients$Id
  est$model$intercept=ifelse(0 %in% idx,
                             est$coefficients$Estimate[idx==0], NA)
  est$model$treatment=ifelse(1 %in% idx,
                             est$coefficients$Estimate[idx==1], NA)
  GId=estimate=elabel=NULL
  if(any(idx>1)){
    est$model$biomarker=data.frame(GId=g[idx[idx>1]],
                                   estimate=est$coefficients$Estimate[idx>1],
                                   elabel=label[idx[idx>1]])
    est$model$biomarker=tidyr::spread(est$model$biomarker,elabel,estimate)
    est$model$biomarker=as.data.frame(apply(est$model$biomarker,1:2,
                                            function(t) ifelse(is.na(t),0,t)))

    if("prog" %in% colnames(est$model$biomarker) & 
       (!("pred" %in% colnames(est$model$biomarker)))){
      est$model$biomarker$pred=0
    }
    
    if("pred" %in% colnames(est$model$biomarker) & 
       (!("prog" %in% colnames(est$model$biomarker)))){
      est$model$biomarker$prog=0
    }
    
    est$weight=data.frame(t(apply(est$model$biomarker[,c("prog","pred")],1,
                                  function(t) lambda[1]/sqrt(sum(t^2))+
                                    c(0,ifelse(sign(t[2]),lambda[3]/abs(t[2]),0)))))
    colnames(est$weight)=c("prog","pred")
    
    if(!is.null(colnames(x))){
      est$model$biomarker$marker=colnames(x)[label=="prog" & (1:ncol(x) %in% idx[idx>1])]
      est$model$biomarker=est$model$biomarker[,c("GId","marker","prog","pred")]
      est$weight=cbind(est$model$biomarker[,c("GId","marker")],est$weight)
    }else{
      est$model$biomarker=est$model$biomarker[,c("GId","prog","pred")]
      est$weight=cbind(GId=est$model$biomarker$GId,est$weight)
    }
    
    Weight=sapply(idx[idx>=1],function(t) ifelse(t==1,0,
                                                 est$weight[est$weight$GId==g[t],label[t]]))
    if(family == "gaussian"){
      est$DF=sum(diag(solve(as.matrix(t(x[,idx[idx>=1]]))%*%as.matrix(x[,idx[idx>=1]])+diag(Weight))
                      %*%as.matrix(t(x[,idx[idx>=1]]))%*%as.matrix(x[,idx[idx>=1]])))
    }else{
      est$DF = sum(1/(1+Weight))
    }
    
  }else{
    if(!is.null(colnames(x))){
      est$model$biomarker=data.frame(GId=NA,marker=NA,prog=NA,pred=NA)
    }else{
      est$model$biomarker=data.frame(GId=NA,prog=NA,pred=NA)
    }
    
    est$weight=NULL
    est$DF=1
  }
  
  # AIC, BIC, GCV criteria by using whole data 
  n <- nrow(x)
  k = est$DF
  
  aic1 = tryCatch({n/2*log(abs(2*est$llikelihood)) + n/2*((1+k/n)/(1-k+2/n))},
                  error=function(e) NA)
  aic2 = tryCatch({log(abs(2*est$llikelihood)/n) + 2*k/n}, error=function(e) NA)
  bic = tryCatch({log(abs(2*est$llikelihood)/n) + log(n)*k/n}, error=function(e) NA)
  gcv = tryCatch({abs(2*est$llikelihood)/(n*(1-k/n)^2)}, error=function(e) NA)

 
  est$criteria=list(aic1,aic2,bic,gcv)
  names(est$criteria)=c("cAIC","AIC","BIC","GCV")
  
  est$call <- match.call()
  class(est) <- "smog"
  est
}

#' @rdname smog.default
#' 
#' @seealso \code{\link{cv.smog}}, \code{\link{predict.smog}}, \code{\link{plot.smog}}.
#' @author Chong Ma, \email{chongma8903@@gmail.com}.
#'  
#' @export
smog.formula <- function(formula, data=list(), g, v, label, lambda1, lambda2, lambda3, ...){
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf,"terms"),data = mf)
  y <- model.response(mf)
  
  est <- smog.default(x,y,g,v,label,lambda1,lambda2,lambda3,...)
  est$call <- match.call()
  est$formula <- formula
  est
}


#' Cross-valiation for smog 
#' 
#' \code{cv.smog} conducts a greedy-search for optimal lambda's and yields a sparse
#' model given a provided model selection criterion. When type is ``nloglike'', 
#' the method allows the \code{nfolds} to be processed in parallel for speeding up 
#' the cross-validation.
#' 
#' @inheritParams smog.default
#' @param type model selction criterion, should choose from ``nloglike'', 
#'             ``cAIC'', ``AIC'', ``BIC'', and ``GCV'', respectively. 
#' @param lambda.max the maximum value for lambda's. If \code{NULL}, the default \code{lambda.max}
#'                   is \eqn{1/\lambda_{min}(x'x)}.
#' @param nlambda.max the maximum number of lambdas' shrunk down from the maximum lambda \code{lambda.max}.
#'                    Default is 20.
#' @param delta the damping rate for lambda's such that \eqn{\lambda_k = \delta^k\lambda_0}. Default is 0.9.
#' @param nfolds number of folds. One fold of the observations in the data are used
#'               as the testing, and the remaining are fitted for model training. 
#'               Default is 5.
#' @param parallel Whether or not process the \code{nfolds} cross-validations in
#'                 parallel. If \code{TRUE}, use \code{\link[foreach]{foreach}} to do each 
#'                 cross-validation in parallel. Default is \code{FALSE}.
#' @param ncores number of cpu's for parallel computing. See
#'               \code{\link[parallel]{makeCluster}} and \code{\link[doParallel]{registerDoParallel}}.
#'               Default is \code{NULL}. 
#' @param ... other arguments that can be supplied to \code{smog}.
#'
#' @return includes the profile containing a path of lambda's and the corresponding model 
#'         selectio criterion value, the optimal lambda's, and the optimal model, respectively.
#'         The \code{type} comes from a list of model selection criteria values, includes the 
#'         average of the negative log-likelihood values and the correction AIC for each fold of the data.
#'         
#'         \item{cvfit}{the fitted model based on the optimal lambda's.}
#'         \item{lhat}{the optimal lambda's which has the minimum model selection criterion.}
#'         \item{profile}{a data frame contains the path of lambda's and the corresponding model selection
#'                        criterion, which is determined by the \code{type}.}
#'         
#' @details When the \code{type} is ``nloglike'', it requires doing \code{nfolds} cross-validations. 
#'          For each cross-validation, evenly split the whole data into \code{nfolds}, and one fold of 
#'          the observations are used as the testing data, and the remaining are used for model training. 
#'          After calculating the AIC for each fold of testing data, return the average of the 
#'          AICs. Note that we keep \code{lambda2}\eqn{=0} during the greedy search for lambda's. 
#'          
#' @examples
#' 
#' # generate design matrix x
#' set.seed(2018)
#' n=50;p=100
#' s=10
#' x=matrix(0,n,1+2*p)
#' x[,1]=sample(c(0,1),n,replace = TRUE)
#' x[,seq(2,1+2*p,2)]=matrix(rnorm(n*p),n,p)
#' x[,seq(3,1+2*p,2)]=x[,seq(2,1+2*p,2)]*x[,1]
#' 
#' g=c(p+1,rep(1:p,rep(2,p)))  # groups 
#' v=c(0,rep(1,2*p))           # penalization status
#' label=c("t",rep(c("prog","pred"),p))  # type of predictor variables
#' 
#' # generate beta
#' beta=c(rnorm(13,0,2),rep(0,ncol(x)-13))
#' beta[c(2,4,7,9)]=0
#' 
#' # generate y
#' data=x%*%beta
#' noise=rnorm(n)
#' snr=as.numeric(sqrt(var(data)/(s*var(noise))))
#' y=data+snr*noise
#' 
#' cvfit=cv.smog(x,y,g,v,label,type="GCV",family="gaussian")
#' 
#' @seealso \code{\link{smog.default}}, \code{\link{smog.formula}}, \code{\link{predict.smog}}, \code{\link{plot.smog}}.
#' @author Chong Ma, \email{chongma8903@@gmail.com}.
#' 
#' @references \insertRef{ma2019structural}{smog}
#' 
#' @export
cv.smog <- function(x, y, g, v, label, type, family = "gaussian", lambda.max = NULL, 
                    nlambda.max = 20, delta = 0.9, nfolds = 5, parallel = FALSE, ncores = NULL, ...){
  # define the maximum lambda value 
  typelist=c("nloglike","cAIC","AIC","BIC","GCV")
  n = nrow(x)
  p = ncol(x)
  
  if(!(type %in% typelist)){
    stop(paste(c(type,"is not one of",typelist,"!"),collapse = " "))
  }
  
  if(family == "gaussian"){
    wx <- cbind(rep(1,nrow(x)),x)
  }
  
  if(type %in% c("nloglike")){
    setlist = as.integer(seq(0,nrow(x),length.out = nfolds+1))
    ncores = ifelse(parallel,ncores,1)
    
    if(parallel){
      `%mypar%` =  `%dopar%`
      ncores = ifelse(is.null(ncores),1,ncores)
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
    }else{
      `%mypar%` = `%do%`
    }
  }
  
  if(family == "gaussian"){
    lambda.max = ifelse(is.null(lambda.max),10*max(t(x)%*%y)/n*max(log(max(p,n)/n),1),lambda.max)
  }
  
  if(family == "binomial"){
    lambda.max = ifelse(is.null(lambda.max),1/min(svd(x)$d)*max(log(max(p,n)/n),1),lambda.max) 
  }
  
  if(family == "coxph"){
    lambda.max = ifelse(is.null(lambda.max),1/min(svd(x)$d)*p/n,lambda.max) 
  }
  
  temp.fit0 = smog::smog(x,y,g,v,label,lambda1=lambda.max,lambda2=0,lambda3=lambda.max,family,...)
  
  if(nrow(coef(temp.fit0))<4){
    repeat{
      lambda.max = lambda.max*delta
      temp.fit1=smog::smog(x,y,g,v,label,lambda1=lambda.max,lambda2=0,lambda3=lambda.max,family,...)
      if(nrow(coef(temp.fit1))>=4)  break
    }
  }else{
    repeat{
      lambda.max = lambda.max/delta
      temp.fit1=smog::smog(x,y,g,v,label,lambda1=lambda.max,lambda2=0,lambda3=lambda.max,family,...)
      if(nrow(coef(temp.fit1))<=4)  break
    }
  }
  
  l1 = l2 = lambda.max
  profile = list(NULL)
  lhat = cvfit = NULL
  
  it = 1    # initialize the search step
  temp=rep(NA,3)
  repeat{
    temp.lambda=list(c(l1*delta,0,l2),c(l1,0,l2*delta),c(l1*delta,0,l2*delta))
    
    if(!(type %in% c("nloglike"))){
      # non-cross-validation type model selection
      for(k in 1:3){
        temp[k] = smog::smog(x,y,g,v,label,lambda1=temp.lambda[[k]][1],lambda2=temp.lambda[[k]][2],lambda3=temp.lambda[[k]][3],family,...)$criteria[[type]]
      }
    }else{
      # cross-validation to calculate the maximum log-likelihood 
      for(k in 1:3){
        i = 1
        cv_res <- foreach::foreach(i=1:nfolds,.combine = c,
                                   .packages = c("foreach","smog"))%mypar%{
                                     tset <- (setlist[i]+1):setlist[i+1]
                                     lset <- c(1:nrow(x))[-tset]
                                     cvmodel <- smog::smog(x,y,g,v,label,lambda1=temp.lambda[[k]][1],lambda2=temp.lambda[[k]][2],lambda3=temp.lambda[[k]][3],family,subset=lset,...)
                                     nt <- length(tset)
                                     kt = cvmodel$DF
                                     
                                     if(family == "gaussian"){
                                       tx = wx[tset,coef(cvmodel)$Id+1]
                                       ty = as.vector(as.matrix(tx) %*% as.matrix(coef(cvmodel)$Estimate))
                                       tloglike = -1/2*sum((y[tset] - ty)^2)
                                     }
                                     
                                     if(family == "binomial"){
                                       probTAB = exp(as.matrix(x[tset,coef(cvmodel)$Id])%*%as.matrix(coef(cvmodel)$Estimate))/
                                         (1+exp(as.matrix(x[tset,coef(cvmodel)$Id])%*%as.matrix(coef(cvmodel)$Estimate)))
                                       probTAB = cbind(1-rowSums(as.matrix(probTAB)),probTAB)
                                       ty = match(y[tset],cvmodel$levels)
                                       tloglike = sum(log(apply(cbind(ty,probTAB),1,function(t) t[-1][t[1]]))) 
                                     }
                                     
                                     if(family == "coxph"){
                                       tx=x[tset,coef(cvmodel)$Id]
                                       ty=y[tset,]
                                       
                                       ttheta = exp(as.matrix(tx)%*%as.matrix(coef(cvmodel)$Estimate))
                                       tloglike = as.numeric(colSums(as.matrix(tx))%*%as.matrix(coef(cvmodel)$Estimate))
                                       
                                       t=NULL
                                       for(t in ty$time[ty$status>0]){
                                         tid0=which(ty$time == t & ty$status>0)
                                         tid1=which(ty$time >= t)
                                         ti = NULL
                                         for(ti in 0:(length(tid0)-1)){
                                           tloglike=tloglike-log(sum(ttheta[tid1])-ti/length(tid0)*sum(ttheta[tid0]))
                                         }
                                       }
                                      }
                                    
                                     c(-tloglike)
                                   }
        temp[k] = mean(cv_res,na.rm = TRUE)
      }
    }
    
    if(1 %in% which(temp==min(temp)) | 3 %in% which(temp==min(temp)) 
       | min(temp[c(1,3)]) - temp[2] < 0.001*abs(min(temp[c(1,3)]))
       ){
      profile[[it]] = c(temp.lambda[[3]][c(1,3)],temp[3])
    }else{
      profile[[it]] = c(temp.lambda[[2]][c(1,3)],temp[2])
    }
    
    l1=profile[[it]][1]
    l2=profile[[it]][2]
    
    # determine the stopping criterion
    if(it>=5){
        if((abs(profile[[it-1]][3] - profile[[it]][3]) > 0.01*abs(profile[[it-1]][3]) & 
            abs(profile[[it-2]][3] - profile[[it-1]][3]) > 0.01*abs(profile[[it-2]][3]) &
            (profile[[it-2]][3] - profile[[it-1]][3] > 10*(profile[[it-1]][3] - profile[[it]][3]) |
            (profile[[it-2]][3] - profile[[it-1]][3] > 0.01*abs(profile[[it-2]][3]) & 
             profile[[it-1]][3] - profile[[it]][3] > 10*(profile[[it-2]][3] - profile[[it-1]][3])))
            ) | it >= nlambda.max){
          break
        } 
    }
    it = it + 1
  }
  if(parallel)  parallel::stopCluster(cl)
  
  profile = data.frame(do.call(rbind,profile))
  colnames(profile)=c("lambda1","lambda2",type)
  lhat = as.numeric(profile[which.min(profile[,3]),1:2])
  cvfit = smog(x,y,g,v,label,lambda1=lhat[1],lambda2=0,lambda3=lhat[2],family,...)$model
  
  return(list(cvfit=cvfit,lhat=lhat,profile=profile))
}


#' predict method for objects of the class smog
#' 
#' \code{predict.smog} can produce the prediction for user-given new data, based on the
#' provided fitted model (\code{object}) in the S3method of \code{smog}. If the \code{newdata} omitted,
#' it would output the prediction for the fitted model itself. The yielded result should
#' match with the family in the provided model. See \code{\link{smog}}.
#' 
#' @param object a fitted object of class inheriting from smog.
#' @param newdata a data frame containing the predictor variables, which are
#'                used to predict. If omitted, the fitted linear predictors 
#'                are used. 
#' @param family  a description of distribution family for which the response 
#'                variable is to be predicted.  
#' @param ... additional arguments affecting the predictions produced.
#' 
#' @details If \code{newdata = NULL}, the fitted.value based on the \code{object}
#'          is used for the prediction. 
#' 
#' @return If \code{family} = ``gaussian'', a vector of prediction for the response is returned.
#'         For \code{family} = ``coxph'', a vector of predicted survival 
#'         probability is returned. When \code{family} = ``binomial'', it outputs a data
#'         frame containing the predicted group labels and the corresponding 
#'         probabilies. 
#' 
#' @seealso \code{\link{smog.default}}, \code{\link{smog.formula}}, \code{\link{cv.smog}}, \code{\link{plot.smog}}.
#' @author Chong Ma, \email{chongma8903@@gmail.com}.
#' 
#' @references \insertRef{ma2019structural}{smog}
#' 
#' @export
predict.smog <- function(object, newdata = NULL, family = "gaussian",...){
  if(is.null(newdata)){
    y <- fitted(object)
  }else{
    if(!is.null(object$formula)){
      x = model.matrix(object$formula, newdata)
    }else{
      x = newdata
    }
    
    if(nrow(coef(object))){
      if(family == "gaussian"){
        wx = cbind(rep(1,nrow(x)),x)
        y <- as.vector(as.matrix(wx[,coef(object)$Id+1]) %*% as.matrix(coef(object)$Estimate))
      }
      
      if(family == "binomial"){
        probTAB = exp(as.matrix(x[,coef(object)$Id])%*%as.matrix(coef(object)$Estimate))/
          (1+exp(as.matrix(x[,coef(object)$Id])%*%as.matrix(coef(object)$Estimate)))
        probTAB = cbind(1-rowSums(probTAB),probTAB)
        predClass = apply(probTAB,1,which.max)
        predProb = apply(probTAB,1,max)
        
        y = data.frame(Class = object$levels[predClass],
                       Prob = predProb)
      }
      
      if(family == "coxph"){
        y = round(exp(-exp(as.matrix(x[,coef(object)$Id])%*%as.matrix(coef(object)$Estimate))),2)
      }
      
      return(y);
    }else{
      stop("the model is not appropriate for making predictions")
    }
  }
}


#' plot method for objects of the class smog
#' 
#' \code{plot.smog} can produce a panel of plots for the primal errors, dual errors, 
#' and the penalized log-likelihood values, based on the provided fitted model 
#' (\code{x}) in the S3method of \code{smog}. 
#' 
#' @param x a fitted object of class inheriting from smog.
#' @param type,xlab default line types and x axis labels for the panel of plots.
#' @param caption a list of y axes labels for the panel of plots. 
#' @param ... additional arguments that could be supplied to \code{\link[graphics]{plot}}
#'            and \code{\link[graphics]{par}}.
#' 
#' @details For the panel of three plots, the \code{xlab} is ``iterations'' and the
#'          \code{type} is ``l'', by default. The \code{ylab} are ``primal error'',
#'          ``dual error'',``log-likelihood'', respectively. This panel of plots can
#'          reflect the convergence performance for the algorithm used in \code{\link{smog}}.
#' 
#' @seealso \code{\link[graphics]{par}}, \code{\link[graphics]{plot.default}}, \code{\link{predict.smog}},
#'          \code{\link{smog.default}}, \code{\link{smog.formula}}, \code{\link{cv.smog}}.
#'          
#' @author Chong Ma, \email{chongma8903@@gmail.com}.
#' 
#' @references \insertRef{ma2019structural}{smog}
#' 
#' @export
plot.smog <- function(x,type = "l",xlab="iteration",caption=list("primal error","dual error","log-likelihood"),...){
  op<-graphics::par(mfrow=c(2,2),mar=c(2.1,4.1,2.1,2.1),no.readonly = TRUE,...)
  graphics::plot(x$PrimalError, type=type, xlab = xlab, ylab = caption[[1]],...)
  graphics::plot(x$DualError, type=type, xlab = xlab, ylab = caption[[2]],...)
  graphics::plot(x$loglike, type=type, xlab = xlab, ylab = caption[[3]],...)
  on.exit(par(op))
}


#' smog generic
#'
#' @param x an object of class from``smog''.
#' @param ... further arguments passed to or from other methods.
#' @keywords internal
#' 
#' @seealso \code{\link{cv.smog}}, \code{\link{smog.default}}, \code{\link{smog.formula}}, 
#'          \code{\link{predict.smog}}, \code{\link{plot.smog}}.
#' @return the cofficients table of the object x of class from ``smog''.
#'         See \code{\link[base]{print}}.
#' 
#' @references \insertRef{ma2019structural}{smog}
#' 
smog <- function(x, ...) UseMethod("smog")

#' @rdname smog
#' @keywords internal
#' 
print.smog <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  
  cat("\n Coefficients:\n")
  print(x$coefficients,row.names=FALSE)
}








