#naive bayes for numeric features
#################################

# useful functions

####### dsnorm ###########
# skewed normal probality density function
#______________________________________________________________________
dsnorm <- function(x,mu=0,sigma=1,shape=1)
{
  x <- (x-mu)/sigma
  2*dnorm(x)*pnorm(shape*x)/sigma
}

######## dsnorm_params ###############
#______________________________________________________________________
dsnorm_params <- function(x,type = "vector")
{
  N <- length(x)
  # second central moment
  m2 <- sum((x-mean(x))**2)/N
  # third central moment
  m3 <- sum((x-mean(x))**3)/N
  # skewness
  gamma <- m3/sqrt(m2**3)
  absgamma <- abs(gamma)
  
  # delta
  pm <- pi/2
  e <- 2/3
  delta <- sqrt(pm*(absgamma**e)/(absgamma**e+((4-pi)/2)**e))
  delta <- ifelse(delta >= 0,min(0.999,delta),max(-0.999,delta))
  # shape
  shape <- sqrt(delta/sqrt(1-delta**2))
  if (gamma < 0) shape <- -shape
  # sigma
  sigma <- sqrt(var(x)*(1-(2*delta**2/pi)))
  # mu
  if (gamma < 0)
    mu <- mean(x) + sigma*delta*sqrt(2/pi)
  else
    mu <- mean(x) - sigma*delta*sqrt(2/pi)
  # return
  if (type == "list")
    list(mu=mu,sigma=sigma,shape=shape)
  else if (type == "vector")
    c(mu=mu,sigma=sigma,shape=shape)
}

########## prior_prob ##################
#______________________________________________________________________
prior_prob <- function(train,test,class,k)
{
  class <- as.factor(class)
  lev <- levels(class)
  # initialize pp list
  pp <- list()
  for (l in lev)
  {
    pp <- c(pp,list(c()))
  }
  names(pp) <- lev
  # knn
  for (rt in 1:nrow(test))
  {
    dist_rt <- 0
    # euclidian distance to every train point, for each test entry
    for (n in names(test))
    {
      dist_rt <- dist_rt +(test[rt,n]-train[,n])**2
    }
    dist_rt <- sqrt(dist_rt)
    # get first k nearest neighbours
    first_k <- order(dist_rt)[1:k]
    # prior probabilities refered to the k-nn + 1 (laplace estimator)
    tab <- prop.table(table(class[first_k])+1)
    for (l in lev)
      pp[[l]] <- c(pp[[l]],tab[[l]])
  }
  pp
}

####### build model ##################
#______________________________________________________________________
nbknn <- function(train,test,class,k)
{
  require(stats)
  # sanity cheks
  # check if train is a data frame
  if (! is.data.frame(train))
  {
    train <- as.data.frame(train)
    warning("train coerced to data frame")
  }
  # check if nrow train = nclass
  if (nrow(train) != length(class))
    stop("the numbers of train and class observation doesn't macht")
  # check if features are numeric
  for (i in ncol(train))
  {
    if (! is.numeric(train[,i]))
      stop(paste("feature",train[,i],"isn't numeric"))
  }
  # check if class is a factor
  if (! is.factor(class))
    stop("class must be a factor")
  # get names of features
  fea_names <- names(train)
  # get class levels
  clevels <- levels(class)
  # prior probability for each level of the class
  pp <- prior_prob(train,test,class,k)
  # pdf's parameters for each level of the class
  # for each level, loop over the features
  pdf_par <- list()
  for (l in clevels)
  {
    params <- list()
    for (f in fea_names)
    {
      data <- train[,f][class == l]
      # get pdf parameters for likelihood computation
      params <- c(params,list(c(mu=mean(data),sigma=sd(data))))
    }
    names(params) <- fea_names
    # store pdf params in a list named with levels name
    pdf_par <- c(pdf_par,list(params))
  }
  names(pdf_par) <- clevels
  # output the model object of class nbnum
  # a list with:
  # train data feature names
  # class levels names
  # all parameters for each level
  # prior probability for each level
  model <- structure(list(fnames=fea_names,
                       clevels=clevels,
                       pdfPar=pdf_par,
                       priProbs=pp), class=c("nbknn","list"))
  model
}

########### predict method for nbknn class ################
#______________________________________________________________________
predict.nbknn <- function(model,test,type="class",laplace=0)
{
  # sanity cheks
  if (laplace != 0)
    laplace <- laplace/100
  # if test is not a data frame coerce it
  if (! is.data.frame(test))
  {
    test <- as.data.frame(test)
    warning("test data coerced to data frame")
  }
  # check if names in test are equal to the names in model
  if (! identical(model$fnames,names(test)))
    stop("test and training features differ")
  # check if test's features are numeric
  for (i in ncol(test))
  {
    if (! is.numeric(test[,i]))
      stop("test features must be numeric")
  }
  # for each level
  post_prob <- list()
  scale_fact <- 0  # scaling factor
  for (l in model$clevels)
  {
    # calculate likelihood
    lparams <- model$pdfPar[[l]]
    li <- 1
    for (f in model$fnames)
    {
      # p(Fi|level) for Fi features
      # multiply all p (likelihoods)
      pdf_par <- lparams[[f]]
      #li <- li * (dsnorm(test[,f],pdf_par["mu"],pdf_par["sigma"],pdf_par["shape"])+laplace)
      li <- li * (dnorm(test[,f],pdf_par["mu"],pdf_par["sigma"])+laplace)
    }
    # posterior probability p(level|Fi)
    # likelihood X prior probability
    postP <- li * model$priProbs[[l]]
    scale_fact <- scale_fact + postP
    post_prob <- c(post_prob,list(postP))
  }
  names(post_prob) <- model$clevels
  # normalized probability
  post_prob <- lapply(post_prob,function(x) x/scale_fact)
  #post_prob <- as.data.frame(post_prob)
  # calculate predicted levels
  pred <- c()
  for (i in 1:nrow(test))
  {
    maxp <- post_prob[[1]][i]
    maxl <- names(post_prob)[1]
    for (l in model$clevels)
    {
      if (post_prob[[l]][i] > maxp)
      {
        maxp <- post_prob[[l]][i]
        maxl <- l
      }
    }
    pred <- c(pred,maxl)
  }
  # if type==class return level with max post. prob.
  if (type == "class")
  {
    return(pred)
  }
  else if (type == "raw")
  {
    # if type==raw return probs
    post_prob <- as.data.frame(post_prob)
    names(post_prob) <- paste("prob_",model$clevels,sep="")
    post_prob$class <- pred
    return(post_prob)
  }
  else
    stop(paste("unrecognized type '",type,"'",sep=""))
}
