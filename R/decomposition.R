#' This function calculates the explained or endowment part of a Kitawaga-Oaxaca-Blinder (KOB) decomposition
#' generalized to a continuous grouping variable. Rarely used directly by the user.
#'
#' @param m a correctly specified model for an KOB decomposition
#' @param xes a character vector giving the names of the independent (or control) variables
#' @param ges a character vector giving the name(s) of the grouping variable and any polynomial transformations. The list should be ordered by increasing polynomial.
#' @param {i, j} values of the group variable to compare
#' @return The amount of difference attributable to the explained part in a KOB decomposition
#' @export
getExplained<-function(m,xes,ges,i=1,j=0) {
  totalE<-0
  for(x in xes) {
    a<-estimateX(m$model[[x]],m$model[[ges[1]]],i)
    b<-m$coefficients[x]
    for( n in c(1:length(ges)) ) {
      g<-ges[n]
      b<-b+(m$coefficients[paste0(x,":",g)]*i^n)
    }
    totalE<-totalE + a*b
  }

  totalE2<-0
  for(x in xes) {
    a<-estimateX(m$model[[x]],m$model[[ges[1]]],j)
    b<-m$coefficients[x]
    for( n in c(1:length(ges)) ) {
      g<-ges[n]
      b<-b+(m$coefficients[paste0(x,":",g)]*i^n)
    }
    totalE2<-totalE2 + a*b
  }
  return(totalE-totalE2)
}

#' This function calculates the unexplained part of a Kitawaga-Oaxaca-Blinder (KOB) decomposition
#' generalized to a continuous grouping variable. This function is rarely directly used by the user.
#'
#' @param m a correctly specified model for an KOB decomposition
#' @param xes a character vector giving the names of the independent (or control) variables
#' @param ges a character vector giving the name(s) of the grouping variable and any polynomial transformations. It is assumed that the first element is the original grouping variable.
#' @param {i, j} values of the group variable to compare
#' @return The amount of difference attributable to the unexplained part of a KOB decomposition.
#' @export
getUnexplained<-function(m,xes,ges,i=1,j=0) {
  a<-m$coefficients[ges[1]]*(i-j)
  for(x in xes) {
    b<-estimateX(m$model[[x]],m$model[[ges[1]]],j)
    c<-m$coefficients[x]
    d<-m$coefficients[x]
    for( n in c(1:length(ges)) ) {
      g<-ges[n]
      c<-c+(m$coefficients[paste0(x,":",g)]*i^n)
      d<-d+(m$coefficients[paste0(x,":",g)]*j^n)
    }

    a<-a + b*(c-d)
  }
  return(a)
}

#' Convenience function
#'
#' Prints the summary of the passed model, as well as the
#' percent explained and percent unexplained

decompose<-function(m,xes,ges,i=1,j=0) {
  ex<-unname(getExplained(m,xes=xes,ges=ges,i=i,j=j))
  un<-unname(getUnexplained(m,xes=xes,ges=ges,i=i,j=j))
  yvar<-all.vars(m$terms)[1]
  total<-estimateX(m$model[[yvar]],m$model[[ges[1]]],i)-estimateX(m$model[[yvar]],m$model[[ges[1]]],j)
  de.list<-list()
  de.list[["model"]]<-m
  de.list[["total"]]<-total
  de.list[["explained"]]<-ex
  de.list[["unexplained"]]<-un
  return(de.list)
}

#' print summary decompose
#'

prDecomp<-function(de) {
  #print(summary(de[["model"]]))
  print(paste0("Total difference: ",round(de[["total"]],digits=4)))
  print(paste0("Explained difference: ",round(de[["explained"]],digits=4)))
  print(paste0("Explained percent: ",round(100*de[["explained"]]/de[["total"]],digits=4)))
  print(paste0("Unexplained difference: ",round(de[["unexplained"]],digits=4)))
  print(paste0("Unexplained percent: ",round(100*de[["unexplained"]]/de[["total"]],digits=4)))
}

#' Generate simulated residuals for a latent variable
#'
#' Based on Wolff (2012), this will generate simulated residuals from a probit model
#' such that the generated residual plus the linear prediction corresponds to the observed
#' outcome (i.e. Y==1:Y*>0; Y==0:Y*<=0)

simulateResiduals<-function(m,df,outcome="y") {
  df$trueres<-NA
  df$works<-FALSE
  df$yhat<-predict(m,type="link")
  n<-1
  while( !all(df$works) ) {
    df$res<-ifelse(df$works,df$trueres,rnorm(dim(df)[1])*n)
    df$works<-(ifelse(df$yhat+df$res<=0,0,1)==df[[outcome]])
    df$trueres<-ifelse(df$works==TRUE&is.na(df$trueres),df$res,df$trueres)
    n<-n+0.1
  }
  return(df$trueres)
}

#' Generate simulated residuals for a latent variable for mlogit
#'
#' Based on Wolff (2012), this will generate simulated residuals from a mlogit model
#' such that the generated residual plus the linear prediction corresponds to the observed
#' outcome

simMlogitRes<-function(m,df,outcome="y") {
  #assumes m is mlogit
  fit<-fitted(m,outcome=FALSE)
  fit<-apply(fit,c(1,2),logit)
  num.cols<-dim(fit)[2]
  num.rows<-dim(fit)[1]
  trueres<-matrix(NA,nrow=dim(fit)[1],ncol=dim(fit)[2])
  truey<-df[[outcome]]
  ynames<-colnames(fit)
  works<-rep(FALSE,length.out=dim(df)[1])
  n<-1
  while( !all(works) ) {
    res<-matrix(n*rlogis(num.cols*num.rows),nrow=num.rows,ncol=num.cols)
    ystar<-fit+res
    predy<-apply(ystar,1,function(x){which(x==max(x))})
    works<-predy==truey
    for( i in c(1:length(predy)) ) {
      if( (ynames[predy[i]] == truey[i]) & all(is.na(trueres[i,])) ) {
        trueres[i,]<-res[i,]
      }
    }
    works<-apply(trueres,1,function(x) all(!is.na(x)))
    n<-n+0.1
  }
  return(trueres)
}

logit.reverse<-function(x) {
  return( 1/(1+exp(-1*x)))
}

logit<-function(x) {
  log(x/(1-x))
}

#' This function calculates the average value of x at values of g=i
estimateX<-function(x,g,i=0) {
  a<-x[which(g==i)]
  b<-mean(a,na.rm = TRUE)
  return(b)
}

#' Given a properly expressed probit regression, decompose the results into an
#' explained and unexplained parts.
#'
#' @param fitted.mdl A probit regression object
#' @param list_of_IVs A character vector indicator the predictor variables included in the probit model
#' @param list_of_grouping_vars A character vector of grouping variables included in the probit model
#' @param values_of_grouping_vars A vector containing all possible values of the grouping variable.
#' @return A list containing two objects. One, the fitted OLS model based
#'   on the estimated log-odds continuous variable. Two, a data frame containing the decomposition
#'   comparing every possible combination of values for the grouping variable.
#'   @export
gen_decomposed_results_probit <- function(fitted.mdl,
                                          list_of_IVs,
                                          list_of_grouping_vars,
                                          values_of_grouping_vars) {
  # save off the data frame used to fit the model
  df.an <- fitted.mdl$data

  # get name of y varaible
  yname <- all.vars(fitted.mdl$formula)[1]

  # -- Calculate continuous log-odds DV --
  # step 1: generate simulated residuals
  df.an$simres <- simulateResiduals(fitted.mdl, df.an, yname)

  # step 2: calculate ystar - the predicted log-odds plus the simulated residuals
  df.an$ystar <- predict(fitted.mdl, type = "link") + df.an$simres

  # generate new formula to use in OLS regression
  char.x <- attr(terms(fitted.mdl$formula), which = "term.labels")
  char.x <- paste0(char.x, collapse = " + ")
  char.x <- paste0("ystar ~ ", char.x, collapse = "")
  new.form <- formula(char.x)

  # run an OLS regression using the new ystar as the response variable
  fitted.lm <- lm(new.form, data = df.an)

  # now generate a data frame containing the decomposition of the OLS regression
  decomp.df <- data.frame()
  for( y1 in values_of_grouping_vars ) {
    for( y2 in values_of_grouping_vars ) {
      if( y1 != y2 ) {
        decom <- decompose(fitted.lm,
                           xes = list_of_IVs,
                           ges = list_of_grouping_vars,
                           i = y2,
                           j = y1)
        df <- data.frame(
          value1=y1,
          value2=y2,
          explained=decom[["explained"]],
          unexplained=decom[["unexplained"]])
        decomp.df <- decomp.df %>% bind_rows(df)
      }
    }
  }

  decomp.df<-decomp.df %>% mutate(total=explained+unexplained,
                                  percent.explained=100*explained/total,
                                  percent.unexplained=100*unexplained/total)

  # return fitted OLS and decomp results data frame
  list(
    "fitted.ols" = fitted.lm,
    "decomp.df" = decomp.df
  )

}
