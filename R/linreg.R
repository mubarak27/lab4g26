#' A class for linear regression with many methods for corresponding values calculations
#'
#' @field y vector The orginial labels of observations in data
#' @field X matrix The features of all observations in data
#' @field beta vector The coefficients in the sample
#' @field y_lb matrix The estimated labels
#' @field x_lb matrix The standard errors between orginial and estimated labels
#' @field dof numeric The degree of freedom of data
#' @field var_rv numeric The residual variance
#' @field var_rcf matrix The variance of the regression coefficients
#' @field t_vcf vector The t-values for each coefficient
#' @field p_vcf vector The p-values for each coefficient
#' @field formula character The formula used in the model
#' @field data data.frame The sample data
#'
#' @import methods
#' @importFrom ggplot2 ggplot
#' @importFrom plyr is.formula
#' @importFrom gridExtra grid.arrange
#' @export linreg
#' @exportClass linreg

linreg <- setRefClass("linreg",
  fields = list(
    y = "numeric",
    X = "matrix",
    y_lb = "matrix",
    beta = "vector",
    x_lb = "matrix",
    dof = "numeric",
    var_rv = "numeric",
    var_rcf = "matrix",
    t_vcf = "vector",
    p_vcf = "vector",
    formula = "formula",
    data = "data.frame",
    name = "character",
    f = "character"
  ),
  methods = list(
    # The initialize methods
    initialize = function(formula, data){
      stopifnot(plyr::is.formula(formula),is.data.frame(data) )
      cat("User::initialize")
      f <<- Reduce(paste, deparse(formula))
      formula <<- formula
      name <<- deparse(substitute(data))
      data <<- data
      X <<- model.matrix(formula,data)
      y <<- data[,all.vars(formula)[1]]
      beta <<- round(as.vector(solve(t(X)%*%X)%*%t(X)%*%y),3)
      y_lb <<- round(as.matrix(X%*%beta),3)
    x_lb <<- y-y_lb
      dof <<- length(y) - ncol(X)
      var_rv <<- as.numeric((t(x_lb)%*%x_lb)/dof)
      var_rcf <<- var_rv * solve(t(X) %*% X)
      t_vcf <<- beta/(sqrt(diag(var_rcf)))
      p_vcf <<- 2*pt(abs(t_vcf), dof,lower.tail = FALSE)
      names(p_vcf) <<- colnames(X)
      names(beta) <<- colnames(X)
    },
    #Estimation of beta and their variances by QR decomposition
    qr_beta = function(){
      Q <- qr.Q(qr(X))
      R <- qr.R(qr(X))
      beta <<- solve(R)%*%t(Q)%*%y
      var_rcf <<- var_rv*solve(R)%*%t(solve(R))
    },
    #Regressions coefficients
    coef = function(){
      return(beta)
    },
    #The fitted values
    pred = function(){
      return(y_lb)
    },
    #The residuals
    resid = function(){
      return(x_lb)
    },
    #The degrees of freedom
    freedomdegree = function(){
      return(dof)
    },
    #The residual variance
    residualvariance = function(){
      return(var_rv)
    },
    #The variance of the regression coefficients
    coefvariance = function(){
      return(var_rcf)
    },
    #The t-values for each coefficient
    t_vcfues = function(){
      return(t_vcf)
    },
    #The p-values for each coefficient
    p_vcfues = function(){
      return(p_vcf)
    },
    #The print function
    print = function(){
      cat(paste("linreg(formula = ",f,", data = ", name,")\n",sep=""))
      base::print(beta)
    },
    #The summary function
    summary = function(){
      std_e <- round((sqrt(diag(var_rcf))),3)
      m <- as.data.frame(matrix(c(round(beta,3), std_e, round(t_vcf,3), format(p_vcf, scientific = 2)), ncol = 4))
      for(i in 1:ncol(X)){
        if(p_vcf[i]>=0&&p_vcf[i]<0.001){
          m[i,5]<-"***"
        }else if(
          p_vcf[i]>=0.001&&p_vcf[i]<0.01){
          m[i,5]<-"**"
        }else if(
          p_vcf[i]>=0.01&&p_vcf[i]<0.05){
          m[i,5]<-"*"
        }else if(
          p_vcf[i]>=0.05&&p_vcf[i]<0.1){
          m[i,5]<-"."
        }else if(
          p_vcf[i]>=0.1&&p_vcf[i]<1){
          m[i,5]<-""
        }
      }
      colnames(m) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)","")
      rownames(m) <- colnames(X)
      base::print(m)
      base::print(paste("Residual standard error:", round(sqrt(var_rv),3), "on", dof, "degrees of freedom"))
    },

    #The plot function
    plot=function(){
      d <- as.data.frame(cbind(y_lb,x_lb, sqrt(abs((x_lb-mean(x_lb))/sd(x_lb)))))

      g1<-ggplot2::ggplot(data=d,ggplot2::aes(x=d[,1],y=d[,2])) +
        ggplot2::geom_point(shape=1, size=4)+
        ggplot2::stat_summary(ggplot2::aes(x=y_lb,group=1),fun.y=median, color="red", geom="line")+
        ggplot2::geom_hline(yintercept = 0, col="grey", linetype="dotted")+
        ggplot2::labs(x="Fitted valies or Predictions",y="Residuals")+
        ggplot2::ggtitle("Residuals vs Fitted Plot")+
        ggplot2::scale_y_continuous()+
        ggplot2::scale_color_discrete("Regression")+
        ggplot2::theme(
          panel.background  = ggplot2::element_blank(),
          plot.background = ggplot2::element_rect(fill="lightskyblue1", color=NA),
          legend.background = ggplot2::element_rect(fill="transparent", color=NA),
          legend.key = ggplot2::element_rect(fill="transparent", color=NA)
        )


      g2<-ggplot2::ggplot(d,ggplot2::aes(y=d[,3],x=d[,1]))+
        ggplot2::geom_point(shape=1, size=4)+
        ggplot2::stat_summary(ggplot2::aes(x=y_lb,group=1),fun.y=median, color="red", geom="line")+
        ggplot2::labs(x="Fitted valies or Predictions",y=" sqrt(|Standardized Residuals|)")+
        ggplot2::ggtitle("sqrt(|Standardized Residuals|) vs Fitted Plot")+
        ggplot2::scale_y_continuous()+
        ggplot2::scale_color_discrete("Regression")+
        ggplot2::theme(
          panel.background  = ggplot2::element_blank(),
          plot.background = ggplot2::element_rect(fill="lightskyblue1", color=NA),
          legend.background = ggplot2::element_rect(fill="transparent", color=NA),
          legend.key = ggplot2::element_rect(fill="transparent", color=NA)
        )

      gg<-gridExtra::grid.arrange(g1, g2)

    }
  )

)

