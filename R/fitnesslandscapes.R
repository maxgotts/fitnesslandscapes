# Easy dot product
dotprod <- function(x,y) {
  x <- as.vector(x)
  y <- as.vector(y)
  z <- c()
  l <- max(length(x),length(y))
  for (i in 1:l) {
    z <- c(z,x[i]*y[i])
  }
  z<-z[!is.na(z)]
  r <- 0
  for (e in z) {
    r <- r+e
  }
  return(r)
}


# Performs PPR (and LDA) on a data set
PPR_LDA <- function(DF=df, exclude=c("Identifier"), include=FALSE, LDA=FALSE, FITNESS="Fitness") {
  # Exclude items
  if (LDA != FALSE) {
    groups <- DF[,LDA]
    DF[,LDA] <- NULL
  }
  if (exclude != FALSE && include != FALSE) {
    print("Ignoring exclude")
    return(PPR_df(DF, FALSE, include, LDA, FITNESS))
  } else if (exclude != FALSE && include == FALSE) {
    exclusion_df <- DF[,exclude]
    DF[,exclude] <- NULL
  } else if (exclude == FALSE && include != FALSE) {
    return(PPR_df(DF, setdiff(colnames(DF),include), FALSE, LDA, FITNESS))
  } else {
    exclusion_df <- data.frame()
  }

  # Perform LDA
  if (LDA != FALSE) {
    DF[,LDA] <- groups
    FORMULA <- as.formula(paste(LDA,paste(setdiff(colnames(DF),c(FITNESS, LDA)),collapse=" + "), sep=" ~ "))
    df.lda <- MASS::lda(FORMULA, DF)
    DF[,LDA] <- NULL
  }

  # Perform PPR
  df.norm <- as.data.frame(scale(DF))
  explanatory <- df.norm
  explanatory[,FITNESS] <- NULL
  explanatory <- as.matrix(explanatory)
  response <- as.matrix(df.norm[,FITNESS])
  df.ppr <- stats::ppr(explanatory,response,nterms=2,maxterms=5)
  pprdirections <- df.ppr$alpha

  # Add directions
  fitness <- DF[,FITNESS]
  DF[,FITNESS] <- NULL
  PP1 <- PP2 <- LD1 <- LD2 <- c()
  for (row in 1:nrow(DF)) {
    PP1 <- c(PP1, dotprod(pprdirections[,"term 1"],DF[row,]))
    PP2 <- c(PP2, dotprod(pprdirections[,"term 2"],DF[row,]))
    if (LDA != FALSE) {
      LD1 <- c(LD1, dotprod(df.lda$scaling[,"LD1"],DF[row,]))
      if (length(df.lda$scaling)>1) {
        LD2 <- c(LD2, dotprod(df.lda$scaling[,"LD2"],DF[row,]))
      }
    }
  }
  DF$PP1 <- PP1
  DF$PP2 <- PP2
  DF$LD1 <- LD1
  if (LDA != FALSE) {
    if (length(df.lda$scaling)>1) {
      DF$LD2 <- LD2
    }
  }
  PP1_vector <- pprdirections[,"term 1"]
  PP2_vector <- pprdirections[,"term 2"]

  # Add excluded items
  New_df <- cbind(DF,exclusion_df)
  New_df[,FITNESS] <- fitness

  # Case-wise return statements
  if (LDA != FALSE) {
    New_df[,LDA] <- groups
    return(list(
      New_df,
      df.ppr,
      df.lda,
      pprdirections,
      df,lda$scaling
    ))
  } else if (LDA == FALSE) {
    return(list(
      New_df,
      df.ppr,
      pprdirections
    ))
  }
}



# Returns a fitness landscape given specific parameters
TPS_landscape <- function(DF=df, VarName1="PP1", VarName2="PP2", output="contour", Theta=30, Phi=30, FitnessMetric="Fitness", DispVar1=VarName1, DispVar2=VarName2, DispFitness=FitnessMetric) {
  Var1 <- DF[,VarName1]
  Var2 <- DF[,VarName2]
  Fitness <- DF[,FitnessMetric]
  tp.m <- as.matrix(data.frame(v1=Var1, v2=Var2))
  t <- fields::Tps(x=tp.m, Y=Fitness,lambda=0.02691373)

  if (output=="plotly") return(plotly::plot_ly(z=~fields::predictSurface(t)$z) %>% plotly::add_surface())
  if (output=="contour") return(fields::surface(t, xlab=DispVar1, ylab=DispVar2,
                                        zlab=DispFitness))
  if (output=="wireframe") return(fields::surface(t, xlab=DispVar1, ylab=DispVar2,
                                          zlab=DispFitness,type="p", theta=Theta, phi=Phi))
  if (output=="matrix") return(fields::predictSurface(t)$z)
  if (output=="model") return(t)
}

# Creates a 2D frequency-binning
binCounts <- function(x,y,increment_x,increment_y, pdf=FALSE) {
  x_seq <- seq(min(x),max(x),increment_x)
  y_seq <- seq(min(y),max(y),increment_y)
  counts <- expand.grid(x=x_seq, y=y_seq)
  rownames(counts) <- paste(counts$x,counts$y,sep=":")
  counts$z <- 0
  for (xv in x) {
    for (yv in y) {
      xi <- floor((xv-min(x))/increment_x)*increment_x+min(x)
      yi <- floor((yv-min(y))/increment_y)*increment_y+min(y)
      counts[paste(xi,yi,sep=":"),"z"] <- counts[paste(xi,yi,sep=":"),"z"]+1
    }
  }
  colnames(counts) <- c("x","y","counts")
  if (pdf) {
    counts$counts <- counts$counts/sum(counts$counts)
  }
  return(counts)
}

# Creates a TPS density surface based on a binning
TPS_distribution <- function(DF=df,x="PP1",y="PP2",output,x_divisor=2,y_divisor=2,Theta=30,Phi=30,pdf=FALSE) {
  x_axis <- DF[,x]
  y_axis <- DF[,x]
  return(TPS_landscape(binCounts(x_axis, y_axis, increment_x=sd(x_axis)/3, increment_y=sd(y_axis)/3, pdf=pdf),
                       "x", "y", output, DispVar1=x, DispVar2=y,FitnessMetric="counts", DispFitness="Frequency"))
}

# Returns a data-frame of an ellipse-segment of a TPS model given specific parameters
ellipse_at <- function(center, radius, increment, tps_model, pdf) {
  pdf(1,1) # To check PDF function works before investing runtime
  p_x <- center[1]
  p_y <- center[2]
  radius_x <- radius[1]
  radius_y <- radius[2]
  increment_x <- increment[1]
  increment_y <- increment[2]
  min_x <- p_x-radius_x
  max_x <- p_x+radius_y
  min_y <- p_y-radius_x
  max_y <- p_y+radius_y
  x_seq <- seq(min_x, max_x, increment_x)
  y_seq <- seq(min_y, max_y, increment_y)
  ellipse <- expand.grid(x=x_seq, y=y_seq)
  ellipse$z <- NA
  for (r in 1:nrow(ellipse)) {
    # cat("sum",(ellipse[r,"x"]-p_x)^2+(ellipse[r,"y"]-p_y)^2,"radius",radius^2,"\n")
    if ((ellipse[r,"x"]-p_x)^2/radius_x^2+(ellipse[r,"y"]-p_y)^2/radius_y^2 <= 1) {
      ellipse[r,"z"] <- 0
    } else {
      ellipse[r,"z"] <- NA
    }
  }
  ellipse <- dyplr::filter(ellipse, !is.na(z))
  ellipse$z <- stats::predict(tps_model, ellipse[,c("x","y")])
  ellipse$p <- NA
  for (r in 1:nrow(ellipse)) {
    ellipse[r,"p"] <- pdf(ellipse[r,"x"]-p_x,ellipse[r,"y"]-p_y)
  }
  return(ellipse)
}
