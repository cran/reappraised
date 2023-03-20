# functions --------------------------------------------------------
stop_mb <- function(x) {
  blnk <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "))
  a <- paste0(blnk,x)
  stop(simpleError(a))
}

compress <- function(x, o= "obs", e= "exp", sfx = "", size=5) {
  #compress down into groups >5- suffix sm refers to small
  x$obs1 <- x [, o]; x$exp1 <- x [, e]
  count <- 1; ob <- 0; ex <- 0; start <- 0; x$obsm <- NA; x$exsm <- NA;  x$start <- NA; x$end <- NA; mx <- nrow(x)

  for (i in 1:mx) {
    xz <- 0
    if (ex + x$exp1 [i] > (size- 0.05)) {
      x$exsm [count] <- ex +x$exp1 [i]
      x$obsm [count] <- ob +x$obs1 [i]
      x$start [count] <- start
      x$end [count] <- i-1
      start <- i
      ex <- 0; ob <- 0; count <- count +1;
      xz <- 1
    }
    if (xz != 1) {
      xz <- 0
      ex <- ex+ x$exp1 [i]
      ob <- ob +x$obs1 [i]
    }
    if(i == mx) {
      if(count == 2) {
        x$exsm[count] <-ex
        x$obsm[count] <- ob
        x$start[count] <- start
        x$end[count] <- i
      } else {
        x$exsm [count-1] <- ex + x$exsm[count-1]
        x$obsm [count-1] <- ob + x$obsm[count-1]
        x$end [count-1] <- i}
    }
  }
  x$range <- ifelse(x$start == x$end, x$start, ifelse(x$end == max(x$end, na.rm = TRUE), paste(x$start,"+", sep=""),
                                                      paste(x$start,"-",x$end,sep="")))
  #chisquare test
  x$p_XS  <- NA
  x$p_XS [1] <-suppressWarnings(chisq.test(x$obsm[!is.na(x$obsm)], p= x$exsm[!is.na(x$exsm)]/sum(x$exsm, na.rm = TRUE)) [[3]])

  x <- x %>% dplyr::select("exsm", "obsm", "start", "end", "range", "p_XS")

  colnames(x) <- purrr::map(colnames(x), ~ paste0(.x,"_",sfx,"_",size))

  return (x)
}

#expected probability functions
expect_p <- function (df, mx_n, mx_g, nm_n, rw_g) {
  a <- df;
  exp_p <- matrix(0, nrow = mx_n *mx_g, ncol = nrow(a))
  for (r in 1:NROW(a)) {# for each row
    for (i in 1:a$group[r]) {#number of groups per row
      s1 <- a[r, nm_n[i]] #n1
      p1 <- a$p[r] #p
      exp_p [rw_g [[i]], r] <- pbinom(0,s1,p1) # prob of 0 and fits into matrix
      exp_p [(rw_g [[i]]+1):(rw_g [[i]]+s1), r] <- pbinom(1:s1,s1,p1) - pbinom((1:s1)-1,s1,p1) #prob of 1:s1
    }
  }
  return(exp_p)}

expect_p_diff <- function (df, ref_p, mx_n, nm_n) {
  a.2 <- df;
  exp_p <- matrix(0, nrow = mx_n, ncol = nrow(a.2))
  for (r in 1:NROW(a.2)) {# for each row
    s1 <- a.2[r, nm_n[1]] #n1
    s2 <- a.2[r, nm_n[2]] #n2
    p1 <- ref_p[1:mx_n, r] #exp grp1
    p2 <- ref_p[(mx_n+1):(mx_n*2), r] #exp grp2

    #to speed it up limit only to those without p=0
    p_min <- min(which.max(p1>0),  which.max(p2>0))
    p_max <- max( max(which(p1>0)),max(which(p2>0)))
    s <- p_max-p_min
    t <-data.table::as.data.table(data.frame(s3 = rep((p_min:p_max)-1, each = s+1),
                                             s4 = rep((p_min:p_max)-1, times =s+1)) %>%
                        dplyr::mutate (s5= abs(s3-s4)))
    #create all possible values and differences for non-zeros
    data.table::setorder(t,s5)
    t$prob <- p1[t$s3+1] * p2[t$s4+1]
    t <- t %>% dplyr::select(s5,prob) %>% dplyr::group_by (s5) %>% dplyr::summarise(sm = sum(prob), .groups = 'drop')
    exp_p [1:nrow(t) , r] <- t$sm
  }
  return(exp_p)}


#graph functions

graph_template <- function() {
  g1 <- list(
    ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black")),
    ggplot2::theme(text = ggplot2::element_text(size = 10, face= "bold")),
    ggplot2::theme(axis.line = ggplot2::element_line(size = 1), axis.ticks = ggplot2::element_line(size = 1)),
    ggplot2::theme(legend.background = ggplot2::element_blank(), legend.key = ggplot2::element_blank(),
                   legend.position= "top"),
    ggplot2::theme(legend.text = ggplot2::element_text(size=8, face="bold")),
    ggplot2::theme(axis.text.y = ggplot2::element_text(size=10, face="bold")),
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=10, face="bold")))
  return (g1)
}

#pval_graph
pval_graph <- function (x = p, round = "y", l=labels, t= title, ref="n", rp= "") {
  g1 <- graph_template()

  p <- x
  p$exp <- sum(p$Freq)/10
  p$p_XS  <- NA
  p$p_XS [1] <-suppressWarnings(chisq.test(p$Freq, p= rep(0.1,10))) [[3]]

  #compare to carlisle/Auck ref
  refp <- c(0.116, 0.097, 0.091, 0.096, 0.093, 0.095, 0.091, 0.086, 0.072, 0.163)
  p$p_XS_ref <- NA
  p$p_XS_ref [1] <-suppressWarnings(chisq.test(p$Freq, p= refp)) [[3]]

  p$prop <- p$Freq/sum(p$Freq)
  if(round == "y") {p$ref <- refp
  } else {p$ref <-  rep(0.1,10)}

  if(round != "y") {p$p_XS_ref [1] <- p$p_XS [1]}

  #for cat p vals
  if (ref == "y") {
    refp <- rp$prop
    p$ref <- refp
    p$p_XS_ref [1] <-suppressWarnings(chisq.test(p$Freq, p= refp)) [[3]]
  }

  p$x <- seq(0.05,0.95,0.1)

  #xaxis every 2nd value
  x.a <- as.character(seq(0,1,0.1)); x.a [1] <- "0.0"
  x.a <- ifelse(x.a %in% as.character(seq(0.1,0.9,0.2)), "", x.a)

  if (p$p_XS_ref [1] <0.001) {p1 <- format(p$p_XS_ref [1], scientific = TRUE, digit = 2)
  } else {p1 <- p$p_XS_ref [1]}
  p.e <- gregexpr("e",p1) [[1]]
  pval <- ifelse(p.e != -1,
                 paste("'p=",substring(as.character(p1),1,p.e-1),"*'*",
                       "10^",substring(as.character(p1),p.e+1,nchar(p1)),sep=""),
                 paste("'p='~'",signif(p1,2),"'", sep=""))

  #set y limits
  zx=1.05
  zy=max(p$prop)
  zy1 <- ifelse(zy < 0.2, 0.25, ifelse(zy < 0.4, 0.5, ifelse(zy < 0.6, 0.75, ifelse(zy < 0.8, 0.9, 1))))
  zz1 <- c(.05,.1,0.15,0.15,0.2) [zy1 ==c(0.25,0.5,0.75,0.9,1)]

  g <-
    ggplot2::ggplot (p) +
    g1 +
    ggplot2::geom_line(ggplot2::aes (x= x , y= ref, group = 1), linetype = "dashed", size = 1)+
    ggplot2::geom_bar(ggplot2::aes (x= x, y=prop), stat = "identity", fill = NA, colour = "black", size=1) +
    ggplot2::labs (x = "p-value", y = "Proportion") +  #axis labels
    ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black")) +
    ggplot2::scale_x_continuous(limits= c(0,1), breaks = seq (0, 1, by = 0.1), labels = x.a) +  #axis values
    ggplot2::scale_y_continuous(limits= c(0,zy1), breaks = seq(0, zy1, zz1)) +
    ggplot2::coord_fixed(ratio = zx/zy1*.67) +
    ggplot2::annotate(geom="text", zx*0, zy1 *1, label= t, color="black", size=3.25, hjust=0) +
    ggplot2::annotate(geom="text", zx*0, zy1 *.95, label= l [[1]], color="black", size=3.25, hjust=0) +
    ggplot2::annotate(geom="text", zx*0, zy1 *.9, label= l [[2]], color="black", size=3.25, hjust=0) +
    ggplot2::annotate(geom="text", zx*0, zy1 *0.85, label= l [[3]], color="black", size=3.25, hjust=0) +
    ggplot2::annotate(geom="text", zx*0, zy1 *0.80, label= pval, color="black", size=3.25,  hjust=0, parse = TRUE)

  return(list(g, pval))
}


#AUC graph function
auc_grph <- function (x, l = labels, round = "y", ref = "n", ra = "") {
  g1 <- graph_template()
  auc <- x

  #for graph using 100 p-values
  carl <- data.frame(p  =  c(0, 0.00038, 0.00601, 0.01417, 0.02368, 0.03376, 0.0438, 0.05406, 0.06497, 0.07489, 0.0849, 0.09427, 0.10405, 0.11462, 0.12479,
                           0.13417, 0.145, 0.15598, 0.16544, 0.17595, 0.18687, 0.19666, 0.20705, 0.21843, 0.22844, 0.23898, 0.25041, 0.2627, 0.27344,
                           0.28421, 0.29579, 0.30521, 0.31544, 0.32652, 0.33694, 0.34751, 0.35781, 0.36674, 0.377, 0.38858, 0.39979, 0.41114, 0.42245,
                           0.43415, 0.44467, 0.45471, 0.46568, 0.47651, 0.48781, 0.4973, 0.50724, 0.51737, 0.5272, 0.53633, 0.54752, 0.55795, 0.56982,
                           0.58058, 0.59111, 0.60231, 0.61261, 0.62385, 0.63426, 0.64564, 0.65624, 0.6673, 0.67913, 0.69026, 0.70139, 0.71389, 0.72632,
                           0.73741, 0.74794, 0.75936, 0.77035, 0.78226, 0.79351, 0.80455, 0.8176, 0.82962, 0.84425, 0.8573, 0.87191, 0.88838, 0.90368,
                           0.91921, 0.93448, 0.95094, 0.97011, 0.99465, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                   y = c(0, 0.00997, 0.01997, 0.02998, 0.03998, 0.04999, 0.05999, 0.06999, 0.08, 0.09, 0.09997, 0.10998, 0.11998, 0.12998, 0.13999,
                         0.14999, 0.16, 0.17, 0.18001, 0.18998, 0.19998, 0.20998, 0.21999, 0.22999, 0.24, 0.25, 0.26, 0.27001, 0.27998, 0.28998,
                         0.29999, 0.30999, 0.31999, 0.33, 0.34, 0.35001, 0.36001, 0.36998, 0.37998, 0.38999, 0.39999, 0.41, 0.42, 0.43004, 0.44001,
                         0.45001, 0.45998, 0.46999, 0.47999, 0.49, 0.5, 0.51, 0.52001, 0.53001, 0.54002, 0.54999, 0.55999, 0.56999, 0.58, 0.59, 0.60001,
                         0.61001, 0.62002, 0.63002, 0.63999, 0.64999, 0.66, 0.67, 0.68001, 0.69001, 0.70001, 0.71002, 0.72002, 0.72999, 0.74, 0.75,
                         0.76, 0.77001, 0.78001, 0.79002, 0.80002, 0.81002, 0.81999, 0.83, 0.84, 0.85001, 0.86001, 0.87002, 0.88002, 0.89002, 0.90003,
                         0.91, 0.92, 0.93001, 0.94001, 0.95001, 0.96002, 0.97009, 0.98006, 0.99013, 1))

  if (round == "y") {y1 <- carl
  } else {y1 <- auc; y1$p <- y1$y}

  if (ref == "y") {y1 <- ra}

  x.a <- as.character(seq(0,1,0.1)); x.a [1] <- "0.0"
  x.a <- ifelse(x.a %in% as.character(seq(0.1,0.9,0.2)), "", x.a)

  auc.g <-
    ggplot2::ggplot () +
    g1 +
    ggplot2::geom_line(data= auc, ggplot2::aes (x= p , y= y, group = 1), linetype = "solid", size = 1)+
    ggplot2::geom_line(data= y1, ggplot2::aes (x= p, y= y, group=1), linetype= "dashed", size=1) +
    ggplot2::labs (x = "p-value", y = "Cumulative proportion") +  #axis labels
    ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black")) +
    ggplot2::scale_y_continuous(limits= c(0.0,1), breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_x_continuous(limits= c(0,1), breaks = seq (0, 1, by = 0.1), labels = x.a) +  #axis values
    ggplot2::coord_fixed(ratio=1/1*0.67)+
    ggplot2::annotate(geom="text", x=0.0, y=1, label= l [[1]], color="black", size=3.25, hjust=0) + #,  fontface = "bold",
    ggplot2::annotate(geom="text", x=0.0, y=0.95, label= l [[2]], color="black", size=3.25,  hjust=0) +
    ggplot2::annotate(geom="text", x=0.0, y=0.9, label= l [[3]], color="black", size=3.25,  hjust=0) +
    ggplot2::annotate(geom="text", x=0.0, y=0.85, label= l [[4]], color="black", size=3.25,  hjust=0, parse = TRUE)+
    ggplot2::annotate(geom="text", x=0.0, y=0.79, label= l [[5]], color="black", size=3.25,  hjust=0)
  return (auc.g)
}


cat_graph <- function (gph= oe,  xtitle = "", ytitle = "", size = 5, sfx = "",
                       text = list(), fn = "", ti = "No", top = "No", t= title) {

  g1 <- graph_template()

  gph$leg <- "Expected"

  vars <- paste0(c("obsm","exsm","range","p_XS"),"_",sfx,"_",size)

  gph <- within(gph, {
    gp <- gph [, vars[4]]
    grange <- gph [, vars[3]]
    gexp <- gph [, vars[2]]
    gobs <- gph [, vars[1]]
  })

  if (gph$gp [1] <0.001) {p <- format(gph$gp [1], scientific = TRUE, digit = 2)
  } else {p <- gph$gp [1]}
  p.e <- gregexpr("e",p) [[1]]

  pval <- ifelse(p.e != -1,
                 paste("'p=",substring(as.character(p),1,p.e-1),"*'*",
                       "10^",substring(as.character(p),p.e+1,nchar(p)),sep=""),
                 paste("'p='~'",signif(p,2),"'", sep=""))

  abs <-
    ggplot2::ggplot (gph[!is.na(gph$gobs),]) + g1+
    ggplot2::geom_bar(ggplot2::aes (x= factor(grange, levels = grange) , y= gobs, fill = leg),
                      stat = "identity", colour = "black", width = 0.2, size=1) +
    ggplot2::geom_point(ggplot2::aes (x= factor(grange, levels = grange) , y= gexp),size=2) +
    ggplot2::geom_line(ggplot2::aes (x= factor(grange, levels = grange) , y= gexp, group = 1,
                                     linetype = "Expected"), size=.75)+
    ggplot2::labs (x = xtitle, y = ytitle, linetype= "", fill="") +  #axis labels
    ggplot2::scale_x_discrete(limits = gph$range[!is.na(gph$range)]) +
    ggplot2::scale_fill_manual(values =NA, labels= "Observed") +
    ggplot2::theme(legend.position = c(.9,.85)) +
    ggplot2::theme(legend.spacing.y = ggplot2::unit(-.05, "npc")) +
    ggplot2::theme(legend.key.width = ggplot2::unit(0.4, 'cm'), legend.key.height = ggplot2::unit(0.4, 'cm'))

  zx=max(ggplot2::ggplot_build(abs)$layout$panel_params[[1]]$x.range)
  zy=max(gph$gobs, gph$gexp, na.rm = TRUE)

  zy1 <- ifelse(zy <= 6.5, 8,
                ifelse(zy < 10, 15,
                       ifelse(zy < 20, 25,
                              ifelse(zy < 30, 35,
                                     ifelse(zy < 40, 50,
                                            ifelse(zy < 60, 75,
                                                   ifelse(zy < 75, 90,
                                                          ifelse(zy < 100, 125,150))))))))

  zz1 <- c(8,15,25,35,50,75,90,125,150)
  zz2 <- c(2,3,5,5,10,15,15,25,25)
  zz3 <- which(zz1 == zy1)
  abs <- abs +  ggplot2::scale_y_continuous(limits= c(0,zy1), breaks = seq(0, zy1, zz2 [zz3])) +
    ggplot2::coord_fixed(ratio = zx/zy1*.67)

  if (sum(!is.na(gph$grange))>5) {abs <- abs + ggplot2::theme(axis.text.x= ggplot2::element_text(angle=90, vjust=0.5))}

  #add text
  abs_all <- abs

  #add title
  if (tolower(substr(ti,1,1)) == "y") {abs <- abs +
    ggplot2::annotate(geom="text", x=zx*.2, y=zy1*.95, label= t, color="black", size=3.25, hjust=0) }

  if (fn == "cohort") {
    #individual graphs have all text
    abs <- abs +
      ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.95, label= paste(text[[1]], "cohorts"), color="black", size=3.25, hjust=0) +
      ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.85, label= paste(text[[2]], "variables"), color="black", size=3.25,  hjust=0) +
      ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.75, label= pval, color="black", size=3.25,  hjust=0, parse = TRUE) +
      ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.65, label= paste(text[[3]], "simulations"), color="black", size=3.25,  hjust=0)

    #now create all graphs- top graph has title plus everything. Otherwise just p-values
    if (tolower(substr(top,1,1)) =="y") {
      if (tolower(substr(ti,1,1)) == "y") {abs_all <- abs_all +
        ggplot2::annotate(geom="text", x=zx*.2, y=zy1*.95, label= t, color="black", size=3.25, hjust=0) }

      abs_all <- abs_all +
        ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.95, label= paste(text[[1]], "cohorts"), color="black", size=3.25, hjust=0) +
        ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.85, label= paste(text[[2]], "variables"), color="black", size=3.25,  hjust=0) +
        ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.75, label= pval, color="black", size=3.25,  hjust=0, parse = TRUE) +
        ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.65, label= paste(text[[3]], "simulations"), color="black", size=3.25,  hjust=0)
    } else {
      abs_all <- abs_all +
        ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.95, label= pval, color="black", size=3.25,  hjust=0, parse = TRUE)
    }
  }

  if (fn == "cat") {
    #individual graphs have all text
    if (length (text) == 1) {
      abs <- abs +
        ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.95, label= paste(text[[1]], "trials"), color="black", size=3.25, hjust=0) +
        ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.85, label= pval, color="black", size=3.25,  hjust=0, parse = TRUE)
    } else {
      abs <- abs +
        ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.95, label= paste(text[[1]], "trials"), color="black", size=3.25, hjust=0) +
        ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.85, label= paste(text[[2]], "trial arms"), color="black", size=3.25,  hjust=0) +
        ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.75, label= pval, color="black", size=3.25,  hjust=0, parse = TRUE)}

    #now create all graphs- top graph has title plus everything. Otherwise just p-values
    if (tolower(substr(top,1,1)) =="y") {
      if (tolower(substr(ti,1,1)) == "y") {abs_all <- abs_all +
        ggplot2::annotate(geom="text", x=zx*.2, y=zy1*.95, label= t, color="black", size=3.25, hjust=0) }

      if (length (text) == 1) {
        abs_all <- abs_all +
          ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.95, label= paste(text[[1]], "trials"), color="black", size=3.25, hjust=0) +
          ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.85, label= pval, color="black", size=3.25,  hjust=0, parse = TRUE)
      } else {
        abs_all <- abs_all +
          ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.95, label= paste(text[[1]], "trials"), color="black", size=3.25, hjust=0) +
          ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.85, label= paste(text[[2]], "trial arms"), color="black", size=3.25,  hjust=0) +
          ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.75, label= pval, color="black", size=3.25,  hjust=0, parse = TRUE)}

      } else {
        abs_all <- abs_all +
          ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.95, label= pval, color="black", size=3.25,  hjust=0, parse = TRUE)
      }
    }

  if (fn == "cat_all") {
    #individual graphs have all text
    abs <- abs +
      ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.95, label= paste(text[[1]], "trials"), color="black", size=3.25, hjust=0) +
      ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.85, label= text[[2]], color="black", size=3.25,  hjust=0) +
      ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.75, label= pval, color="black", size=3.25,  hjust=0, parse = TRUE)

    #now create all graphs- top graph has title plus everything. Otherwise just p-values
    if (tolower(substr(top,1,1)) =="y") {
      if (tolower(substr(ti,1,1)) == "y") {abs_all <- abs_all +
        ggplot2::annotate(geom="text", x=zx*.2, y=zy1*.95, label= t, color="black", size=3.25, hjust=0) }

      abs_all <- abs_all +
        ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.95, label= paste(text[[1]], "trials"), color="black", size=3.25, hjust=0) +
        ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.85, label= text[[2]], color="black", size=3.25,  hjust=0) +
        ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.75, label= pval, color="black", size=3.25,  hjust=0, parse = TRUE)

    } else {
      abs_all <- abs_all +
        ggplot2::annotate(geom="text", x=zx*.5, y=zy1*.95, label= pval, color="black", size=3.25,  hjust=0, parse = TRUE)
    }
  }

  return(list(abs, abs_all))
}

rr_graph <- function(rr= oe,  xtitle = "", ytitle = "", size = 5, sfx = "",
                     text = list(), fn = "", ti = "No", top = "No", t= title) {

  g1 <- graph_template()

  vars <- paste0(c("obsm","exsm","range"),"_",sfx,"_",size)

  rr <- within (rr, {
    range <- rr [,vars[3]]
    e <- rr [, vars[2]]
    o <- rr [, vars[1]]
    n <- sum(o, na.rm = TRUE)
    estimate <- (o/n/(e/n))
    n1 <- e+o
    total <- n*2
    n2 <- total - n1
    conf.level <- 0.95
    norm.pp <- qnorm(1 - (1 - conf.level)/2)
    p.v <- 2 * (1 - pnorm(abs((o - n1 * n/total)/sqrt(n1 * n2 * n * n/total/total/(total - 1)))))
    RRL <- estimate * exp(-norm.pp * sqrt(1/o - 1/n + 1/e - 1/n))
    RRU <- estimate * exp(norm.pp * sqrt(1/o - 1/n + 1/e - 1/n))
  })

  zy <- c(0.01, 0.05, 0.1, .25,0.4, 1, 2, 4, 8, 16)
  zyl <- ifelse(min(rr$RRL[rr$RRL >0], na.rm=TRUE) <0.01, 0.01, max(zy[zy- min(rr$RRL[rr$RRL >0], na.rm = TRUE)<0]))
  zyh <- ifelse(max(rr$RRU, na.rm = TRUE) >16, 16, min (zy[zy- max(rr$RRU, na.rm = TRUE)>0]))

  oe <- ggplot2::ggplot (rr[!is.na(rr$o),]) + g1+
    ggplot2::geom_point(ggplot2::aes (x= range , y= estimate),size=2)+
    ggplot2::geom_errorbar(ggplot2::aes(x= range, ymin=RRL, ymax=RRU), width=.1,
                           position= ggplot2::position_dodge(0.05), size=0.75) +
    ggplot2::geom_line(ggplot2::aes (x= range , y= estimate, group = 1), size=.75)+
    ggplot2::geom_line(ggplot2::aes (x= range , y= 1, group = 1), linetype = "dashed", size = 0.5)+
    ggplot2::labs (x = xtitle, y = ytitle) +  #axis labels
    ggplot2::scale_x_discrete(limits =rr$range[!is.na(rr$range)]) +
    ggplot2::scale_y_continuous(limits= c(zyl,zyh), breaks = zy) +
    ggplot2::coord_trans(y="log2")

  zx=max(ggplot2::ggplot_build(oe)$layout$panel_params[[1]]$x.range)
  zyr = max(ggplot2::ggplot_build(oe)$layout$panel_params[[1]]$y.range)-
    min(ggplot2::ggplot_build(oe)$layout$panel_params[[1]]$y.range)

  if (sum(!is.na(rr$range))>5) {oe <- oe + ggplot2::theme(axis.text.x= ggplot2::element_text(angle=90, vjust= 0.5))}

  oe <- oe + ggplot2::theme(aspect.ratio=.67)

  #add text
  oe_all <- oe

  #add title
  if (tolower(substr(ti,1,1)) == "y") {oe <- oe +
    ggplot2::annotate(geom="text", x=zx*.2, y=zyl*2^0.85, label= t, color="black", size=3.25, hjust=0) }

  if (fn == "cohort") {
    #individual graphs have all text
    oe <- oe +
      ggplot2::annotate(geom="text", x=zx*.2, y=zyl*1.05, label= paste(text[[1]], "cohorts"), color="black", size=3.25, hjust=0) +
      ggplot2::annotate(geom="text", x=zx*.4, y=zyl*1.05, label= paste(text[[2]], "variables"), color="black", size=3.25,  hjust=0) +
      ggplot2::annotate(geom="text", x=zx*.4, y=zyl*2^.85, label= paste(text[[3]], "simulations"), color="black", size=3.25,  hjust=0)

    #all graphs have no text
  }

  if (fn == "cat") {
    #individual graphs have all text
    if (length (text) == 1) {
      oe <- oe +
        ggplot2::annotate(geom="text", x=zx*.2, y=zyl*1.05, label= paste(text[[1]], "trials"), color="black", size=3.25,  hjust=0)
    } else {
      oe <- oe +
        ggplot2::annotate(geom="text", x=zx*.2, y=zyl*1.05, label= paste(text[[1]], "trials"), color="black", size=3.25,  hjust=0) +
        ggplot2::annotate(geom="text", x=zx*.2, y=zyl*2^.5, label= paste(text[[2]], "trial arms"), color="black", size=3.25,  hjust=0)}
    }
    return(list(oe, oe_all))
}

digit_g <- function (x = p, l=labels, t= title) {
  g1 <- graph_template()

  p <- x

  p$exp <- sum(p$Freq)/10
  p$p_XS  <- NA
  p$p_XS [1] <-chisq.test(p$Freq, p= rep(0.1,10)) [[3]]

  p$prop <- p$Freq/sum(p$Freq)
  p$ref <-  rep(0.1,10)
  p$x <- seq(0.05,0.95,0.1)

  if (p$p_XS [1] <0.001) {p1 <- format(p$p_XS [1], scientific = TRUE, digit = 2)
  } else {p1 <- p$p_XS [1]}
  p.e <- gregexpr("e",p1) [[1]]
  pval <- ifelse(p.e != -1,
                 paste("'p=",substring(as.character(p1),1,p.e-1),"*'*",
                       "10^",substring(as.character(p1),p.e+1,nchar(p1)),sep=""),
                 paste("'p='~'",signif(p1,2),"'", sep=""))

  #set y limits
  zy=max(p$prop)
  zy1 <- ifelse(zy < 0.125, 0.15, ifelse(zy < 0.2, 0.25, ifelse(zy < 0.4, 0.5, ifelse(zy < 0.8, 0.9, 1))))
  zz1 <- c(.03,.05,0.1,0.15,0.2) [zy1 ==c(0.15,0.25,0.5,0.9,1)]

  g <- ggplot2::ggplot (p) +
    g1 +
    ggplot2::geom_line(aes (x= Var1 , y= ref, group = 1), linetype = "dashed", size = 1)+
    ggplot2::geom_bar(aes (x= Var1, y=prop), stat = "identity", fill = NA, colour = "black", size=1) +
    ggplot2::labs (x = "final digit", y = "Proportion") +  #axis labels
    ggplot2::theme(panel.border = ggplot2::element_blank(), panel.grid.major = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(colour = "black")) +
    ggplot2::scale_y_continuous(limits= c(0,zy1), breaks = seq(0, zy1, zz1)) +
    ggplot2::coord_fixed(ratio = 10/zy1*.67) +
    ggplot2::annotate(geom="text", 1, zy1 *1, label= t, color="black", size=3.25, hjust=0) +
    ggplot2::annotate(geom="text", 1, zy1 *.95, label= l [[1]], color="black", size=3.25, hjust=0) +
    ggplot2::annotate(geom="text", 1, zy1 *.9, label= l [[2]], color="black", size=3.25, hjust=0) +
    ggplot2::annotate(geom="text", 1, zy1 *0.85, label= pval, color="black", size=3.25,  hjust=0, parse = TRUE)

  return(list(g, pval))
}








