optionsWindow <- function (env) {
  
  env$optsWindow <- gwindow("Options", visible=FALSE, height=50, width=50)
  env$optsGroup <- ggroup(container=env$optsWindow, horizontal=FALSE)
  env$l <- glayout(container=env$optsGroup, homogeneous=TRUE, spacing=30)
  env$l[1, 1] <- (env$l1 <- glayout(container=env$l, spacing=2))
  env$l[1, 2] <- (env$l2 <- glayout(container=env$l, spacing=2))
  
  env$l1[1, 1:3] <- (env$opts$headings$t1 <- glabel("Graph sizes", container=env$l1))
  env$l1[2, 2, anchor=c(0, 0)] <- (env$opts$text$t2 <- glabel("height", container=env$l1))  
  env$l1[2, 3, anchor=c(0, 0)] <- (env$opts$text$t1 <- glabel("width", container=env$l1))
  env$l1[3, 2] <- (env$opts$spHeight <- gspinbutton(from=250, to=500, by=10, 
                                                    value=250, digits=0, 
                                                    container=env$l1))
  env$l1[3, 3] <- (env$opts$width <- gspinbutton(from=500, to=800, by=10, 
                                                 value=500, digits=0, 
                                                 container=env$l1))
  env$l1[3, 1, anchor=c(-1, 0)] <- (env$opts$textBf$t1 <- glabel("Spectrum graph", container=env$l1))
  env$l1[4, 1, anchor=c(-1, 0)] <- (env$opts$textBf$t2 <- glabel("Chromatogram", container=env$l1))
  env$l1[4, 2] <- (env$opts$chHeight <- gspinbutton(from=250, to=500, by=10, 
                                                    value=250, digits=0, 
                                                    container=env$l1))  
  i <- 5
  
  env$l1[i, 1:3] <- (env$opts$separator$t1 <- glabel("", container=env$l1))
  env$l1[i + 1, 1:3] <- (env$opts$headings$t3 <- glabel("Labels", container=env$l1))
  env$l1[i + 2, 1:2] <- (env$opts$text$t1 <- glabel("Number of peaks to label", container=env$l1))
  env$l1[i + 2, 3] <- (env$opts$labels <- gedit(5, coerce.with=as.numeric, width=2, container=env$l1))
  
  env$l2[1, 1:3] <- (env$opts$headings$t2 <- glabel("MS mode", container=env$l2))
  env$l2[2, 2, anchor=c(1, 0)] <- (env$opts$text$t4 <- glabel("MS1     ", container=env$l2))
  env$l2[2, 3, anchor=c(1, 0)] <- (env$opts$text$t5 <- glabel("MS2     ", container=env$l2))
  env$l2[3:4, 2] <- (env$opts$ms1 <- gradio(items=c("", " "), selected=2, container=env$l2))
  env$l2[3:4, 3, anchor=c(1, 0)] <- (env$opts$ms2 <- gradio(items=c("", " "), selected=1, container=env$l2))
  env$l2[3, 1, anchor=c(-1, 0)] <- (env$opts$text$t6 <- glabel("Centroided             ", container=env$l2))
  env$l2[4, 1, anchor=c(-1, 0)] <- (env$opts$text$t7 <- glabel("Profile", container=env$l2))
  i <- 5
  env$l2[i, 1:3] <- (env$opts$separator$t2 <- glabel("", container=env$l2))
  env$l2[i + 1, 1:3] <- (env$opts$headings$t3 <- glabel("Chromatogram mode", container=env$l2))
  env$l2[i + 2, 1:3] <- (env$opts$chromaMode <- gradio(items=c("Total ion count", 
                                                               "Base peak intensity"), 
                                                       selected=1, container=env$l2))
  env$l2[i + 3, 1:3] <- (env$opts$separator$t3 <- glabel("", container=env$l2))
  
  env$optsButtons <- ggroup(container=env$optsGroup, horizontal=TRUE)
  addSpring(env$optsButtons)
  env$opts$ok <- gbutton("OK", container=env$optsButtons, width=50)
  env$opts$cancel <- gbutton("Cancel", container=env$optsButtons, width=50, 
                             handler=function(h, ...) dispose(env$optsWindow))
  env$opts$apply <- gbutton("Apply", container=env$optsButtons, width=50)
  addSpace(env$optsButtons, 18)
  visible(env$optsWindow) <- TRUE
}
e <- new.env()
optionsWindow(e)
