length(getHook("before.plot.new"))

setHook("before.plot.new", value=list(function() {
  devices <- names(dev.list())
  if("RStudioGD"%in%devices) 
    dev.set(which(devices=="RStudioGD"))
  else if(length(devices[devices!="Cairo"])==0) dev.new()
}))

x <- dev.list()
class(x)
names(x)


setHook("before.plot.new", value=list)