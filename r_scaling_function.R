## This function is used to normalize the length of some vector of numerical values.
## Mainly used to have density data from transcribed regions of the same length.

# Data      : input vector of numerical values
# FinalLength  : [100] the length of the resulting vector
# Extention : [FALSE] should the input vector be considered in three parts (extension-main-extention) :  i.e gene body + flanking regions
# Ext_length   : [c(300, 300)] size of the extentions (base-pair) from the input vector
# Ext_value : [c(10, 10)] size to which the extentions should be reduced
# method : ["mean"] should the resulting vector be the mean or the median of the values

Scale_Vector=function(Data, FinalLength=100, Extention=FALSE, Ext_length=c(300, 300), Ext_value=c(10,10), method="mean"){
  if(length(Ext_length)==1){tmp_val=Ext_length;Ext_length=rep(tmp_val, 2)}
  if(length(Ext_value)==1){tmp_val=Ext_value;Ext_value=rep(tmp_val, 2)}
  # browser()
  BPnorm=c(1,ceiling(c(1:FinalLength)*length(Data)/FinalLength))
  if(Extention){
    #FinalLength=100; Extention=TRUE; Ext_length=300; Ext_value=10;method="mean"
    (step=c(Ext_length[1]/Ext_value[1], (length(Data)-(sum(Ext_length)))/(FinalLength-(sum(Ext_value))), Ext_length[2]/Ext_value[2]))
    (start=ceiling(seq(1, Ext_length[1]+1-step[1], step[1])))
    (middle=ceiling(seq(start[length(start)]+step[1], length(Data)-Ext_length[2]+(1*step[2]), step[2])))

    if(length(middle)<FinalLength-(sum(Ext_value))+1){
      (middle=ceiling(seq(start[length(start)]+step[1], length(Data)-Ext_length[2]+step[2], length.out=FinalLength-(sum(Ext_value))+1)))
    }

    (end=ceiling(seq(middle[length(middle)]+step[3], length(Data)+step[3], step[3])))
    length(c(start, middle, end))
    (BPnorm=c(start, middle, end))
  }
  if(method=="mean"){(Normed=sapply(c(1:FinalLength), function(x) mean(Data[BPnorm[x]:BPnorm[x+1]], na.rm=T)))}
  if(method=="median"){(Normed=sapply(c(1:FinalLength), function(x) mean(Data[BPnorm[x]:BPnorm[x+1]], na.rm=T)))}
  return(Normed)
}
