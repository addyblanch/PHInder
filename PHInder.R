#MATCHING ANALYSIS single piece of code that will look for correlations between two data sets
#that come from the same samples. We are particularly thinking about bacteria and gene data from
#the same set of individuals.The code only looks at presence/absence rather than abundance.
library(reshape2)
#See the example data for appropriate input Both are data-frames with column 1
#representing the bacteria/gene/whatever is being tested. The remaining columns should
#represent each sample and should have the same labels to allow for matchings. Any columns
#whose name is unique to one of the data frames will be deleted. min.presence states how often
#an element should appear in order to be considered for analysis. The default setting of 1
#will remove all non-present/ever-present elements. Increasing this will reduce the number of elements
#considered but will increase the power of corellation tetss. Max presence is automatically set 
#to the number of samples-min.presence for symmetry but can be changed.
matchings.analysis = function(x,y,min.presence = 1,max.presence = NULL){
 
  #GET DATA IN APPROPRIATE FORMAT
  #Make Col 1 the row name
  rownames(x) = x[,1]
  x = x[,-1]
  rownames(y) = y[,1]
  y = y[,-1]
  
  #Remove samples that are not repeated across both datasets
  y = y[,colnames(x)]
  x = x[,colnames(y)]
  
  #Make logical and remove any elements that are (n)ever-present.Transpose matrix,
  x = as.data.frame(t(x) > 0) 
  y = as.data.frame(t(y) > 0) #Logical for presence and transpose 
  
  x = x[,colSums(x) >= min.presence]
  y = y[,colSums(y) >= min.presence] #Remove elements that are not present enough
  
  if (is.null(max.presence) == T){
    max.presence = nrow(x)-min.presence
  }
  x = x[,colSums(x) <= max.presence]
  y = y[,colSums(y) <= max.presence] #Remove elements that are not present enough
  
  #CREATE A DATA FRAME TO SHOW NUMBER OF MATCHES
  cors = vector(mode = "list", length = nrow(x))
  for (i in 1:nrow(x)){
    #Find all matchings
    cors[[i]] = as.data.frame(as.numeric(x[i,]) %*% t(as.numeric(y[i,]))) #All T/T
    rownames(cors[[i]]) = colnames(x)
    colnames(cors[[i]]) = colnames(y)
  }
  #Have now created a matrix for each sheep which gives T if both present or both absent
  matchings = Reduce('+', cors)
  rm(i,cors)
  
  #Create matchings data frame
  matchings$x_element <- rownames(matchings)
  matchings <- melt(matchings, id.vars = c("x_element"), variable.name = "y_element", value.name = "matches")
  matchings$x_count = colSums(x)[match(matchings$x_element,names(colSums(x)))] #How often does x occur
  matchings$y_count = colSums(y)[match(matchings$y_element,names(colSums(y)))] #How often does truth occur
  

  #Find the probability of these counts resulting in these number of matches or more. this is the raw p-value
  #Matching probs follow a hypergeometric distribution
  matchings$raw.p = phyper(matchings$matches-1,matchings$x_count,nrow(x)-matchings$x_count,matchings$y_count,lower.tail=FALSE)
  matchings=matchings[order(matchings$raw.p),]
  matchings$p.holm = p.adjust(matchings$raw.p,method = "holm")
  matchings$p.BH = p.adjust(matchings$raw.p,method = "BH")
  
  return(matchings) #A matrix giving every correlation in order of significance including
  #p-value adjustments for the Holm correction and Benjamini Hochberg procedure
  
}