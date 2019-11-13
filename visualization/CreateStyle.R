#########################################################################################################################
#########################################################################################################################
###### PROJECT:        DBNs
###### NAME:           CreateStyle.R
###### AUTHOR:         Daniel Ruiz-Perez, PhD Student
###### AFFILIATION:    Florida International University
###### 
###### DESCRIPTION:    This file creates the style.xml file needed to visualize it in Cytoscape. It takes as input a network
######                , a file with the mean abundance and a base style file to update.
#########################################################################################################################
#########################################################################################################################


library(scales)
library(stringr)
options("scipen"=100, "digits"=4)


setwd("/Users/elyonyi/Documents/GitHub/longitudinal_microbiome_analysis_public/visualization/")

nameOfNetwork = "infant_gut_microbiota_dbn_sample_alignment_filtered_sr3d_dbnIntraBoot_noPeriodOfStudy.graphml"#list.files(pattern =".*\\.graphml")
abundances = unlist(read.csv('meanAbundanceInfantFilteredAlignment.txt',sep="\t"))

# infant gut data set
names = c("Day of life sample obtained","Gestational age at birth","Postconceptional age sample obtained",
          "Gender","Mode of birth","Room type","Human milk used",
          "Days of antibiotics","Actinobacteria","Alphaproteobacteria","Bacilli","Bacteroidia",
          "Betaproteobacteria","Clostridia","Cyanobacteria","Epsilonproteobacteria","Erysipelotrichi","Flavobacteria",
          "Fusobacteria","Gammaproteobacteria","Holophagae","Unclassified")


# vaginal data set
#names = c("SubjectID","Day_Period","Aerococcus","Anaerococcus","Ureaplasma","Parvimonas","L. iners",
#          "Finegoldia","Staphylococcus","Porphyromonas","Atopobium","Gardnerella","L. crispatus","Peptostreptococcus",
#          "Sneathia","Streptococcus","Prevotella","Peptoniphilus","Incertae_Sedis_XI.1","L. gasseri","Corynebacterium",
#          "Dialister","L.otu5","L.otu3","Ruminococcaceae.3")


# For tooth gum
#names = c("Day sample obtained","Gestational day of delivery", "Veillonella","Veillonella dispar","Gemellaceae","Streptococcus1",
#          "Streptococcus2","Granulicatella","Porphyromonas","Prevotella","Haemophilus parainfluenzae","Eikenella","Kingella")


#non-taxa for infant gut data set
attr = c("Subject","Day of life sample obtained","Gestational age at birth","Postconceptional age sample obtained","Gender",
         "Mode of birth","Room type","Human milk used","Days of antibiotics")

#non-taxa for vaginal data set
#attr = c("Nugent_Category",	"Day_Period")

#non-taxa for tooth gum data set
#attr = c("Day sample obtained","Gestational day of delivery")


f = readChar("styleBase.xml", file.info("styleBase.xml")$size)

split = strsplit(f, "<dependency name=\"nodeSizeLocked\" value=\"true\"/>", fixed = T)
before = unlist(split)[1]
after = unlist(split)[2]





################## NODE COLOR
aux ="\n            <visualProperty name=\"NODE_FILL_COLOR\" default=\"#ff9232\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">"
for (i in 1:length(names)){
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"#4d93a8\" attributeValue=\"",names[i],"_ti+1","\"/>",sep="")
}

aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")
f = paste(before,"<dependency name=\"nodeSizeLocked\" value=\"true\"/>",aux,after,sep= "")


############### SIZE
network = readChar(nameOfNetwork, file.info(nameOfNetwork)$size)
network= unlist(strsplit(x =network, "\r\n"))
network= unlist(strsplit(x =network, "\n"))
maxAcum = 1
for (i in 1:length(names))
  maxAcum = max(maxAcum, length(grep(pattern = paste(".*target=\"",names[i],"_ti+1","\".*",sep=""), x = network, ignore.case = T)))

aux =paste(aux,"\n            <visualProperty name=\"NODE_SIZE\" default=\"40\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
for (i in 1:length(names)){
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"",(20*length(grep(pattern = paste(".*target=\"",names[i],"_ti","\".*",sep=""), x = network, ignore.case = T))/maxAcum+40),"\" attributeValue=\"",names[i],"_ti","\"/>",sep="")
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"",(20*length(grep(pattern = paste(".*target=\"",names[i],"_ti\\+1","\".*",sep=""), x = network, ignore.case = T))/maxAcum+40),"\" attributeValue=\"",names[i],"_ti+1","\"/>",sep="")
}
aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")
f = paste(before,"<dependency name=\"nodeSizeLocked\" value=\"true\"/>",aux,after,sep= "")



############### TRANSPARENCY
abundancesScaled = rescale(abundances,c(90,255))
abundancesScaled = c(rep(200,length(attr)),abundancesScaled)
aux =paste(aux,"\n            <visualProperty name=\"NODE_TRANSPARENCY\" default=\"255\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
for (i in 1:length(names)){
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"",round(abundancesScaled[i],1),"\" attributeValue=\"",names[i],"_ti","\"/>",sep="")
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"",round(abundancesScaled[i],1),"\" attributeValue=\"",names[i],"_ti+1","\"/>",sep="")
}
aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")
f = paste(before,"<dependency name=\"nodeSizeLocked\" value=\"true\"/>",aux,after,sep= "")



############### SHAPE
aux =paste(aux,"\n            <visualProperty name=\"NODE_SHAPE\" default=\"ELLIPSE\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
for (i in 1:length(attr)){
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"DIAMOND\" attributeValue=\"", paste(attr[i],"_ti",sep=""),"\"/>",sep="")
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"DIAMOND\" attributeValue=\"", paste(attr[i],"_ti+1",sep=""),"\"/>",sep="")
}
aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")
f = paste(before,"<dependency name=\"nodeSizeLocked\" value=\"true\"/>",aux,after,sep= "")


split = strsplit(f, "<dependency name=\"arrowColorMatchesEdge\" value=\"false\"/>", fixed = T)
before = unlist(split)[1]
after = unlist(split)[2]

################## EDGE LINE TYPE
aux = "\n            <visualProperty name=\"EDGE_LINE_TYPE\" default=\"SOLID\">\n                <discreteMapping attributeType=\"string\" attributeName=\"shared name\">"
for (i in 1:length(names)){
  for (j in 1:length(names)){
    aux = paste(aux,"\n                    <discreteMappingEntry value=\"EQUAL_DASH\" attributeValue=\"",names[i],"_ti+1 (-) ", names[j],"_ti+1","\"/>",sep="")
    aux = paste(aux,"\n                    <discreteMappingEntry value=\"EQUAL_DASH\" attributeValue=\"",names[i],"_ti (-) ", names[j],"_ti","\"/>",sep="")
  }
}
aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")


################## EDGE PAINT
# aux =paste(aux,"\n            <visualProperty name=\"EDGE_UNSELECTED_PAINT\" default=\"#CC0033\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
# for (i in 1:length(names)){
#     aux = paste(aux,"\n                    <discreteMappingEntry value=\"#FF9999\" attributeValue=\"",names[i],"_ti (-) ", names[i],"_ti+1","\"/>",sep="")
#   }
# 
# aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")


# ##### Make the self loops more transparent
# aux =paste(aux,"\n            <visualProperty name=\"EDGE_TRANSPARENCY\" default=\"255\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
# for (i in 1:length(names)){
#     aux = paste(aux,"\n                    <discreteMappingEntry value=\"70\" attributeValue=\"",names[i],"_ti (-) ", names[i],"_ti+1","\"/>",sep="")
#   }
# 
# aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")
# 


################ Now transparency is based on the boot score
aux =paste(aux,"\n            <visualProperty name=\"EDGE_TRANSPARENCY\" default=\"255\">\n                <discreteMapping attributeType=\"string\" attributeName=\"name\">",sep="")
for (i in 1:length(names)){
  aux = paste(aux,"\n                    <discreteMappingEntry value=\"70\" attributeValue=\"",names[i],"_ti (-) ", names[i],"_ti+1","\"/>",sep="")
}

aux = paste(aux,"\n                </discreteMapping>\n            </visualProperty>" ,sep="")



################## EDGE TRANSPARENCY
weights = str_extract(network[grep("key=\"key_bootScore",network)], "\\-*\\d+\\.+\\d+")
weights = weights[!is.na(weights)]
maxAbsWeight =  max(abs(as.numeric(weights)))
minAbsWeight =  min(abs(as.numeric(weights)))
# maxAbsWeight = 150 #when we hide intra edges we need this

######## For edge coefficient
aux =paste(aux,"\n<visualProperty name=\"EDGE_TRANSPARENCY\" default=\"2.0\">\n")
aux =paste(aux,"  <continuousMapping attributeType=\"float\" attributeName=\"bootScore\">\n")
aux =paste(aux,"  <continuousMappingPoint lesserValue=\"130.0\" greaterValue=\"130.0\" equalValue=\"130.0\" attrValue=\"",minAbsWeight,"\"/>\n",sep="")
aux =paste(aux,"  <continuousMappingPoint lesserValue=\"130.0\" greaterValue=\"255.0\" equalValue=\"255.0\" attrValue=\"",maxAbsWeight,"\"/>\n",sep="")
aux =paste(aux,"  </continuousMapping>\n")
aux =paste(aux,"  </visualProperty>\n")



################## EDGE WIDTH
weights = str_extract(network[grep("key=\"key_weight",network)], "\\-*\\d+\\.+\\d+")
weights = weights[!is.na(weights)]

maxAbsWeight =  max(abs(as.numeric(weights)))
medianAbsWeight =  median(abs(as.numeric(weights)))
# maxAbsWeight = 150 #when we hide intra edges we need this



######## For edge coefficient normalized
aux =paste(aux,"\n<visualProperty name=\"EDGE_WIDTH\" default=\"2.0\">\n")
aux =paste(aux,"  <continuousMapping attributeType=\"float\" attributeName=\"weight\">\n")
aux =paste(aux,"  <continuousMappingPoint lesserValue=\"15.0\" greaterValue=\"15.0\" equalValue=\"15.0\" attrValue=\"-",1,"\"/>\n",sep="")
aux =paste(aux,"  <continuousMappingPoint lesserValue=\"1.0\" greaterValue=\"1.0\" equalValue=\"1.0\" attrValue=\"0.0\"/>\n")
aux =paste(aux,"  <continuousMappingPoint lesserValue=\"15.0\" greaterValue=\"15.0\" equalValue=\"15.0\" attrValue=\"",1,"\"/>\n",sep="")
aux =paste(aux,"  </continuousMapping>\n")
aux =paste(aux,"  </visualProperty>\n")



# ######## For edge coefficient
# aux =paste(aux,"\n<visualProperty name=\"EDGE_WIDTH\" default=\"2.0\">\n")
# aux =paste(aux,"  <continuousMapping attributeType=\"float\" attributeName=\"weight\">\n")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"10.0\" greaterValue=\"15.0\" equalValue=\"15.0\" attrValue=\"-",maxAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"1.0\" greaterValue=\"10.0\" equalValue=\"10.0\" attrValue=\"-",medianAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"1.0\" greaterValue=\"1.0\" equalValue=\"1.0\" attrValue=\"0.0\"/>\n")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"1.0\" greaterValue=\"10.0\" equalValue=\"10.0\" attrValue=\"",medianAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"10.0\" greaterValue=\"15.0\" equalValue=\"15.0\" attrValue=\"",maxAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  </continuousMapping>\n")
# aux =paste(aux,"  </visualProperty>\n")
  
########### For edge confidence
# aux =paste(aux,"<visualProperty name=\"EDGE_WIDTH\" default=\"2.0\">\n")
# aux =paste(aux,"  <continuousMapping attributeType=\"float\" attributeName=\"weight\">\n")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"16.0\" greaterValue=\"20.0\" equalValue=\"20.0\" attrValue=\"-",maxAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"8.0\" greaterValue=\"16.0\" equalValue=\"16.0\" attrValue=\"-",80.0,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"8.0\" greaterValue=\"8.0\" equalValue=\"8.0\" attrValue=\"-",20.0,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"1.0\" greaterValue=\"1.0\" equalValue=\"1.0\" attrValue=\"0.0\"/>\n")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"8.0\" greaterValue=\"8.0\" equalValue=\"8.0\" attrValue=\"",20.0,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"8.0\" greaterValue=\"16.0\" equalValue=\"16.0\" attrValue=\"",80.0,"\"/>\n",sep="")
# aux =paste(aux,"  <continuousMappingPoint lesserValue=\"16.0\" greaterValue=\"20.0\" equalValue=\"20.0\" attrValue=\"",maxAbsWeight,"\"/>\n",sep="")
# aux =paste(aux,"  </continuousMapping>\n")
# aux =paste(aux,"  </visualProperty>\n")


fileConn<-file("style.xml")
writeLines(paste(before,"<dependency name=\"arrowColorMatchesEdge\" value=\"false\"/>",aux,after,sep=""), fileConn)
close(fileConn)
