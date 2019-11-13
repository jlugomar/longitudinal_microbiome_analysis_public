#########################################################################################################################
#########################################################################################################################
###### PROJECT:        DBNs
###### NAME:           ReplaceNamesInNetwork.R
###### AUTHOR:         Daniel Ruiz-Perez, PhD Student
###### AFFILIATION:    Florida International University
###### 
###### DESCRIPTION:    This file is for when the names of the taxa inside the graphs are incorrect
#########################################################################################################################
#########################################################################################################################



library(scales)
library(stringr)
library(igraph)
options("scipen"=100, "digits"=4)

setwd("/Users/elyonyi/Documents/GitHub/longitudinal_microbiome_analysis_public/visualization/")

folder = "normalization" #"Demo" #"finalFiguresBootscore" #"finalFiguresBootscore"
#nameOfNetwork = "DemoFigure.graphml"#list.files(pattern =".*\\.graphml")
# nameOfNetwork = "human_pregnancy_microbiota_tooth_gum_dbn_sample_alignment_sr21d_dbnIntraBoot.graphml"#list.files(pattern =".*\\.graphml")
# nameOfNetwork = paste(folder,nameOfNetwork,sep="/")




files = list.files(path = folder, pattern = ".graphml")

for (file in files){
  
  nameOfNetwork = paste(folder,file,sep="/")

  
  f = readChar(nameOfNetwork, file.info(nameOfNetwork)$size)
  
  f = gsub("1023075", "Veillonella", f)
  f = gsub("4475758", "Veillonella dispar", f)
  f = gsub("4453535", "Gemellaceae", f)
  f = gsub("4479989", "Streptococcus1", f)
  # f = gsub("4425214", "Streptococcus2", f)
  f = gsub("4483015", "Granulicatella1", f)
  f = gsub("4409545", "Porphyromonas1", f)
  f = gsub("1008348", "Prevotella", f)
  f = gsub("4477696", "Haemophilus parainfluenzae", f)
  # f = gsub("4401186", "Eikenella", f)
  # f = gsub("4480775", "Kingella", f)
  f = gsub("4442130", "Streptococcus2", f)
  f = gsub("4483174", "Leptotrichia", f)
  f = gsub("4310398", "Prevotella nanceiensis", f)
  f = gsub("4467992", "Streptococcus3", f)
  f = gsub("4447416", "Neisseria", f)
  f = gsub("4396235", "Neisseria subflava", f)
  f = gsub("4452538", "Fusobacterium", f)
  f = gsub("4479603", "Prevotella melaninogenica", f)
  f = gsub("4443574", "Haemophilus", f)
  f = gsub("4465803", "Porphyromonas2", f)
  f = gsub("1696853", "Granulicatella2", f)
  
  
  
  #save it and take care of the null character that appears
  writeChar(f,nameOfNetwork,nchars = nchar(f))
  r = readBin(nameOfNetwork, raw(), file.info(nameOfNetwork)$size)
  r[r==as.raw(0)] = as.raw(0x20) ## replace with 0x20 = <space>
  writeBin(r, nameOfNetwork)

}