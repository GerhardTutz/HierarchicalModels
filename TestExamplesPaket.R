library(OrdHierarchical)
data(GLES)

# Model 1
hier <- OrdHierarchical(dat=GLES, k=7, nameresp="Terrorism", fittype="disp", tests="LR", modeltype="split", pred=c("Age","Gender","Abitur"))
hier
round(hier$pargroupedresponse[3:5,] ,digits=3)
round(hier$parnonneutra[5:10,] ,digits=3)

test_equality(hier, tests="LR")
test_dispersion(hier, tests="LR")

hier$loglikfullmodel
hier$AIC

# Model 2
hier <- OrdHierarchical(dat=GLES, k=7, nameresp="Terrorism", fittype="disp", tests="LR", modeltype="split", pred=c("Age","Gender","Abitur"), preddisp="Age")
hier
round(hier$pardispersionmodel[5:9,], digits=3)

hier$loglikdispersionmodel
hier$AICdisp

# Model 3
hier <- OrdHierarchical(dat=GLES, k=7, nameresp="Terrorism", fittype="disp", tests="LR", modeltype="split", pred=c("Age","Gender","Abitur"), preda=c("Age","Gender","Abitur"), preddis=c("Age"))
hier

hier$loglikfullmodel
hier$AIC

# Separation of neutral category
hier <- OrdHierarchical(dat=GLES, k=7, nameresp="Terrorism", fittype="disp", tests="LR", modeltype="neutral", cum=FALSE, pred=c("Age","Gender","Abitur"))
hier

round(hier$parneutral[2:4,], digits = 3)
round(hier$fitnonneutral$parnonneutral[5:10,], digits=3)

# reduced model
hier <- OrdHierarchical(dat=GLES, k=7, nameresp="Terrorism", fittype="disp", tests="LR", modeltype="neutral", cum=FALSE, pred=c("Age","Gender","Abitur"),
                        predneutral=c("Age","Gender"), preddis=c("Age"))
hier

# With cumulative model
hier <- OrdHierarchical(dat=GLES, k=7, nameresp="Terrorism", fittype="disp", tests="LR", modeltype="neutral", cum=TRUE, pred=c("Age","Gender","Abitur"))
hier

hier$loglikfullmodel
hier$AIC

# Alternative fitting tools
hier <- OrdHierarchical(dat=GLES, k=7, nameresp="Terrorism", fittype="clm", pred=c("Age","Gender","Abitur"))
hier

hier <- OrdHierarchical(dat=GLES, k=7, nameresp="Terrorism", fittype="vglm", pred=c("Age","Gender","Abitur"))
hier

# Separation of neutral category
hier <- OrdHierarchical(dat=GLES, k=7, nameresp="Terrorism", fittype="vglm", modeltype="neutral", cum=FALSE, pred=c("Age","Gender","Abitur"), family="cum")
hier

hier <- OrdHierarchical(dat=GLES, k=7, nameresp="Terrorism", fittype="vglm", modeltype="neutral", cum=TRUE, pred=c("Age","Gender","Abitur"), family="cum")
hier

# Fitting general hierarchical models
nodes <- list(list(c(1,2,3),4,c(5,6,7)),list(1,2,3),list(5,6,7))

hier <- OrdHierarchical(dat=GLES, k=7, nameresp="Terrorism", fittype="clm", pred=c("Age","Gender","Abitur"), modeltype="general", general_structure=nodes)
hier

structure_hier <- buildTree(hier)
print(structure_hier, "level")
plot(structure_hier)

nodes <- list(list(c(1,2,3),4,c(5,6,7)),list(1,2,3),list(5,6,7))
predlist <- list(c("Age","Gender","Abitur"),"Age","Age")

hier <- OrdHierarchical(dat=GLES, k=7, nameresp="Terrorism", fittype="clm", pred=c("Age","Gender","Abitur"), modeltype="general", general_structure=nodes,
                        general_pred=predlist)
hier

# Happiness
data(Happiness)

hier <- OrdHierarchical(dat=Happiness, k=10, nameresp="happicat", fittype="disp", tests="LR", modeltype="split", pred=c("Age","Gender"))
hier
round(hier$pargroupedresponse[2:3,] ,digits=3)
round(hier$parnonneutra[9:12,] ,digits=3)

test_equality(hier, tests="LR")
test_dispersion(hier, tests="LR")

hier$loglikfullmodel
hier$AIC

# Model with dispersion effects
hierdisp <- OrdHierarchical(dat=Happiness, k=10, nameresp="happicat", fittype="disp", tests="LR", modeltype="split", pred=c("Age","Gender"),
                        preddisp=c("Age","Gender"))
hierdisp



