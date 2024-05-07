HPAIboth <- read.csv2("./data/HPAIboth.csv", sep = ",")
save(HPAIboth, file = "./data/HPAIboth.rda")

HPAIcloaca <- read.csv2("./data/HPAIcloaca.csv", sep = ",")
save(HPAIcloaca, file = "./data/HPAIcloaca.rda")

HPAItrachea <- read.csv2("./data/HPAItrachea.csv", sep = ",")
save(HPAItrachea, file = "./data/HPAItrachea.rda")

HPAIfs <- read.csv2("./data/HPAIfs.csv", sep = ",")
save(HPAIfs, file = "./data/HPAIfs.rda")

library(tidyverse)
HPAIcloaca%>%reframe(.by = c(Vac, Experiment),
                     fs = 5-min(s))

HPAItrachea%>%reframe(.by = c(Vac, Experiment),
                     fs = 5-min(s))

HPAIfs

