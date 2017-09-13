library(synapseClient)
library(openxlsx)
library(dplyr)
synapseLogin()
sch <- read.table(synGet("syn9884855")@filePath, sep = "\t", header = T)

write.xlsx(sch %>% filter(comparison == "HS01CUDC907vsHS11CUDC907"), "HS01CUDC907vsHS11CUDC907.xlsx")
write.xlsx(sch %>% filter(comparison == "HS01DMSOvsHS11DMSO"), "HS01DMSOvsHS11DMSO.xlsx")
write.xlsx(sch %>% filter(comparison == "HS01GSK2126458vsHS11GSK2126458"), "HS01GSK2126458vsHS11GSK2126458.xlsx")
write.xlsx(sch %>% filter(comparison == "HS01panobinostatvsHS11panobinostat"), "HS01panobinostatvsHS11panobinostat.xlsx")


mgh <- read.table(synGet("syn6166525")@filePath, sep = "\t", header = T)
write.xlsx(mgh %>% filter(diffExptest == "Syn5.CUDC907-Syn1.CUDC907"), "Syn5.CUDC907-Syn1.CUDC907.xlsx")
write.xlsx(mgh %>% filter(diffExptest == "Syn5.GSK2126458-Syn1.GSK2126458"), "Syn5.GSK2126458-Syn1.GSK2126458.xlsx")
write.xlsx(mgh %>% filter(diffExptest == "Syn5.Panobinostat-Syn1.Panobinostat"), "Syn5.Panobinostat-Syn1.Panobinostat.xlsx")
write.xlsx(mgh %>% filter(diffExptest == "Syn5.DMSO-Syn1.DMSO"), "Syn5.DMSO-Syn1.DMSO.xlsx")

ms <- read.table(synGet("syn7437782")@filePath, sep = "\t", header = T) %>% filter(cellLine1 == "MS03" & cellLine2 == "MS12")

write.xlsx(ms %>% filter(treatment1 == "Pano" & treatment2 == "Pano"), "MS03-MS12.Pano.xlsx")
write.xlsx(ms %>% filter(treatment1 == "DMSO" & treatment2 == "DMSO"), "MS03-MS12.DMSO.xlsx")
write.xlsx(ms %>% filter(treatment1 == "CUDC" & treatment2 == "CUDC"), "MS03-MS12.CUDC.xlsx")
write.xlsx(ms %>% filter(treatment1 == "GSK458" & treatment2 == "GSK458"), "MS03-MS12.GSK.xlsx")

