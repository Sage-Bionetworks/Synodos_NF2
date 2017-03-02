library(synapseClient)
library(ggplot2)
library(dplyr)
library(reshape2)
library(viridis)
synapseLogin()

x<-synGet("syn8314523")@filePath
drugdat<-read.table(x, sep = '\t', header = TRUE, fill = TRUE)

colnames(drugdat) <- c('Protocol.Name', 'Sample.ID', 'Cell.Line', 'Cell.Type', 
                       'NF2.Status', 'AC50.uM', 'CRC', 'Max.Resp', 'Log.AC50.uM', 
                       'Compound.Name', 'AUC', 'AUC.Fit', 'Gene.Name')

sub<-select(drugdat, Cell.Line, Max.Resp, Compound.Name, Gene.Name)

sub2 <- sub %>% 
  group_by(Compound.Name) %>% 
  nest() %>% 
  mutate(num_cells = map(data, nrow)) %>% 
  filter(num_cells == 6) %>% 
  mutate(data = map(data, function(x) {
    spread(x, Cell.Line, Max.Resp)
    })) %>% 
  unnest(data)

ggplot(sub2, aes(x= HS01, y=HS12)) +
    geom_point(stat = "identity", aes(color=(HS01/HS12))) + 
    coord_fixed() +
    scale_color_viridis() +
    geom_text(data=subset(sub2, HS01/HS12>4), aes(label = Compound.Name), size = 2)

intersect(filter(sub2, HS01/HS12>1.1)$Compound.Name, filter(sub4, HS12/HS01>1.1)$Compound.Name)


ggplot(sub2, aes(x= MS02, y=HS12)) +
  geom_point(stat = "identity", aes(color=(MS02/HS12))) + 
  coord_fixed() +
  scale_color_viridis() +
  geom_text(data=subset(sub2, MS02/HS12>4), aes(label = Compound.Name), size = 2)

ggplot(sub2, aes(x= Syn5, y=Syn1)) +
  geom_point(stat = "identity", aes(color=(Syn5/Syn1))) + 
  coord_fixed() +
  scale_color_viridis() +
  geom_text(data=subset(sub2, Syn5/Syn1>4), aes(label = Gene.Name), size = 2)

ggplot(sub2, aes(x= Syn6, y=Syn1)) +
  geom_point(stat = "identity", aes(color=(Syn6/Syn1))) + 
  coord_fixed() +
  scale_color_viridis() +
  geom_text(data=subset(sub2, Syn6/Syn1>4), aes(label = Compound.Name), size = 2)


sub3 <- select(drugdat, Cell.Line, AUC, Compound.Name)
sub4 <- sub3 %>% 
  group_by(Compound.Name) %>% 
  nest() %>% 
  mutate(num_cells = map(data, nrow)) %>% 
  filter(num_cells == 6) %>% 
  mutate(data = map(data, function(x) {
    spread(x, Cell.Line, AUC)
  })) %>% 
  unnest(data)

ggplot(sub4, aes(x= HS01, y=HS12)) +
  geom_point(stat = "identity", aes(color=(HS12/HS01))) + 
  coord_fixed() +
  scale_color_viridis() +
  geom_text(data=subset(sub4, HS12/HS01>1.5), aes(label = Compound.Name), size = 2)

ggplot(sub4, aes(x= MS02, y=HS12)) +
  geom_point(stat = "identity", aes(color=(HS12/MS02))) + 
  coord_fixed() +
  scale_color_viridis() +
  geom_text(data=subset(sub4, HS12/MS02>1.5), aes(label = Compound.Name), size = 2)

ggplot(sub4, aes(x= Syn5, y=Syn1)) +
  geom_point(stat = "identity", aes(color=(Syn1/Syn5))) + 
  coord_fixed() +
  scale_color_viridis() +
  geom_text(data=subset(sub4, Syn1/Syn5>1.5), aes(label = Compound.Name), size = 2)

ggplot(sub4, aes(x= Syn6, y=Syn1)) +
  geom_point(stat = "identity", aes(color=(Syn1/Syn6))) + 
  coord_fixed() +
  scale_color_viridis() +
  geom_text(data=subset(sub4, Syn1/Syn6>1.5), aes(label = Compound.Name), size = 2)
