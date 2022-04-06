library(readr)
library(dplyr)

#IMPORT ALL SETS OF DATA, RENAME COLUMNS AND DROP FIRST COLUMN
SLC22A2 <- read_csv("./data/network/SLC22A2.csv") %>%
  select(-"type")
SLC22A2<-data.frame(append(SLC22A2, c(transporter="SLC22A2"))) 

SLC22A12 <- read_csv("./data/network/SLC22A12.csv") %>%
  select(-"type")
SLC22A12<-data.frame(append(SLC22A12, c(transporter="SLC22A12"))) 

SLC22A6 <- read_csv("./data/network/SLC22A6.csv") %>%
  select(-"type")
SLC22A6<-data.frame(append(SLC22A6, c(transporter="SLC22A6"))) 

SLC22A8 <- read_csv("./data/network/SLC22A8.csv") %>%
  select(-"type")
SLC22A8<-data.frame(append(SLC22A8, c(transporter="SLC22A8"))) 

SLC22A3 <- read_csv("./data/network/SLC22A3.csv") %>%
  select(-"type")
SLC22A3<-data.frame(append(SLC22A3, c(transporter="SLC22A3"))) 

SLC22A16 <- read_csv("./data/network/SLC22A16.csv") %>%
  select(-"type")
SLC22A16<-data.frame(append(SLC22A16, c(transporter="SLC22A16"))) 

SLC22A7 <- read_csv("./data/network/SLC22A7.csv") %>%
  select(-"type")
SLC22A7<-data.frame(append(SLC22A7, c(transporter="SLC22A7"))) 

SLC22A11 <- read_csv("./data/network/SLC22A11.csv") %>%
  select(-"type")
SLC22A11<-data.frame(append(SLC22A11, c(transporter="SLC22A11"))) 

SLC22A10 <- read_csv("./data/network/SLC22A10.csv") %>%
  select(-"type")
SLC22A10<-data.frame(append(SLC22A10, c(transporter="SLC22A10"))) 

SLC22A1 <- read_csv("./data/network/SLC22A1.csv") %>%
  select(-"type")
SLC22A1<-data.frame(append(SLC22A1, c(transporter="SLC22A1"))) 

SLC22A4 <- read_csv("./data/network/SLC22A4.csv") %>%
  select(-"type")
SLC22A4<-data.frame(append(SLC22A4, c(transporter="SLC22A24"))) 

SLC22A5 <- read_csv("./data/network/SLC22A5.csv") %>%
  select(-"type")
SLC22A5<-data.frame(append(SLC22A5, c(transporter="SLC22A5")))


#BIND ALL COLUMNS TOGETHER
total <- bind_rows(SLC22A1,SLC22A2,SLC22A3,SLC22A4,SLC22A5)

library(tidyverse)

#get distinct transporters
sources <- total %>%
  distinct(transporter)%>%
  rename(label=transporter)

#get distinct drug names
destinations <- total %>%
  distinct(properties.generic.name)%>%
  rename(label=properties.generic.name)

# Join the two data to create node
# Add unique ID
nodes <- full_join(sources, destinations, by = "label") 
nodes <- nodes %>%
  mutate(id = 1:nrow(nodes)) %>%
  select(id, everything())
head(nodes, 3)

nodes$group <- c(rep("A",5), rep("B",169))

per_route <- total %>%
  group_by(transporter, properties.generic.name) %>%
  summarise(weight = n()) %>%
  ungroup()

#THIS FORMULA TURN COLUMNS INTO VARIABLES
set_lists_to_chars <- function(x) {
  if(class(x) == 'list') {
    y <- paste(unlist(x[1]), sep='', collapse=', ')
  } else {
    y <- x 
  }
  return(y)
}

#CREATE A NEW DATA FRAME
total <- data.frame(lapply(per_route, set_lists_to_chars), stringsAsFactors = F)


# (a) Join nodes id for source column
edges <- total %>%
  left_join(nodes,by = c("transporter"= "label")) %>%
  rename(from=id)

# (b) Join nodes id for destination column
edges <- edges %>%
  left_join(nodes, by = c("properties.generic.name" = "label")) %>%
  rename(to=id)

edges <- select(edges,from,to,weight)
head(edges,3)



library(igraph)

#BUILD THE EDGEs AND VERTICES
net.igraph <- graph_from_data_frame(
  d= edges, vertices= nodes,
  directed= FALSE)
net.igraph
net.igraph[]

#SEE THE FIRST PLOT
plot(net.igraph,
     vertex.label.family="Helvetica",
     vertex.label.color="black",
     edge.arrow.size=1,
     edge.arrow.width=2,
     vertex.label.cex=0.5,
     remove.multiple = F, remove.loops = T)

#add a layout for it to stop moving everytime I run the code again
V(net.igraph)$size <- log(strength(net.igraph)) * 4 + 3
l= layout_with_fr(net.igraph)

par(mar=c(0,0,0,0)); plot(net.igraph,
                          vertex.label.family="Helvetica",
                          vertex.label.color="black",
                          vertex.color= "orange",
                          edge.arrow.size=1,
                          edge.arrow.width=2,
                          vertex.label.cex=0.5,
                          remove.multiple = F, remove.loops = T,
                          layout=l)

####spread the graph
l <- layout_with_fr(net.igraph)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
plot(net.igraph, rescale=F, layout=l*1.2)

##
l <- layout_with_kk(net.igraph)
plot(net.igraph, layout=l)
plot(net.igraph, layout=layout_with_dh, vertex.label.cex=0.5,)

# Community detection based on label propagation:
clp <- cluster_label_prop(net.igraph)
class(clp)

# Community detection returns an object of class "communities" # which igraph knows how to plot:
plot(clp, net.igraph)

# We can also plot the communities without relying on their built-in plot:
V(net.igraph)$community <- clp$membership
colrs <- adjustcolor( c("gray50", "tomato", "gold", "yellowgreen"), alpha=.6)
plot(net.igraph, vertex.color=colrs[V(net)$community])


#### ANIMATION!!
library('visNetwork')

vis.nodes <- nodes
vis.edges <- edges

graph <- visNetwork(vis.nodes, vis.edges,
                    width="100%", height="800px", background="#eeefff",
                    main="SLC22 family drug network",
                    submain="Ñ=",
                    footer="source:Uniprot") %>%
  visGroups(groupname = "A", color= "tomato", shape = "diamond") %>%
  visLayout(randomSeed=123) %>%
  visPhysics(maxVelocity = 20) %>%
  visEdges(arrows = "to")
graph


# We'll start by adding new node and edge attributes to our dataframes. 
#THIS PART SHOULD CREATE A DROPDOWN MANU TO SELECT THE LABEL
#changing some details 
vis.links$width <- 1+edges$weight/8 # line width
vis.links$color <- "gray"    # line color  
vis.links$arrows <- "middle" # arrows: 'from', 'to', or 'middle'
vis.links$smooth <- FALSE    # should the edges be curved?
vis.links$shadow <- FALSE    # edge shadow

visnet <- visNetwork(vis.nodes, vis.links) %>%
  visPhysics(maxVelocity = 20)
visnet

visOptions(visnet, highlightNearest = TRUE, selectedBy = "label")

#CLUSTER BY LOUVAIN
plot(net.igraph, 
     vertex.size = 5,
     mark.groups = cluster_louvain(net.igraph), 
     vertex.color = "grey80",
     vertex.border.color = "grey60",
     vertex.label.cex = 0.5,
     vertex.label.color = "black")

# SAVE NETWORK GRAPH
visSave(graph, file = "network.html", background = "white")