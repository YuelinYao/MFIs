### UI
MBUI <- function(){
  tagList(
    tags$h3(paste0("Markov blanket"), style = "color: steelblue;"),
    plotOutput(outputId ="MarkovBlanket", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("MarkovBlanket_plot","Download as .pdf")
  )}


### Input function
MBInput<- function(){
  tagList( 
    textInput("NodeGene", "Input Gene name:",value = "IGHG4"),
    actionButton("action_mb","Submit",icon("paper-plane"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  )
}



## Function
# The markov blanket should include the the parent, children, and the spouse
mb <- function(graph, nodes = igraph::V(graph)) {
  lapply(
    seq_along(nodes),
    function(i) {
      nodes_out <- igraph::neighborhood(
        graph,
        nodes[[i]],
        order = 1,
        mode = "out"
      )[[1]] %>% unlist()%>%names()%>%
        setdiff(nodes[[i]]) %>% ## children 
        
        #print()%>% 
        igraph::neighborhood(graph, ., order = 1, mode = "in") %>% #print()%>% #the parent of your children
        unlist()%>%names()#%>% print()
      #parents
      igraph::neighborhood(graph, nodes[[i]], order = 1, mode = "in")[[1]] %>%unlist()%>%names()%>%
        union(nodes_out) %>% 
        setdiff(nodes[[i]]) 
    }
  )
}
