### UI
UPsetUI <- function(){
  tagList(
    tags$h3(paste0("Upset Plot"), style = "color: steelblue;"),
    plotOutput(outputId ="UPsetPlot", width = "90%") %>% withSpinner(color="#4682B4"),
    downloadButton("upset_plot","Download as .pdf")
  )}


### Input function
UPsetInput<- function(){
  tagList( 
    textInput("selected_cluster_upset", "Input Stator state(s):",value = "19,46,47"),
    radioButtons(inputId = "UpsetMode", label = div("Mode:",actionButton("mode_upset",  label = NULL, 
                                                                         icon = icon("info"),  size = "extra-small")),
                 choices = c("Intersect" = "intersect", "Distinct" = "distinct", "Union" = "union"), 
                 selected = "intersect", inline = TRUE),
    
    bsTooltip("mode_upset","Explaination of each mode.",placement = "bottom", trigger = "hover",
              options = NULL), #bsTooltip display a text when the cursor trigger the button (trigger option) The target is the 
    
    actionButton("action_upset","Submit",icon("paper-plane"), 
                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
  )
}




Up_set<-function(selected_cluster,List,mode){
  selected_cluster<-strsplit(selected_cluster, ",\\s*")[[1]]
  names(List)<-gsub("cluster_C:","",names(List))
  #selected_cluster<-paste0("cluster_C:",selected_cluster)
  m=make_comb_mat(List[selected_cluster],mode = mode)  
  return(m)

}
