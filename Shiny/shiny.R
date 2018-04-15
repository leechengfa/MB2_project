library(shiny)

##Background codes (once-off loading)
#libraries
library(raster)
library(rgdal)
library(RStoolbox)
library(lattice)
library(ggplot2)
#dplyr and graphics --> perhaps do a "if require loop"

#data from github --> pre-generated band ratios, spectral indices, map, etc.
ras <-brick("data/p224r63_2011.grd") 
train <- readOGR("data/training_2011.shp")
valid <- readOGR("data/validation_2011.shp")
#User values
#basic superClass --> get benchmark timing and accuracy
ptm <- proc.time() #start timer
sc<- superClass(ras,trainData=train,valData=valid, 
                responseCol="class_name")
norm.timing <-  proc.time() - ptm #end timer
new.timing <- c(0,0,0)
sc.acc <- sc$validation$performance$overall[1]
#Downloading pre-rendered stacks
#band ratio stack
br.stack <- brick("raster/bandratio.tif",geoinfo=geoinfo) #band ratio stack from GitHub
br.max <- length(names(br.stack))
#spectral Indices stack
si.stack <- brick("raster/specras.tif",geoinfo=geoinfo) #spectral indices stack from GitHub
si.max <- length(names(si.stack))
#tasseled cap stack
tc.stack <- brick("raster/tasseled.tif",geoinfo=geoinfo) #tasseled cap stack from GitHub


##User Interface
ui <- fluidPage(
  #App title
  titlePanel("Is increasing the amount of spectral features always better for classification?"),
  
  #Panels
  sidebarLayout(
    #Side panel
    sidebarPanel(
      #Option to select additional amount of spectral bands to add
      helpText("Please choose the amount of additional spectral features to add.
               Layers will be randomly chosen from these options."),
      fluidRow(
        column(10,
               sliderInput("bandratio",label="Band ratio layers",
                           h3("Sliders"),
                           min=0,max=br.max,value=0,step=1), #max value = length of band ratio stack
               sliderInput("specIndices",label="Spectral indices",
                           h3("Sliders"),
                           min=0,max=si.max,value=0,step=1), #max value = length of spectral indices stack
               sliderInput("tasseledcap",label="Tasseled Cap layers",
                           h3("Sliders"),
                           min=0,max=3,value=0,step=1), #max value = length of tasseled cap stack
               #text output
               textOutput("band.selection"),
               br(),
               actionButton("update",label="Run classification"),
               br(),
               textOutput("exe_msg")
        )
      )
    ),
    #Main panel
    mainPanel(
      fluidPage(
        #Output (compare with original)
        h4("Classification parameters for a normal stack"),
        h6(paste("Processing time =",round(norm.timing[3],digits=3),"seconds")),
        h6(paste("Overall accuracy =",round(sc.acc,digits=4)*100,"%")),
        h4("Classification parameters for current stack"),
        h6(textOutput("time.msg")),
        h6(textOutput("acc.msg")),
        h4(paste("Classification map")),
        #Map output (classification map)
        plotOutput("new.map")
      )
    )
  )
)


##Server logic
server <- function(input,output){
  
  #text output for band selection
  output$band.selection <- renderText({
    paste("You have selected a total of",
          input$bandratio+input$specIndices+input$tasseledcap,
          "additional bands."
          )
  })
  
  
  #Creating reactive values
  r.values <- reactiveValues(t=0,a=0)
  map <- reactiveValues()
  c.msg <- reactiveValues(o="Let's get started.")
  
  #Execute band stacking (action button)
 observeEvent(input$update,
              isolate({
                c.msg$o <- "Proceeding to classify."
                run.stack <- ras
                  if(input$bandratio > 0){
                    br.select <- br.stack[[sort(sample(1:br.max,input$bandratio))]]
                    run.stack <- stack(run.stack,br.select)
                    }
                  if(input$specIndices > 0){
                    si.select <-  si.stack[[sort(sample(1:si.max,input$specIndices))]]
                    run.stack <- stack(run.stack,si.select)}
                  if(input$tasseledcap > 0){
                    tc.select <- tc.stack[[sort(sample(1:3,input$tasseledcap))]]
                    run.stack <- stack(run.stack,tc.select)}
                htm <- proc.time() #start timer
                  sc<- superClass(run.stack,trainData=train,valData=valid,
                                  responseCol="class_name")
                  new.timing <-  proc.time() - htm #end timer
                  r.values$t <- new.timing[3]
                  r.values$a <- sc$validation$performance$overall[1]
                  map$m <- sc$map
                c.msg$o <- "Classification complete."
  })
  )

  
  #text output for execution message
  output$exe_msg <- renderText(paste(c.msg$o))
  
   #output timing and overall accuracy
  output$time.msg <- renderText({
    paste("Processing Time =",
          round(r.values$t,digits=3),
          "seconds")
  })
  output$acc.msg <- renderText({
    paste("Overall accuracy =",
          round(r.values$a,digits=4)*100,
          "%")
  })

  #output map
  #reactive for output map here
  output$new.map <- renderPlot({plot(map$m)})
}


##Run application
shinyApp(ui=ui,server=server)
