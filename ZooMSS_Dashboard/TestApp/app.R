## app.R ##
library(shinydashboard)
library(tidyverse)

header <- dashboardHeader(title = "ZooMSS Dashboard")
# dat <- readRDS("../RawOutput/DATE_JOBNAME_0350.RDS")

## Sidebar content
sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Overview", tabName = "dashboard", icon = icon("dashboard"))# ,
    # menuItem("Temporal Charts",  tabName = "charts", icon = icon("bar-chart-o")),
    # menuItem("Widgets", tabName = "widgets", icon = icon("th"))
  )
)


## Body content
body <- dashboardBody(
  tabItems(
    # First tab content
    tabItem(tabName = "dashboard", h2("ZooMSS Final Static Output"),
            fluidRow(
              valueBoxOutput("SST", width = 3),
              valueBoxOutput("Chl", width = 3),
              valueBoxOutput("Lat", width = 3),
              valueBoxOutput("Lon", width = 3)),
            
            fluidRow(
              fileInput("ModelFile", "Choose RDS File", accept = c(".RDS"), width = "80%", ),
              tableOutput('GroupsTable')
            )#,
            # fluidRow(
            #   box(plotOutput("speciesPlot", height = 400)),
            #   box(plotOutput("PPMRPlot", height = 400)),
            #   # box(title = "Controls", sliderInput("slider", "Number of observations2:", 1, 100, 50))
            # ),
            # fluidRow(
            #   box(plotOutput("growthPlot", height = 400)),
            #   box(plotOutput("predationPlot", height = 400)),
            #   # box(plotOutput("growthPlot", height = 400)),
            #   # box(title = "Controls", sliderInput("slider", "Number of observations2:", 1, 100, 50))
            # )
    )#,
    
    # # Second tab content
    # tabItem(tabName = "charts",
    #         h2("Time varying plots"),
    #         fluidRow(box(plotOutput("timeSpeciesPlot"), width = "100%")),
    #         fluidRow(box(plotOutput("timeGrowthPlot"), width = "100%")),
    #         fluidRow(box(plotOutput("timePredationPlot"), width = "100%"))
    # ),
    
    # # Third tab content
    # tabItem(tabName = "widgets",
    #         h2("Widgets tab content")
    # )
  )
)

ui <- dashboardPage(header, sidebar, body)

server <- function(input, output) {
  
  # dat_out <- readRDS(input$ModelFile$datapath)
  # dat <- reactiveFileReader(1000, session, input$ModelFile$datapath, readRDS)
  dat <- reactive({
    # browser()
    if (is.null(input$ModelFile))
      return(NULL)
    
    readRDS(input$ModelFile$datapath)
    })
  
  
  output$GroupsTable <- renderTable({
    if (is.null(dat()))
      return(NULL)
    dat()$model$param$Groups %>% 
      select(-c(Repro, PlotColour, GrossGEscale))}, 
    striped = TRUE, bordered = TRUE, label = "Test Group Parameters")
  
}

shinyApp(ui, server)
