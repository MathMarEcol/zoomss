library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Old Faithful Geyser Data"),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput("shift",
                  "shift of 2nd set",
                  min = -1,
                  max = 1,
                  value = 0, 
                  step = 0.1)  # I added a step
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot")
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  # This reactive environment can be accessed using data().
  data <- reactive({
    d1 <- matrix(nrow=100, ncol=512)
    for(i in 1:70){
      d1[i,] <- density(rnorm(1000),from = -3, to = 3)$y
    }
    for(i in 71:100){
      d1[i,] <- density(rnorm(1000, input$shift),from = -3, to = 3)$y
    }
    x <- density(rnorm(1000), from = -3, to = 3)$x
    list(x = x, d1 = d1)  # make sure that all objects are returned
  })
  
  output$distPlot <- renderPlot({
    matplot(data()$x, t(data()$d1), type = "l", lty = 1, col = c(rep(1,70),rep(2,30)))
  })
}

# Run the application
shinyApp(ui = ui, server = server)