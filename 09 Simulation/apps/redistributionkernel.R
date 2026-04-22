library(shiny)
library(shinydashboard)
library(amt)
library(circular)
library(sf)
library(terra)

ui <- dashboardPage(
  
  dashboardHeader(title = "Redistribution kernel"),
  
  dashboardSidebar(
    
    actionButton("reset", "Reset parameters"),
    
    hr(),
    
    h4("Movement parameters"),
    
    sliderInput("shape","Gamma shape",
                min=0.1,max=50,value=20),
    
    sliderInput("scale","Gamma scale",
                min=0.1,max=50,value=2),
    
    sliderInput("kappa","Turn angle concentration",
                min=0,max=20,value=1),
    
    hr(),
    
    h4("Selection"),
    
    sliderInput("beta_sel","Selection coefficient",
                min=-10,max=10,value=0,step=0.5),
    
    hr(),
    
    h4("Landscape"),
    
    sliderInput("dist_2","Inner buffer radius",
                min=1,max=99,value=20),
    
    sliderInput("dist_3","Barrier width",
                min=0,max=100,value=10),
    
    hr(),
    
    h4("Starting position"),
    
    sliderInput("start_x","Start x",
                min=-200,max=200,value=0),
    
    sliderInput("start_y","Start y",
                min=-200,max=200,value=0),
    
    sliderInput("ta_start","Start turning angle",
                min=-pi,max=pi,value=0)
    
  ),
  
  dashboardBody(
    
    fluidRow(
      
      box(plotOutput("step_plot"),
          width=6,
          title="Step-length distribution",
          status="primary"),
      
      box(plotOutput("angle_plot"),
          width=6,
          title="Turn-angle distribution",
          status="primary")
      
    ),
    
    fluidRow(
      
      box(plotOutput("hab_plot",
                     height=400,
                     click="map_click"),
          width=6,
          title="Habitat landscape (click to set start location)",
          status="primary"),
      
      box(plotOutput("rdk_plot",height=400),
          width=6,
          title="Redistribution kernel",
          status="primary")
      
    )
    
  )
)

server <- function(input, output, session) {
  
  model <- reactive({
    
    make_issf_model(
      coefs = c("x_end" = input$beta_sel),
      sl = make_gamma_distr(shape=input$shape, scale=input$scale),
      ta = make_vonmises_distr(kappa=input$kappa)
    )
    
  })
  
  landscape <- reactive({
    
    r1 <- rast(xmin=-200, xmax=200,
               ymin=-200, ymax=200, res=1)
    
    p1 <- st_buffer(st_point(c(0,0)), dist=100)
    p2 <- st_buffer(st_point(c(0,0)), dist=input$dist_2)
    
    l1 <- st_linestring(cbind(c(0,0),c(-200,200))) |>
      st_buffer(dist=input$dist_3)
    
    hab <- st_as_sfc(
      list(st_difference(p1,p2) |> st_difference(l1))
    ) |> st_as_sf()
    
    r1 <- rasterize(hab, r1, background=0)
    names(r1) <- "x"
    
    r1
    
  })
  
  start_loc <- reactive({
    make_start(c(input$start_x,input$start_y), ta_ = input$ta_start)
  })
  
  rdk <- reactive({
    
    redistribution_kernel(
      model(),
      map = landscape(),
      start = start_loc(),
      as.rast = TRUE,
      landscape = "discrete"
    )
    
  })
  
  observeEvent(input$map_click, {
    
    updateSliderInput(session, "start_x",
                      value = round(input$map_click$x))
    
    updateSliderInput(session, "start_y",
                      value = round(input$map_click$y))
    
  })
  
  observeEvent(input$reset, {
    
    updateSliderInput(session,"shape",value=20)
    updateSliderInput(session,"scale",value=2)
    updateSliderInput(session,"kappa",value=5)
    
    updateSliderInput(session,"beta_sel",value=10)
    
    updateSliderInput(session,"dist_2",value=20)
    updateSliderInput(session,"dist_3",value=10)
    
    updateSliderInput(session,"start_x",value=0)
    updateSliderInput(session,"start_y",value=0)
    
    updateSliderInput(session,"ta_start",value=0)
    
  })
  
  output$step_plot <- renderPlot({
    
    x <- seq(
      1,
      round(qgamma(0.99,
                   shape=input$shape,
                   scale=input$scale)),
      len=100
    )
    
    plot(
      x,
      dgamma(x,
             shape=input$shape,
             scale=input$scale),
      type="l",
      xlab="distance",
      ylab="density",
      col="steelblue",
      lwd=2
    )
    
  })
  
  output$angle_plot <- renderPlot({
    
    theta <- seq(-pi,pi,len=100)
    
    plot(
      theta,
      circular::dvonmises(theta,
                          mu=0,
                          kappa=input$kappa),
      type="l",
      xlab="relative turn angle",
      ylab="density",
      col="darkgreen",
      lwd=2
    )
    
  })
  
  output$hab_plot <- renderPlot({
    
    r <- landscape()
    
    plot(r)
    
    points(
      input$start_x,
      input$start_y,
      pch=21,
      bg="red",
      col="white",
      cex=2
    )
    
  })
  
  output$rdk_plot <- renderPlot({
    
    k <- rdk()$redistribution.kernel
    
    terra::plot(k)
    
    points(
      input$start_x,
      input$start_y,
      pch=21,
      bg="red",
      col="white",
      cex=2
    )
    
  })
  
}

shinyApp(ui, server)