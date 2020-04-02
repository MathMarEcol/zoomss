## app.R ##
library(shinydashboard)
library(tidyverse)

header <- dashboardHeader(title = "ZooMSS Dashboard")
# dat <- readRDS("../RawOutput/DATE_JOBNAME_0350.RDS")

## Sidebar content
sidebar <- dashboardSidebar(
    sidebarMenu(
        menuItem("Overview", tabName = "dashboard", icon = icon("dashboard")),
        menuItem("Temporal Charts",  tabName = "charts", icon = icon("bar-chart-o")),
        menuItem("Widgets", tabName = "widgets", icon = icon("th"))
    )
)

## Body content
body <- dashboardBody(
    tabItems(
        # First tab content
        tabItem(tabName = "dashboard",
                h2("ZooMSS Final Static Output"),

                fluidRow(
                    valueBoxOutput("SST", width = 3),
                    valueBoxOutput("Chl", width = 3),
                    valueBoxOutput("Lat", width = 3),
                    valueBoxOutput("Lon", width = 3)),

                fluidRow(
                    fileInput("ModelFile", "Choose RDS File", accept = c(".RDS"), width = "80%", ),
                    tableOutput('GroupsTable')
                ),
                br(),
                fluidRow(
                    column(6, plotOutput("speciesPlot", height = 400)),
                    column(6, plotOutput("PPMRPlot", height = 400))
                ),
                br(),
                fluidRow(
                    column(6, plotOutput("growthPlot", height = 400)),
                    column(6, plotOutput("predationPlot", height = 400))
                )
        ),

        # Second tab content
        tabItem(tabName = "charts",
                h2("Time varying plots"),
                fluidRow(box(plotOutput("timeSpeciesPlot"), width = "100%")),
                fluidRow(box(plotOutput("timeGrowthPlot"), width = "100%")),
                fluidRow(box(plotOutput("timePredationPlot"), width = "100%"))
        ),

        # Third tab content
        tabItem(tabName = "widgets",
                h2("Widgets tab content")
        )
    )
)

ui <- dashboardPage(header, sidebar, body)

server <- function(input, output) {
    options(shiny.maxRequestSize=100*1024^2) # Set max file size to 100 MB

    dat <- reactive({
        # browser()
        # if (is.null(input$ModelFile))
        #     return(NULL)
        req(input$ModelFile)
        readRDS(input$ModelFile$datapath)
    })

    output$copepods <- renderValueBox({
        # if (is.null(dat()))
        #     return(NULL)
        NoCopepods <- round(sum(dat()$abundances[dat()$model$param$Groups$Species=="OmniCopepods" | dat()$model$param$Groups$Species=="CarnCopepods"]))
        valueBox(
            value = formatC(NoCopepods, digits = 1, format = "f"),
            subtitle = "Total Number of copepods (m-3)",
            icon = icon("area-chart"),
            color = "yellow"
        )
    })

    output$fish <- renderValueBox({
        # if (is.null(dat()))
        #     return(NULL)
        fish <- colSums(dat()$abundances[which(dat()$model$param$Groups$Type=="Fish"),])
        FishBiomass <- sum(fish * dat()$model$w)

        valueBox(
            value = formatC(FishBiomass, digits = 1, format = "f"),
            subtitle = "Mean Biomass of Fish (m-3)",
            icon = icon("area-chart"),
            color = "red"
        )
    })

    output$SST <- renderValueBox({
        # if (is.null(dat()))
        #     return(NULL)
        valueBox(
            value = formatC(dat()$model$param$sst, digits = 2, format = "f"),
            subtitle = "SST",
            # icon = icon("area-chart"),
            color = "red"
        )
    })

    output$Chl <- renderValueBox({
        # if (is.null(dat()))
        #     return(NULL)
        valueBox(
            value = formatC(dat()$model$param$chlo, digits = 2, format = "f"),
            subtitle = "Chl. a (mg m-3)",
            # icon = icon("area-chart"),
            color = "green"
        )
    })

    output$Lat <- renderValueBox({
        # if (is.null(dat()))
        #     return(NULL)
        valueBox(
            value = formatC(dat()$model$param$Lat, digits = 2, format = "f"),
            subtitle = "Latitude",
            # icon = icon("area-chart"),
            color = "blue"
        )
    })
    output$Lon <- renderValueBox({
        # if (is.null(dat()))
        #     return(NULL)
        valueBox(
            value = formatC(dat()$model$param$Lon, digits = 2, format = "f"),
            subtitle = "Longitude",
            # icon = icon("area-chart"),
            color = "blue"
        )
    })


    output$GroupsTable <- renderTable({
        # if (is.null(dat()))
        #     return(NULL)
        dat()$model$param$Groups %>%
            dplyr::select(-c(Repro, PlotColour, GrossGEscale))},
        striped = TRUE, bordered = TRUE, label = "Test Group Parameters")


    # Plot abundance spectra by species
    output$speciesPlot <- renderPlot({
        # if (is.null(dat()))
        #     return(NULL)
        species <- dat()$abundances
        rownames(species) <- dat()$model$param$Groups$Species
        species <- as_tibble(t(species))
        species <- species %>%
            add_column("Weight" = dat()$model$w)
        species <- species %>%
            pivot_longer(-Weight, names_to = "Species", values_to = "Abundance") %>%
            filter(Abundance > 0) %>%
            mutate(Species = factor(Species, levels = dat()$model$param$Groups$Species))

        ggplot(data = species, mapping = aes(x = log10(Weight), y = log10(Abundance), colour = Species)) +
            geom_line() +
            geom_point() +
            scale_color_manual(values = dat()$model$param$Groups$PlotColour) +
            theme_bw() +
            labs(subtitle = "Abundance Spectrum")
    })

    # Plot growth rates by species
    output$growthPlot <- renderPlot({
        if (is.null(dat()))
            return(NULL)
        growth <- dat()$growth
        rownames(growth) <- dat()$model$param$Groups$Species
        growth <- as_tibble(t(growth))
        growth <- growth %>%
            add_column("Weight" = dat()$model$w)
        growth <- growth %>%
            pivot_longer(-Weight, names_to = "Species", values_to = "GrowthRate") %>%
            filter(GrowthRate > 0) %>%
            mutate(Species = factor(Species, levels = dat()$model$param$Groups$Species))

        ggplot(data = growth, mapping = aes(x = log10(Weight), y = log10(GrowthRate), colour = Species)) +
            geom_line() +
            geom_point() +
            scale_color_manual(values = dat()$model$param$Groups$PlotColour) +
            theme_bw() +
            labs(subtitle = "Growth Rates")
    })

    # Plot predation by species
    output$predationPlot <- renderPlot({
        if (is.null(dat()))
            return(NULL)
        predation <- dat()$predation
        rownames(predation) <- dat()$model$param$Groups$Species
        predation <- as_tibble(t(predation))
        predation <- predation %>%
            add_column("Weight" = dat()$model$w)
        predation <- predation %>%
            pivot_longer(-Weight, names_to = "Species", values_to = "PredationRate") %>%
            mutate(Species = factor(Species, levels = dat()$model$param$Groups$Species))

        ggplot(data = predation, mapping = aes(x = log10(Weight), y = log10(PredationRate), colour = Species)) +
            geom_line() +
            geom_point() +
            scale_color_manual(values = dat()$model$param$Groups$PlotColour) +
            theme_bw() +
            labs(subtitle = "Predation Rates")
    })


    # Plot predation by species
    output$PPMRPlot <- renderPlot({
        if (is.null(dat()))
            return(NULL)

        min_size = min(dat()$model$param$Groups$W0) # smallest size class
        max_size = max(dat()$model$param$Groups$Wmax) # largest size class
        # w = 10^(seq(from = min_size, to = max_size, 0.1)) # all size classes

        # Calculate PPMR (beta) table, where dim1 = group, dim2 = body size with
        # value being PPMR for that body size (this is not realised PPMR - not
        # emergent from diet but calculated from m-values and Wirtz, 2012 equation)
        D.z = 2*(3*(dat()$model$w)*1e12/(4*pi))^(1/3) # convert body mass g to ESD (um)
        zoo_m = dat()$model$param$Groups$PPMRscale # pull out m-values from parameter table
        betas =  log10(t(sapply(zoo_m, function(x){(exp(0.02*log(D.z)^2 - x + 1.832))^3}))) # Convert m to betas, using Wirtz 2012 equation
        betas = betas[-which(dat()$model$param$Groups$Type=="Fish"),] # remove fish rows


        ## Modify beta matrix for larvaceans and salps - all size classes for these groups feed on same prey, so log10PPMR increases by 0.1 for each 0.1 log10 size interval
        betas[which(dat()$model$param$Groups$Species=="Larvaceans"),45:75] <- betas[which(dat()$model$param$Groups$Species=="Larvaceans"),44] + seq(0.1,3.1,0.1) # Larvaceans (index 44 in w vector is smallest size class, 75 is maximum size class)
        betas[which(dat()$model$param$Groups$Species=="Salps"),61:121] <- betas[which(dat()$model$param$Groups$Species=="Salps"),61] + seq(0.1,6.1,0.1) # Larvaceans (index 61 in w vector is smallest size class, 121 is maximum size class

        # Calculate ave abundances across oligo/eutro grid squares,
        # then calculate ave biomass and proportion of total zoo biomass
        # that is from each group size class

        ave_biom = sweep(dat()$abundances, 2, dat()$model$w, "*") # Calculate oligo biomass for zoo groups
        ave_biom = ave_biom[-which(dat()$model$param$Groups$Type=="Fish"),] # remove rows for fish
        beta_props = ave_biom/sum(ave_biom) # Calculate fraction of zoo biomass in each group, in each size class

        out <- list()
        out[[1]] <- betas
        out[[2]] <- beta_props
        names(out) <- c("betas", "beta_props")

        temp <- density(betas, weights = beta_props)

        PPMR <- tibble("x" = temp$x, "y" = temp$y)
        PPMR$mn_beta <- sum(beta_props*betas)

        ggplot() +
            geom_line(data = PPMR, mapping = aes(x = x, y = y), size = 1.2) +
            theme_bw() +
            theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
            labs(x = expression('log' [10] * PPMR),
                 y = "Zoop. Biomass Proportion", subtitle = "PPMR") +
            geom_vline(data = PPMR, mapping = aes(xintercept = mn_beta), colour = 'black')

    })




    # Plot abundance by time
    output$timeSpeciesPlot <- renderPlot({
        if (is.null(dat()))
            return(NULL)
        tspecies <- rowSums(dat()$model$N,dims = 2)
        colnames(tspecies) <- dat()$model$param$Groups$Species
        tspecies <- as_tibble(tspecies)
        tspecies$Time <- seq(dat()$model$param$dt*dat()$model$param$isave,dat()$model$param$tmax, dat()$model$param$dt*dat()$model$param$isave)
        tspecies <- tspecies %>%
            pivot_longer(-Time, names_to = "Species", values_to = "Abundance") %>%
            filter(Abundance > 0) %>%
            mutate(Species = factor(Species, levels = dat()$model$param$Groups$Species))

        ggplot(data = tspecies, mapping = aes(x = Time, y = log10(Abundance), colour = Species)) +
            geom_line(size = 0.2) +
            geom_point(size = 0.2) +
            scale_color_manual(values = dat()$model$param$Groups$PlotColour) +
            theme_bw() +
            labs(subtitle = "Abundance") +
            xlab("Time (Years)")
    })


    # Plot abundance by time
    output$timeGrowthPlot <- renderPlot({
        if (is.null(dat()))
            return(NULL)
        gg <- rowSums(dat()$model$gg,dims = 2)
        colnames(gg) <- dat()$model$param$Groups$Species
        gg <- as_tibble(gg)
        gg$Time <- seq(dat()$model$param$dt*dat()$model$param$isave,dat()$model$param$tmax, dat()$model$param$dt*dat()$model$param$isave)
        gg <- gg %>%
            pivot_longer(-Time, names_to = "Species", values_to = "Growth") %>%
            filter(Growth > 0) %>%
            mutate(Species = factor(Species, levels = dat()$model$param$Groups$Species))

        ggplot(data = gg, mapping = aes(x = Time, y = log10(Growth), colour = Species)) +
            geom_line(size = 0.2) +
            geom_point(size = 0.2) +
            scale_color_manual(values = dat()$model$param$Groups$PlotColour) +
            theme_bw() +
            labs(subtitle = "Growth Rate") +
            xlab("Time (Years)")
    })

    # Plot predation by time
    output$timePredationPlot <- renderPlot({
        if (is.null(dat()))
            return(NULL)
        m2 <- rowSums(dat()$model$M2,dims = 2)
        colnames(m2) <- dat()$model$param$Groups$Species
        m2 <- as_tibble(m2)
        m2$Time <- seq(dat()$model$param$dt*dat()$model$param$isave,dat()$model$param$tmax, dat()$model$param$dt*dat()$model$param$isave)
        m2 <- m2 %>%
            pivot_longer(-Time, names_to = "Species", values_to = "Predation") %>%
            filter(Predation > 0) %>%
            mutate(Species = factor(Species, levels = dat()$model$param$Groups$Species))

        ggplot(data = m2, mapping = aes(x = Time, y = Predation, colour = Species)) +
            geom_line(size = 0.2) +
            geom_point(size = 0.2) +
            scale_color_manual(values = dat()$model$param$Groups$PlotColour) +
            theme_bw() +
            labs(subtitle = "Predation Rate") +
            xlab("Time (Years)")
    })

}

shinyApp(ui, server)
