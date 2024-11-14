library(rhandsontable)
library(shiny)
library(shinythemes)

source("model_fns.R")

## User Interface
##
ui <- shinyUI(fluidPage(
  ##
  tags$head(
    tags$style(HTML("
                    div.MathJax_Display{
                    text-align: left !important;
                    }"))),
  ##
  ##
  navbarPage(title = HTML("Urate-Lowering Therapy<br/>PKPD Simulator"), 
             fluid=FALSE,
             theme = shinytheme("flatly"),
             ##
             ## INTRODUCTION MENU
             ##
             tabPanel(h4("Introduction"),
                      tabsetPanel(
                        tabPanel("Background", fluid = TRUE,
                                 ##
                                 mainPanel(
                                   ##
                                   br(),
                                   p("This R Shiny App is based on a publication describing the development and application of a pharmacokinetic-pharmacodynamic (PK-PD) model of febuxostat:",
                                     (tags$em("Hill-McManus D,Soto E,Marshall S,et al. Impact of non-adherence on the safety and efficacy of uric acid-lowering therapies in the treatment of gout. Br J Clin Pharmacol 2018; 84:142-52."))),
                                   p("It uses a two-compartment PK model with first-order absorption and elimination, and a semi-mechanistic PD model of uric acid production and elimination.
                                   The PK model includes covariate submodels on the systemic clearance parameter and inter-individual variation on clearance and absorption."),
                                   p("The PD model comprises four compartments to characterize the time course of serum uric acid (sUA), urinary uric acid (uUA), xanthine and urinary xanthine. 
                                   It utilizes a zero-order production rate (",HTML(paste0("k",tags$sub("0"))),") governing the formation of xanthine and first-order production rates 
                                       characterizing its biotransformation to UA (",HTML(paste0("k",tags$sub("1"))),") and the elimination of xanthine (",HTML(paste0("k",tags$sub("2"))),") and UA (",HTML(paste0("k",tags$sub("3"))),") into the urine. These, in turn, are parameterized in terms of volumes and 
                                       clearance terms.
                                     There is inter-individual variation on the parameters governing the inhibition of the rate of conversion of hypoxanthine to xanthine and xanthine to uric acid."),
                                   p("Under the 'Model Inputs' tab there are sub-pages where the parameters defining the initial state of the subject(s), the PK model and PD model can be modified. 
                                   These parameters are initially all set to the default taken from the paper cited above. 
                                     Under the 'Simulation' tab there are options to perform dose taking simulations for either a single subject or multiple subjects. 
                                     For a single subject, there is no inter-individual variation and all parameters are fixed at the population typical values, 
                                     while for multiple subjects new parameter values are generated for each subject. "),
                                   p("There are additional inputs that can be entered on the pages for single and multiple subject simulations. 
                                   These include the number of subjects, the simulation time, as well as inputs relating to drug adherence. 
                                   Drug adherence includes the timing of treatment dropout and the fraction of dose implementation. 
                                   For the single subject a fixed dropout time can be used, while for multiple subjects the dropout times are randomly generated from a Weibull distribution for which parameters can be inputted. 
                                   The dose implementation fraction specifies the probability of any individual dose being taken, irrespective of treatment dropout. 
                                     Each dose event is treated independently of whether any other doses were taken. "),
                                   br(),br()
                                 )))),
             ##
             ## SET INPUTS MENU
             ##
             tabPanel(h4("Model Inputs"),
                      tabsetPanel(
                        #
                        # SUBJECT INPUTS
                        #
                        tabPanel("Subject", fluid = TRUE,
                                 ##
                                 mainPanel(
                                   uiOutput('ex0'),
                                   h3("Specify Subject Parameters"),
                                   p("The attributes of the subject, or cohort of subjects, can be modified."),
                                   p("Simulation of a single subject will use the mean values of these attributes, whereas a simulation of multiple subjects will generate new values according to the inter-individual variability parameters (IIV)."),
                                   p("These attributes impact on drug PK via covariate submodels (e.g. creatinine clearance and weight)."),
                                   p("They are also used to define the inital state of the PD system (e.g. mean baseline serum uric acid and xanthine)."),
                                   rHandsontableOutput("hot_subject"),
                                   width=12
                                 )),
                        #
                        # PK INPUTS
                        #
                        tabPanel("PK Model", fluid = TRUE,
                                 ##
                                 mainPanel(
                                   h3("PK Model Diagram"),
                                   p("The diagram below presents the two-compartment PK model with first-order absorption and elimination."),
                                   img(src='pkmodel.png', align = "centre", width = 700),
                                   uiOutput('ex1'),
                                   h3("Specify PK Parameters"),
                                   p("The default PK parameters shown below can be modified."),
                                   rHandsontableOutput("hot_pk"),
                                   width=12,
                                   br(),br(),br()
                                 )),
                        #
                        # PD INPUTS
                        #
                        tabPanel("PD Model", fluid = TRUE,
                                 ##
                                 mainPanel(
                                   h3("PD Model Diagram"),
                                   p("The diagram below presents the semi-mechanistic PD model"),
                                   img(src='pdmodel.png', align = "centre", width = 700),
                                   uiOutput('ex2'),
                                   h3("Specify PD Parameters"),
                                   p("The default PD parameters shown below can be modified."),
                                   rHandsontableOutput("hot_pd"),
                                   width=12,
                                   br(),br(),br()
                                 )))),
             ##
             ## SIMULATIONS MENU
             ##
             tabPanel(h4("Simulation"),
                      tabsetPanel(
                        ##
                        ##
                        tabPanel("Simulate One", fluid = FALSE,
                                 
                                 sidebarLayout(
                                   sidebarPanel(width=4,
                                     ##
                                     h3("Set Dose (mg)"),
                                     numericInput("dose", NULL, 80, min=0, max=500),
                                     ##
                                     h3("Adherence Inputs"),
                                     numericInput("imp_lvl", "Implementation Proportion", 1, min=0, max=1),
                                     numericInput("n_do", "Day of dropout", 100, min=0, max=1000),
                                     selectInput("timing_yn", "Variable timing?",
                                                 c("Yes" = TRUE,
                                                   "No" = FALSE)),
                                     h3("Run Simulation"),
                                     numericInput("ndays", "Number of days", 100, min=0, max=1000),
                                     actionButton("go1", "Simulate Subject")
                                   ),
                                   
                                   mainPanel(
                                     h3("Simulation Results"),
                                     plotOutput(outputId = "adhere1", height = 300),
                                     plotOutput(outputId = "sim1", height = 500),
                                     textOutput("text1"),
                                     textOutput("text2"),
                                     textOutput("text3")
                                   )
                                 )
                        ),
                        ##
                        ##
                        tabPanel("Simulate Many", fluid = FALSE,
                                 
                                 sidebarLayout(
                                   sidebarPanel(width=4,
                                     ##
                                     h3("Set Dose (mg)"),
                                     numericInput("dose2", NULL, 80, min=0, max=500),
                                     ##
                                     h3("Adherence Inputs"),
                                     numericInput("imp_lvl2", "Implementation fraction", 1, min=0, max=1),
                                     numericInput("do_scale", "Weibull scale (dropout)", 3.20E-03, min=0, max=1),
                                     numericInput("do_shape", "Weibull shape (dropout)", 0.8, min=0, max=1),
                                     selectInput("timing_yn2", "Variable timing?",
                                                 c("Yes" = TRUE,
                                                   "No" = FALSE)),
                                     h3("Run Simulation"),
                                     numericInput("ndays2", "Number of days", 30, min=1, max=100),
                                     numericInput("nsubj", "Number of subjects", 2, min=1, max=10),
                                     actionButton("go2", "Simulate Subjects")
                                   ),
                                   
                                   mainPanel(
                                     h3("Simulation Results"),
                                     textOutput("text_2a"),
                                     plotOutput("sim2", width = 500),
                                     br(),
                                     textOutput("text_2b"),
                                     plotOutput("hist", width = 500),
                                     textOutput("text_2c"),
                                     textOutput("text4"),br(),br(),br()
                                   )
                                 )
                        )
                      ))
             
  )))

##
##
server <- shinyServer(function(input, output) {
  
  ## Load PK PD defaults
  df_pk <- read.csv("defaults/pk_inputs.csv")
  df_pd <- read.csv("defaults/pd_inputs.csv")
  df_subj <- read.csv("defaults/subj_inputs.csv")
  
  ##
  output$ex1 <- renderUI({
    withMathJax(
      helpText('$$Q = k_{12}V_{2} = k_{21}V_{3}$$'),
      helpText('$$CL = k_{e}V_{2}$$'),
      helpText("$$CL_{i} = CL_{\\mu}*(CRCL_{i}/CRCL_{ref})^{\\tau_{CRCL}}*(WT_{i}/WT_{ref})^{\\tau_{WT}}*e^{\\eta_{CL_{i}}}$$"),
      helpText('$$VC_{i} = VC_{\\mu}*(WT_{i}/WT_{ref})$$'),
      helpText('$$KA_{i} = KA_{\\mu}*e^{\\eta_{KA_{i}}}$$'),
      helpText('$$Where:\\quad \\eta_{X_i} \\sim N(0, \\omega_X)$$')
    )
  })
  
  ##
  output$ex2 <- renderUI({
    withMathJax(
      helpText('$$INH_1 = 1 - \\frac{IMAX*C(t)}{IC50_{1} + C(t)}$$'),
      helpText('$$INH_2 = 1 - \\frac{IMAX*C(t)}{IC50_{1} + C(t)}$$'),
      helpText("$$STIM_1 = 1 + \\frac{EMAX*C(t)}{EC50 + C(t)}$$")
    )
  })
  
  ## Set reactive values
  values <- reactiveValues(doPlot = FALSE)
  
  ## Observe for HOT SUBJECT
  observe({
    if (!is.null(input$hot_subject)) {
      df_subj = hot_to_r(input$hot_subject)
    } else {
      if (is.null(values[["df_subj"]]))
        df_subj <- df_subj
      else
        df_subj <- values[["df_subj"]]
    }
    values[["df_subj"]] <- df_subj
  })
  
  ## Observe for HOT PK
  observe({
    if (!is.null(input$hot_pk)) {
      df_pk = hot_to_r(input$hot_pk)
    } else {
      if (is.null(values[["df_pk"]]))
        df_pk <- df_pk
      else
        df_pk <- values[["df_pk"]]
    }
    values[["df_pk"]] <- df_pk
  })
  
  ## Observe for HOT PD
  observe({
    if (!is.null(input$hot_pd)) {
      df_pd = hot_to_r(input$hot_pd)
    } else {
      if (is.null(values[["df_pd"]]))
        df_pd <- df_pd
      else
        df_pd <- values[["df_pd"]]
    }
    values[["df_pd"]] <- df_pd
  })
  
  output$hot_subject <- renderRHandsontable({
    ##
    df_subj <- values[["df_subj"]]
    if (!is.null(df_subj))
    ##
    col_highlight = 1
    
    rhandsontable(df_subj, col_highlight = col_highlight) %>%
      hot_cols(colWidths = c(300, 100, 100, 100)) %>%
      hot_col("Parameter", readOnly = TRUE) %>%
      hot_col("Unit", readOnly = TRUE) %>%
      hot_col("Description", readOnly = TRUE) %>%
      hot_cols(renderer = "
            function(instance, td, row, col, prop, value, cellProperties) {
                Handsontable.renderers.NumericRenderer.apply(this, arguments);
                if (instance.params) {
                    hcols = instance.params.col_highlight
                    hcols = hcols instanceof Array ? hcols : [hcols]
                }
                if (instance.params && hcols.includes(col)) td.style.background = '#CCCCFF';
            }")
  })
  
  ## Create table for PK
  output$hot_pk <- renderRHandsontable({
    df_pk <- values[["df_pk"]]
    if (!is.null(df_pk))
      ##
      col_highlight = 1
      rhandsontable(df_pk, col_highlight = col_highlight) %>%
      hot_cols(colWidths = c(300, 100, 100, 100)) %>%
      hot_col("Parameter", readOnly = TRUE) %>%
      hot_col("Unit", readOnly = TRUE) %>%
      hot_col("Description", readOnly = TRUE) %>%
      hot_cols(renderer = "
            function(instance, td, row, col, prop, value, cellProperties) {
                Handsontable.renderers.NumericRenderer.apply(this, arguments);
                if (instance.params) {
                    hcols = instance.params.col_highlight
                    hcols = hcols instanceof Array ? hcols : [hcols]
                }
                if (instance.params && hcols.includes(col)) td.style.background = '#CCCCFF';
            }")
  })
  
  ## Create table for PD
  output$hot_pd <- renderRHandsontable({
    df_pd <- values[["df_pd"]]
    if (!is.null(df_pd))
      ##
      col_highlight = 1
      rhandsontable(df_pd, col_highlight = col_highlight) %>%
      hot_cols(colWidths = c(300, 100, 100, 100)) %>%
      hot_col("Parameter", readOnly = TRUE) %>%
      hot_col("Unit", readOnly = TRUE) %>%
      hot_cols(renderer = "
            function(instance, td, row, col, prop, value, cellProperties) {
                Handsontable.renderers.NumericRenderer.apply(this, arguments);
                if (instance.params) {
                    hcols = instance.params.col_highlight
                    hcols = hcols instanceof Array ? hcols : [hcols]
                }
                if (instance.params && hcols.includes(col)) td.style.background = '#CCCCFF';
            }")
  })
  
  ##
  observeEvent(input$go1, {
    
    ## Get inputs values
    ipk <- isolate(values[["df_pk"]])
    ipd <- isolate(values[["df_pd"]])
    isubj <- isolate(values[["df_subj"]])
    
    ## Run simulations
    simdata <- fn_runsim1(nsubj = 1,
                          ndays = isolate(input$ndays),
                          isolate(input$dose),
                          isolate(input$imp_lvl),
                          isolate(input$n_do),
                          isolate(input$do_scale),
                          isolate(input$do_shape),
                          isolate(input$timing_yn),
                          isubj,
                          ipk,
                          ipd,
                          iiv = 0)
    
    ## Print CP plot
    output$adhere1 <-  renderPlot({
      fn_plot_adherence(simdata,
                    isolate(input$ndays))})
    
    ## Print CP plot
    output$sim1 <-  renderPlot({
      fn_plot_pkpd1(simdata,
                    isolate(input$ndays))})
    
    ##
    output$text1 <- renderText(paste(c("Initial sUA conc: ", simdata$sua[1])," mg/dL", sep=""))
    output$text2 <- renderText(paste(c("Final sUA conc: ", round(simdata$sua[length(simdata$sua)],2), "mg/dL"), sep=""))
    output$text3 <- renderText(paste(c("Mean sUA conc (after 14 days): ", round(fn_meansUA(simdata,14), 2), "mg/dL"), sep=""))
  })
  
  ##
  observeEvent(input$go2, {
    
    ## Get inputs values
    ipk <- isolate(values[["df_pk"]])
    ipd <- isolate(values[["df_pd"]])
    isubj <- isolate(values[["df_subj"]])
    
    ## Run simulations
    simdata <- fn_runsim1(nsubj = isolate(input$nsubj),
                          ndays = isolate(input$ndays2),
                          isolate(input$dose2),
                          isolate(input$imp_lvl2),
                          isolate(input$n_do2),
                          isolate(input$do_scale),
                          isolate(input$do_shape),
                          isolate(input$timing_yn2),
                          isubj,
                          ipk,
                          ipd,
                          iiv=1)
    
    ##
    output$text_2a <- renderText("The time course of plasma concentration and sUA concentration for all subjects.")
    
    ## Print CP plot
    output$sim2 <-  renderPlot({
      fn_plot_pkpd2(simdata,
                    isolate(input$ndays2))})
    
    ##
    output$text_2b <- renderText("The distirbution of final sUA concentration.")
    
    ##
    output$hist <-  renderPlot({
      fn_plot_hist(simdata,
                    isolate(input$ndays2))})
    
    ##
    output$text_2c <- renderText("Results summary.")
    output$text4 <- renderText(paste(c("% responding (final sUA < 6 mg/dL): ", round(fn_respondersUA(simdata,6,isolate(input$nsubj))*100, 1), "%"), sep=""))
    
  })
})

## Run the application 
shinyApp(ui = ui, server = server)

##
##