library(shiny)
library(DT)
library(survival)
library(shinydashboard)
library(ggplot2)

ui <- dashboardPage(
  dashboardHeader(title = "Survival Analyze"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("File Upload", tabName = "file_upload", icon = icon("file-upload"))
    )
  ),
  dashboardBody(
    tags$style(HTML("
      .flex-container {
        display: flex;
        justify-content: between;
        height: 600px;
      }
      .flex-Input-add {
        width: 1100px;
      }
       .flex-table {
        width: 4100px;
      }
    ")),
    tabItems(
      tabItem(tabName = "file_upload",
              div(class = "flex-container",
                  div(class = "flex-Input-add",
                      box(
                        title = "Input File", status = "primary", solidHeader = TRUE,
                        fileInput("file1", "Choose CSV File",
                                  multiple = FALSE,
                                  accept = c("text/csv",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        checkboxInput("header", "Header", TRUE),
                        radioButtons("disp", "Display",
                                     choices = c(Head = "head",
                                                 All = "all"),
                                     selected = "head")
                      )
                  ),
                  div(class = "flex-Input-add",
                      box(
                        title = "Add New Patient", status = "primary", solidHeader = TRUE,
                        numericInput("new_t_event", "t.event", value = 0),
                        selectInput("new_event", "event", choices = c("0", "1"), selected = "0"),
                        numericInput("new_t_death", "t.death", value = 0),
                        selectInput("new_death", "death", choices = c("0", "1"), selected = "0"),
                        selectInput("new_group", "group", choices = c("4", "8", "11", "14"), selected = "4"),
                        numericInput("new_CXCL12", "CXCL12", value = 0.0),
                        actionButton("add_btn", "Add")
                      )
                  ),
                  div(class = "flex-table",
                      box(
                        title = "Data Table", status = "primary", solidHeader = TRUE,
                        DT::dataTableOutput("editable_table")
                      )
                  )
              ),
              fluidRow(
                box(
                  title = "Cox Model Survival Plot", status = "primary", solidHeader = TRUE,
                  plotOutput("cox_plot")
                ),
                box(
                  title = "Exponential Survival Plot", status = "primary", solidHeader = TRUE,
                  plotOutput("exp_plot")
                )
              ),
              fluidRow(
                box(
                  title = "Relative Risk (RR)", status = "primary", solidHeader = TRUE,
                  textOutput("rr_value")
                ),
                box(
                  title = "P-value", status = "primary", solidHeader = TRUE,
                  textOutput("p_value")
                )
              )
      )
    )
  )
)

server <- function(input, output, session) {
  
  df <- reactiveVal(NULL)
  
  observeEvent(input$file1, {
    req(input$file1)
    tryCatch(
      {
        data <- read.csv(input$file1$datapath, header = input$header)
        data <- na.omit(data)
        data$t.event <- as.numeric(data$t.event)
        data$event <- as.numeric(data$event)
        df(data)
        print(head(data))  
        print(summary(data))  
      },
      error = function(e) {
        stop(safeError(e))
      }
    )
  })
  
  observeEvent(input$add_btn, {
    new_row <- data.frame(
      t.event = as.numeric(input$new_t_event),
      event = as.numeric(input$new_event),
      t.death = as.numeric(input$new_t_death),
      death = as.numeric(input$new_death),
      group = as.numeric(input$new_group),
      CXCL12 = as.numeric(input$new_CXCL12)
    )
    df(rbind(df(), new_row))
    write.csv(df(), input$file1$datapath, row.names = FALSE)
  })
  
  output$editable_table <- DT::renderDataTable({
    DT::datatable(df(), editable = TRUE)
  })
  
  observeEvent(input$editable_table_cell_edit, {
    info <- input$editable_table_cell_edit
    str(info)
    df()[info$row, info$col] <- info$value
    write.csv(df(), input$file1$datapath, row.names = FALSE)
  })
  
  output$cox_plot <- renderPlot({
    req(df())
    data <- df()
    median_CXCL12 <- median(data$CXCL12, na.rm = TRUE)
    data$CXCL12_group <- ifelse(data$CXCL12 > median_CXCL12, "High", "Low")
    res <- survfit(Surv(t.event, event) ~ CXCL12_group, data = data)
    plot(res, col = c("red", "blue"), lty = 1:2, xlab = "Time", ylab = "Survival Probability")
    legend("bottomleft", legend = c("Low CXCL12", "High CXCL12"), col = c("blue", "red"), lty = 1:2)
  })
  
  output$exp_plot <- renderPlot({
    req(df())
    data <- df()
    median_CXCL12 <- median(data$CXCL12, na.rm = TRUE)
    data$CXCL12_group <- ifelse(data$CXCL12 > median_CXCL12, "High", "Low")
    
    high_data <- data[data$CXCL12_group == "High",]
    low_data <- data[data$CXCL12_group == "Low",]
    
    lambda_high <- sum(high_data$event) / sum(high_data$t.event)
    lambda_low <- sum(low_data$event) / sum(low_data$t.event)
    
    print(paste("Lambda High:", lambda_high))
    print(paste("Lambda Low:", lambda_low))
    
    time_points <- seq(0, max(data$t.event), length.out = 100)
    survival_prob_high <- exp(-lambda_high * time_points)
    survival_prob_low <- exp(-lambda_low * time_points)
    
    plot(time_points, survival_prob_high, type = "l", col = "red", xlab = "Time", ylab = "Survival Probability",
         main = "Exponential Survival Plot")
    lines(time_points, survival_prob_low, col = "blue")
    legend("bottomleft", legend = c("Low CXCL12","High CXCL12"), col = c("blue", "red"), lty = 1:1)
  })
  
  output$rr_value <- renderText({
    req(df())
    data <- df()
    
    median_CXCL12 <- median(data$CXCL12, na.rm = TRUE)
    data$CXCL12_group <- ifelse(data$CXCL12 > median_CXCL12, "High", "Low")
    
    high_data <- data[data$CXCL12_group == "High",]
    low_data <- data[data$CXCL12_group == "Low",]
    
    risk_high <- sum(high_data$event) / sum(high_data$t.event)
    risk_low <- sum(low_data$event) / sum(low_data$t.event)
    
    rr <- risk_high / risk_low
    paste("Relative Risk (RR) between High CXCL12 and Low CXCL12 groups: ", round(rr, 2))
  })
  
  output$p_value <- renderText({
    req(df())
    data <- df()
    
    median_CXCL12 <- median(data$CXCL12, na.rm = TRUE)
    data$CXCL12_group <- ifelse(data$CXCL12 > median_CXCL12, "High", "Low")
    
    high_data <- data[data$CXCL12_group == "High",]
    low_data <- data[data$CXCL12_group == "Low",]
    
    risk_high <- sum(high_data$event) / sum(high_data$t.event)
    risk_low <- sum(low_data$event) / sum(low_data$t.event)
    rr <- risk_high / risk_low
    
    log_rr <- log(rr)
    var_log_rr <- 1/sum(high_data$event) + 1/sum(low_data$event) # Variance of log(RR)
    z_score <- abs(log_rr) / sqrt(var_log_rr)
    p_value <- 2 * (1 - pnorm(abs(z_score)))
    
    if(p_value < 0.05)
    {
      paste("P-value: ", formatC(p_value, format = "f", digits = 10), ", H1 : RR != 1 is accepted ")
    }else{
      paste("P-value: ", formatC(p_value, format = "f", digits = 10), ", H0 : RR = 1 is accepted ")
    }
  })
}

shinyApp(ui, server)