library(shinydashboard)
library(shiny)
library(shinyBS)

textInputRow<-function (inputId, label, value = "") 
{
    div(style="display:inline-block",
        tags$label(label, `for` = inputId), 
        tags$input(id = inputId, type = "text", class="form-control", value = value))
}

shinyApp(ui=ui, server=server)
