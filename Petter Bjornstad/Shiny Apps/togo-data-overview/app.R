library(shiny)
library(shinymanager)

# 1) credentials
credentials <- data.frame(
  user     = c("shiny", "admin"),
  password = c("azerty", "12345"),
  admin    = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

# 2) your "real" UI (raw)
raw_ui <- fluidPage(
  titlePanel("Protected app"),
  mainPanel(
    textOutput("whoami"),
    textOutput("msg")
  )
)

# 3) wrap UI in secure_app
ui <- secure_app(raw_ui)

# 4) server logic
server <- function(input, output, session) {
  # authentication
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )
  
  # only runs if logged in
  output$whoami <- renderText({
    paste("User:", res_auth$user)
  })
  output$msg <- renderText("It works after login!")
}

shinyApp(ui, server)