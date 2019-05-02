
 ###########################################################################################################

#' Enrichment visualisation
#'
#' @param x Dataframe of enrichment data
#'
#' @return Shiny app.
#' @examples
#' \dontrun{
#' x <- perfFishTestHPOSingle("Iodine")
#' visEnrich(x)}
#'
#'
#'
#'
#' @export
 visEnrich <- function(x){

if("HPO_Name" %in% colnames(x))
  {
   shinyApp(
     ui = fluidPage(
       ###############################
       tags$style("
                  h1{
                  color:#008080;
                  }"),

       ################################################################
       navbarPage(h1("Phexpo"),
           tabPanel(h3("Dataset"),
                    h4(tags$b("Table of Dataset")),
                    br(),
                    h4("Table 1. Enriched HPO Terms"),
                    div(dataTableOutput('table'),style="font-size:110%")),

           tabPanel(h3("Barchart"),
                  h4(tags$b("Phenotype Enrichment Barchart")),
                  sidebarLayout(
                    sidebarPanel(
                    br(),
                    selectInput("valu","P-value selection",names(x[8:10]), selected = names(x[8])),
                    br()),

                    mainPanel(
                    sliderInput("slider1", label = h4("Threshold Slider"), min = -5, max = 2, value = -3, step = 1),
                    plotOutput("plot"),
                    br(),
                    numericInput("threshold",label = h4("Threshold Value"),0.05,min=0.000000000000000000000000005,max=0.5),
                    plotOutput("plot2")
                    )
                   )

         ))),

     server = function(input, output,session) {

       Dataset <- x
       Dataset_bar<- x %>%
         dplyr::mutate(logp = log10(.data$p_value)) %>%
         dplyr::mutate(logb = log10(.data$bonf)) %>%
         dplyr::mutate(logf = log10(.data$FDR))

       Dataset_2 <- reactive({ if (input$valu == "p_value") {
         Dataset_bar$HPO_Name <- factor(Dataset$HPO_Name, levels = Dataset$HPO_Name[order(-Dataset$p_value)])
         Dataset_bar[,c("HPO_Name","logp")] %>% dplyr::filter(.[[2]] < input$slider1)
       } else if (input$valu == "bonf") {
         Dataset_bar$HPO_Name <- factor(Dataset$HPO_Name, levels = Dataset$HPO_Name[order(-Dataset$bonf)])
         Dataset_bar[,c("HPO_Name","logb")] %>% dplyr::filter(.[[2]] < input$slider1)
       } else {
         Dataset_bar$HPO_Name <- factor(Dataset$HPO_Name, levels = Dataset$HPO_Name[order(-Dataset$FDR)])
         Dataset_bar[,c("HPO_Name","logf")] %>% dplyr::filter(.[[2]] < input$slider1)
       }
                })




       ###### Data Table ##########
       output$table <- renderDataTable({Dataset})  #include.rownames = TRUE, spacing = "xs")

       ######### Plots #########


       output$plot <- renderPlot({
         Dataset_2() %>%
           ggplot2::ggplot(ggplot2::aes(x=.[[1]],y=.[[2]])) +
           ggplot2::geom_bar(stat="identity", fill="springgreen3") +
           ggplot2::theme_bw() +
           ggplot2::theme(text = ggplot2::element_text(size=16)) +
           ggplot2::coord_flip() +
           ggplot2::ggtitle("Phenotype Enrichment") +
           ggplot2::xlab("HPO Term") +
           ggplot2::ylab(input$valu) })



         Dataset_3 <- reactive({
           if (input$valu == "p_value") {
           Dataset_bar$HPO_Name <- factor(Dataset$HPO_Name, levels = Dataset$HPO_Name[order(-Dataset$p_value)])
           Dataset_bar[,c("HPO_Name",input$valu)] %>% dplyr::filter(.[[2]] < input$threshold)
         } else if (input$valu == "bonf") {
           Dataset_bar$HPO_Name <- factor(Dataset$HPO_Name, levels = Dataset$HPO_Name[order(-Dataset$bonf)])
           Dataset_bar[,c("HPO_Name",input$valu)] %>% dplyr::filter(.[[2]] < input$threshold)
         } else {
           Dataset_bar$HPO_Name <- factor(Dataset$HPO_Name, levels = Dataset$HPO_Name[order(-Dataset$FDR)])
         Dataset_bar[,c("HPO_Name",input$valu)] %>% dplyr::filter(.[[2]] < input$threshold)
         }
         })


         output$plot2 <- renderPlot({
           Dataset_3() %>%
             ggplot2::ggplot(ggplot2::aes(x=.[[1]],y=.[[2]])) +
             ggplot2::geom_bar(stat="identity", fill="deepskyblue") +
             ggplot2::theme_bw() +
             ggplot2::theme(text = ggplot2::element_text(size=16)) +
             ggplot2::coord_flip() +
             ggplot2::ggtitle("Phenotype Enrichment") +
             ggplot2::xlab("HPO Term") +
             ggplot2::ylab(input$valu) })







     }
   )


} else {

  shinyApp(
    ui = fluidPage(

#########################################################################################
  tags$style("
             h1{
             color:#008080;
             }"),

  ################################################################
  navbarPage(h1("Phexpo"),
             tabPanel(h3("Dataset"),
                      h4(tags$b("Table of Dataset")),
                      br(),
                      h4("Table 1. Enriched Chemicals"),
                      div(dataTableOutput('table'),style="font-size:110%")),

             tabPanel(h3("Barchart"),
                      h4(tags$b("Chemical Enrichment Barchart")),
                      sidebarLayout(
                        sidebarPanel(
                          br(),
                          selectInput("valu","P-value selection",names(x[8:10]), selected = names(x[8])),
                          br()),

                        mainPanel(
                          sliderInput("slider1", label = h4("Threshold Slider"), min = -5, max = 2, value = -3, step = 1),
                          plotOutput("plot"),
                          br(),
                          numericInput("threshold",label = h4("Threshold Value"),0.05,min=0.000000000000000000000000005,max=0.5),
                          plotOutput("plot2")
                        )
                      )

             ))),

      server = function(input, output) {

        Dataset <- x
        Dataset_bar<- x %>%
          dplyr::mutate(logp = log10(.data$p_value)) %>%
          dplyr::mutate(logb = log10(.data$bonf)) %>%
          dplyr::mutate(logf = log10(.data$FDR))

        Dataset_2 <- reactive({ if (input$valu == "p_value") {
          Dataset_bar$X..ChemicalName <- factor(Dataset$X..ChemicalName, levels = Dataset$X..ChemicalName[order(-Dataset$p_value)])
          Dataset_bar[,c("X..ChemicalName","logp")] %>% dplyr::filter(.[[2]] < input$slider1)
        } else if (input$valu == "bonf") {
          Dataset_bar$X..ChemicalName <- factor(Dataset$X..ChemicalName, levels = Dataset$X..ChemicalName[order(-Dataset$bonf)])
          Dataset_bar[,c("X..ChemicalName","logb")] %>% dplyr::filter(.[[2]] < input$slider1)
        } else {
          Dataset_bar$X..ChemicalName <- factor(Dataset$X..ChemicalName, levels = Dataset$X..ChemicalName[order(-Dataset$FDR)])
          Dataset_bar[,c("X..ChemicalName","logf")] %>% dplyr::filter(.[[2]] < input$slider1)
        }
        })




        ###### Data Table ##########
        output$table <- renderDataTable({Dataset})  #include.rownames = TRUE, spacing = "xs")

        ######### Plots #########


        output$plot <- renderPlot({
          Dataset_2() %>%
            ggplot2::ggplot(ggplot2::aes(x=.[[1]],y=.[[2]])) +
            ggplot2::geom_bar(stat="identity", fill="springgreen3") +
            ggplot2::theme_bw() +
            ggplot2::theme(text = ggplot2::element_text(size=16)) +
            ggplot2::coord_flip() +
            ggplot2::ggtitle("Chemical Enrichment") +
            ggplot2::xlab("Chemical") +
            ggplot2::ylab(input$valu) })



        Dataset_3 <- reactive({
          if (input$valu == "p_value") {
            Dataset_bar$X..ChemicalName <- factor(Dataset$X..ChemicalName, levels = Dataset$X..ChemicalName[order(-Dataset$p_value)])
            Dataset_bar[,c("X..ChemicalName",input$valu)] %>% dplyr::filter(.[[2]] < input$threshold)
          } else if (input$valu == "bonf") {
            Dataset_bar$X..ChemicalName <- factor(Dataset$X..ChemicalName, levels = Dataset$X..ChemicalName[order(-Dataset$bonf)])
            Dataset_bar[,c("X..ChemicalName",input$valu)] %>% dplyr::filter(.[[2]] < input$threshold)
          } else {
            Dataset_bar$X..ChemicalName <- factor(Dataset$X..ChemicalName, levels = Dataset$X..ChemicalName[order(-Dataset$FDR)])
            Dataset_bar[,c("X..ChemicalName",input$valu)] %>% dplyr::filter(.[[2]] < input$threshold)
          }
        })


        output$plot2 <- renderPlot({
          Dataset_3() %>%
            ggplot2::ggplot(ggplot2::aes(x=.[[1]],y=.[[2]])) +
            ggplot2::geom_bar(stat="identity", fill="deepskyblue") +
            ggplot2::theme_bw() +
            ggplot2::theme(text = ggplot2::element_text(size=16)) +
            ggplot2::coord_flip() +
            ggplot2::ggtitle("Chemical Enrichment") +
            ggplot2::xlab("Chemical") +
            ggplot2::ylab(input$valu) })

    }
  )


}
 }
