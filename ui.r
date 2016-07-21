# server.R
library(shinydashboard)
library(shiny)
library(RColorBrewer)
library(ggplot2)
library(shinyBS)
library(plotly)
library(FactoMineR)
library(VennDiagram)


families_list = list( "GABAa Receptors (ionotropic)"="GABAa Receptors (ionotropic)" ,
                      "GABAb Receptors (metabotropic)"="GABAb Receptors (metabotropic)",
                      "Not-NMDA Glutamate Rec. (ionotropic)"="Not-NMDA Glutamate Rec. (ionotropic)",
                      "NMDA Glutamate Rec. (ionotropic)"="NMDA Glutamate Rec. (ionotropic)",
                      "Metabotropic Glutamate Receptors"="Metabotropic Glutamate Receptors",
                      "Potassium-Channels" = "Potassium-Channels", 
                      "Sodium-Channels" = "Sodium-Channels", 
                      "Calcium-Channels"="Calcium-Channels","CNG-Channels"="CNG-Channels" ,
                      "Ion channels and receptors" = "Ion channels and receptors",
                      "Synaptic release machinery genes" = "Synaptic release machinery genes")

DE_list<-c(c("All genes"="All genes", "Specific genes"="Specific genes"),families_list )

shinyUI(
  dashboardPage(
    dashboardHeader(title=""),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Upload", tabName = "upload", icon = icon("upload")),
        menuItem("Gene Families", tabName = "Families", icon = icon("pie-chart")),
        menuItem("Individual Genes", tabName = "Indiv", icon = icon("bar-chart")),
        menuItem("PCA", tabName = "PCAtab", icon = icon("mouse-pointer")),
        menuItem("Differential Expression Analysis", tabName = "DEA", icon = icon("area-chart"))
        
        #menuItem("PCA", tabName = "PCAtab", icon = icon("fa fa-spinner fa-spin fa-2x fa-fw"))
      )
    ),
    dashboardBody(
      tabItems(
        # First tab content
        if (TRUE){
          tabItem(tabName = "upload",fluidPage(
            titlePanel("Uploading Files"),
            sidebarLayout(
              sidebarPanel(
                fileInput('file1', 'Choose File',
                          accept=c('text/csv', 
                                   'text/comma-separated-values,text/plain', 
                                   '.csv')),
                tags$hr(),
                
                radioButtons('sep', 'Separator',
                             c(Comma=',',
                               Semicolon=';',
                               Tab='\t'),
                             ','),
                radioButtons('quote', 'Quote',
                             c(None='',
                               'Double Quote'='"',
                               'Single Quote'="'"),
                             '"')
              ),
              mainPanel(
                DT::dataTableOutput("contents")
              )
            )
          ))
        },
        # Second tab content
        if (TRUE){
          tabItem(tabName = "Families", fluidPage(
            sidebarLayout(
              
              
              sidebarPanel (checkboxGroupInput("checklist_family", label = h3("Cell-types"), 
                                               choices = list("CHC1" = "CHC1", "CHC2" = "CHC2", "Pv" = "Pv", "Sst-CR"="Sst-CR",
                                                              "Sst-Nos1"="Sst-Nos1","VIP-CCK"="VIP-CCK","VIP-CR"="VIP-CR"),
                                               selected = c("CHC1", "CHC2" , "Pv", "Sst-CR",
                                                            "Sst-Nos1","VIP-CCK","VIP-CR")),
                            selectInput("selectfam", label = h3("Select gene family"), 
                                        choices = families_list, 
                                        selected = "Sodium-Channels"),
                            
                            selectInput("selectfamscale", label = h3("Select scale of the plot"), 
                                        choices = list("Expression" = "Expression", "Log2" = "Log2"), 
                                        selected = "Expression"),
                            width = 2
              ),
              
              mainPanel(
                tabsetPanel(
                  tabPanel("Boxplot", fluidPage( div(style='height:700px; width:1000px; overflow-x: scroll',
                                                     plotOutput("fam_boxplot", width = "100%")))), 
                  
                  tabPanel("Barplot (mean values)", fluidPage( div(style='height:700px; width:1000px; overflow-x: scroll',
                                                                   plotlyOutput("fam_barplot2", width = "100%"))))
                  , 
                  tabPanel("Pie-chart",fluidPage(conditionalPanel(condition = "input.checklist_family.indexOf('CHC1') > -1" ,plotlyOutput("fam_pie_CHC1")),
                                                 conditionalPanel(condition = "input.checklist_family.indexOf('CHC2') > -1" ,plotlyOutput("fam_pie_CHC2")),
                                                 conditionalPanel(condition = "input.checklist_family.indexOf('Pv') > -1" ,plotlyOutput("fam_pie_Pv")),
                                                 conditionalPanel(condition = "input.checklist_family.indexOf('Sst-CR') > -1" ,plotlyOutput("fam_pie_Sst_CR")),
                                                 conditionalPanel(condition = "input.checklist_family.indexOf('Sst-Nos1') > -1" ,plotlyOutput("fam_pie_Sst_Nos1")),
                                                 conditionalPanel(condition = "input.checklist_family.indexOf('VIP-CCK') > -1" ,plotlyOutput("fam_pie_VIP_CCK")),
                                                 conditionalPanel(condition = "input.checklist_family.indexOf('VIP-CR') > -1" ,plotlyOutput("fam_pie_VIP_CR"))
                                                 
                  ))
                  
                  
                ))
              
            )
            
            
          ))
        },
        # Third tab content
        if (TRUE){
          tabItem(tabName = "Indiv", fluidPage(sidebarLayout(sidebarPanel (
            checkboxGroupInput("checklist_ind", label = h3("Cell-types"), 
                               choices = list("CHC1" = "CHC1", "CHC2" = "CHC2", "Pv" = "Pv", "Sst-CR"="Sst-CR",
                                              "Sst-Nos1"="Sst-Nos1","VIP-CCK"="VIP-CCK","VIP-CR"="VIP-CR"),
                               selected = c("CHC1", "CHC2" , "Pv", "Sst-CR",
                                            "Sst-Nos1","VIP-CCK","VIP-CR")),
            selectizeInput('IndSelectGenes', label = h3('Genes to plot:'), 
                           choices = NULL, multiple = TRUE, options = list(maxItems = 5)),
            
            selectInput("selectIndscale", label = h3("Select scale of the plot"), 
                        choices = list("Expression" = "Expression", "Log2" = "Log2"), 
                        selected = "Expression"),
            
            width = 2), 
            mainPanel(
              tabsetPanel(tabPanel("Expression per cell",fluidPage(conditionalPanel(condition = "output.IndSelectnumber > 0" ,plotlyOutput("indiv_gene1")),
                                                                   conditionalPanel(condition = "output.IndSelectnumber > 1" ,plotlyOutput("indiv_gene2")),
                                                                   conditionalPanel(condition = "output.IndSelectnumber > 2" ,plotlyOutput("indiv_gene3")),
                                                                   conditionalPanel(condition = "output.IndSelectnumber > 3" ,plotlyOutput("indiv_gene4")),
                                                                   conditionalPanel(condition = "output.IndSelectnumber > 4" ,plotlyOutput("indiv_gene5"))
              )),
              tabPanel("Boxplot", fluidPage( div(style='height:700px; width:1000px; overflow-x: scroll',
                                                 plotOutput("IndBox", width = "100%")))),
              
              tabPanel("Barplot (mean values)", fluidPage( div(style='height:700px; width:1000px; overflow-x: scroll',
                                                               plotlyOutput("IndBar", width = "100%")))),
              
              tabPanel("Pie-chart",fluidPage(conditionalPanel(condition = "output.checklist_ind.indexOf('CHC1') > -1" ,plotlyOutput("ind_pie_CHC1")),
                                             conditionalPanel(condition = "output.checklist_ind.indexOf('CHC2') > -1" ,plotlyOutput("ind_pie_CHC2")),
                                             conditionalPanel(condition = "output.checklist_ind.indexOf('Pv') > -1" ,plotlyOutput("ind_pie_Pv")),
                                             conditionalPanel(condition = "output.checklist_ind.indexOf('Sst-CR') > -1" ,plotlyOutput("ind_pie_Sst_CR")),
                                             conditionalPanel(condition = "output.checklist_ind.indexOf('Sst-Nos1') > -1" ,plotlyOutput("ind_pie_Sst_Nos1")),
                                             conditionalPanel(condition = "output.checklist_ind.indexOf('VIP-CCK') > -1" ,plotlyOutput("ind_pie_VIP_CCK")),
                                             conditionalPanel(condition = "output.checklist_ind.indexOf('VIP-CR') > -1" ,plotlyOutput("ind_pie_VIP_CR"))
                                             
              ))
              
              
              
              )))))
          #if(TRUE){fluidPage(conditionalPanel(condition = "output.IndSelectnumber > -1" ,plotlyOutput("fam_pie_CHC1"))#,
          # conditionalPanel(condition = "output.IndSelectnumber > 1" ,print("cucu1")),
          #conditionalPanel(condition = "output.IndSelectnumber > 2" ,print("cucu1")),
          #conditionalPanel(condition = "output.IndSelectnumber > 3" ,print("cucu1")),
          #conditionalPanel(condition = "output.IndSelectnumber > 4" ,print("cucu1"))
          #)))}
          
          #))))
          #verbatimTextOutput("cucu")
        },
        #Fourth
        if (TRUE){
          tabItem(tabName = "PCAtab", fluidPage(sidebarLayout(sidebarPanel (
            checkboxGroupInput("checklist_PCA", label = h3("Cell-types"), 
                               choices = list("CHC1" = "CHC1", "CHC2" = "CHC2", "Pv" = "Pv", "Sst-CR"="Sst-CR",
                                              "Sst-Nos1"="Sst-Nos1","VIP-CCK"="VIP-CCK","VIP-CR"="VIP-CR"),
                               selected = c("CHC1", "CHC2" , "Pv", "Sst-CR",
                                            "Sst-Nos1","VIP-CCK","VIP-CR")),
            selectInput("selectPCAscale", label = h3("Select scale of the expression values"), 
                        choices = list("Expression" = "Expression", "Log2" = "Log2"), 
                        selected = "Expression"),
            selectInput("selectPCAmethod", label = h3("Select transformation"), 
                        choices = list("2DPCA" = "2DPCA", "3D" = "3D"), 
                        selected = "PCA"),
            # selectInput("selectPCAcolor", label = h3("Select coloring factor"), 
            #             choices = list("Cell type" = "Cell type", "Batch" = "Batch"), 
            #             selected = "Batch"),
            
            width = 2), 
            mainPanel(
              tabsetPanel(tabPanel("All genes",fluidPage(plotlyOutput("PCA1", height = "1000px"), plotlyOutput("scree1"),
                                                         plotlyOutput("cont1"))),
                          tabPanel("Per family",fluidPage(selectInput("PCAselectfam", label = h3("Select gene family"), 
                                                                      choices = families_list, 
                                                                      selected = "Sodium-Channels"),
                                                          plotlyOutput("PCA2", height = "1000px"),
                                                          plotlyOutput("scree2"),
                                                          plotlyOutput("cont2"))),
                          tabPanel("Selected genes",fluidPage(selectizeInput('PCASelectGenes', label = h3('Genes to plot:'), 
                                                                             choices = NULL, multiple = TRUE, options = list(maxItems = 50)),
                                                              plotlyOutput("PCA3", height = "1000px"),
                                                              plotlyOutput("scree3"),
                                                              plotlyOutput("cont3")))
                          
                          
                          
              )))))},
        #Fifth
        if(TRUE){
          tabItem(tabName = "DEA", fluidPage(titlePanel("Select your groups and genes and click the button for Differential Expression Analysis!"),sidebarLayout(sidebarPanel (
            
            actionButton("Help_button", "Help!"),
            bsPopover(id="Help_button", title="Differential Expression Analysis", content="Genes differentially expressed will be displayed once you click the Calculate button. Each batch assigned to Group 1 will be tested against the whole set for Group 2 and viceversa using a Wilcoxon test.The null hypothesis is that both samples (batch and set) belong to populations with the same distribution. For each gene, there is a n number of batches that are allowed to fail the test to still consider the gene differentially expressed.",
                      
                      placement = "bottom", trigger = "click",
                      options = NULL),
            checkboxGroupInput("check_Group1", label = h3("Group 1"), 
                               choices = list("CHC1" = "CHC1", "CHC2" = "CHC2", "Pv" = "Pv", "Sst-CR"="Sst-CR",
                                              "Sst-Nos1"="Sst-Nos1","VIP-CCK"="VIP-CCK","VIP-CR"="VIP-CR"),
                               selected = c("CHC1")),
            checkboxGroupInput("check_Group2", label = h3("Group 2"), 
                               choices = list("CHC1" = "CHC1", "CHC2" = "CHC2", "Pv" = "Pv", "Sst-CR"="Sst-CR",
                                              "Sst-Nos1"="Sst-Nos1","VIP-CCK"="VIP-CCK","VIP-CR"="VIP-CR"),
                               selected = c("Pv")),
            selectInput("DEselectfam", label = h3("Select gene family"), 
                        choices =DE_list , 
                        selected = "All genes"),
            conditionalPanel(condition="input.DEselectfam == 'Specific genes'", 
                             #helpText("TAB 2 SELECTED")
                             selectizeInput('DEspecific', label = h3('Select genes:'), 
                                            choices = NULL, multiple = TRUE, options = list(maxItems = 200))
            ),
            bsTooltip(id = "DEselectfam", title = "Select a gene family or specific genes (an input box will appear below)", placement = "left", trigger = "hover"),
            sliderInput("pvalue", "p-value:", 0, 1, 0.01,
                        step = 0.01, animate=
                          animationOptions(interval=300, loop=TRUE)),
            bsTooltip(id = "pvalue", title = "p-value to be used for each Wilcoxon test that will be performed", placement = "left", trigger = "hover"),
            sliderInput("n_value", "Allowed number of batches in disagreement", min=0, max=3, value=0),
            bsTooltip(id = "n_value", title = "This is the number n of batches allowed to fail the test for Differential Expresion per group. A value of 0 means all batches must be differentially expressed for a given gene. Recommended values: 0 or 1", placement = "left", trigger = "hover"),
            actionButton("Calc_DE", "Calculate!")
          ), mainPanel(
            
            
            conditionalPanel(condition="input.Calc_DE > 0", 
                             helpText("Differentially expressed genes:"),
                             verbatimTextOutput("nText"),plotOutput("venn", width="500px"),
                             bsTooltip(id = "venn", title = "Intersection shows the number of genes that are not differentially expressed", placement = "below", trigger = "hover"),
                             #pyramid
                             uiOutput("renderpyr"),
                             #plotlyOutput("DEpyramid2"),
                             DT::dataTableOutput("DE_table")
            )
          ))))
        }
        
      )))
)