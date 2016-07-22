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

ui <- dashboardPage(
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

###################################################################################################################################
########################################SERVER


server <- function(input, output,session) { 
  library(FactoMineR)
  raw_table<-data.frame()
  #First tab
  if (TRUE){
  output$contents <- DT::renderDataTable({
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    raw_table<<-read.csv(inFile$datapath, sep=input$sep, 
                        quote=input$quote)
    #gene_list<-colnames(raw_table[5:length(raw_table)])
    
    
    updateSelectizeInput(session, inputId="IndSelectGenes",
                         choices= colnames(raw_table[5:length(raw_table)]), selected = colnames(raw_table[5:length(raw_table)])[1:5] )
    updateSelectizeInput(session, inputId="PCASelectGenes",
                         choices= colnames(raw_table[5:length(raw_table)]), selected = colnames(raw_table[5:length(raw_table)])[1:5] )
    updateSelectizeInput(session, inputId="DEspecific",
                         choices= colnames(raw_table[5:length(raw_table)]), selected = colnames(raw_table[5:length(raw_table)])[1:5] )
    
    #plotPCA(raw_table)
    
    #output$data<-raw_table
    DT::datatable(raw_table, options = list(scrollX = TRUE, searching= FALSE))
    
    
  })
  
  }
  
  #Second tab
  if(TRUE){
    
    if(TRUE) ##Plotting functions
    { #famBoxplot<-function(input_df){
      getPalette = colorRampPalette(brewer.pal(9, "Set1"))
      famBoxplot<-reactive({
      library(reshape)
      library(ggplot2)
        #come back
        
        if(input$selectfamscale=="Expression"){
      input_df<-getDatafam()
      gg_df<-melt(cbind(celltype = input_df$Cell.type, input_df[5:length(input_df)]))
      gg1<-ggplot(gg_df, aes(variable, value, fill=celltype, dodge=celltype)) +
        geom_boxplot(outlier.size=0, position=position_dodge())+ylab("Expression (TPM)")+
        xlab("Genes")+theme(panel.margin = unit(0.5, "in"),
                            legend.position="left",axis.text=element_text(size=20),
                            axis.title=element_text(size=20))+
        ggtitle(paste0(input$selectfam, " family (Boxplots of gene expression)")) + 
        theme(plot.title = element_text(size=20,lineheight=0.8, face="bold"))
        
      #+ggtitle("Boxplots of expression for"+as.character(input$selectfam) ) 
      #colnames(gg_df)
        }else{
          input_df<-getDatafam()
          gg_df<-melt(cbind(celltype = input_df$Cell.type, input_df[5:length(input_df)]))
          gg_df$logval<-log2(gg_df$value)
          gg1<<-ggplot(gg_df, environment = environment(),aes(variable, logval, fill=celltype, dodge=celltype)) +
            geom_boxplot(outlier.size=0, position=position_dodge())+ylab("Expression log2(TPM+1)")+xlab("Genes")+
            xlab("Genes")+theme(panel.margin = unit(0.5, "in"),
                                legend.position="left",axis.text=element_text(size=20),
                                axis.title=element_text(size=20))+
            ggtitle(paste0(input$selectfam, " family (Boxplots of gene expression)")) + 
            theme(plot.title = element_text(size=20,lineheight=0.8, face="bold"))
          #gg1<-gg1+theme(legend.position="left")
             #theme(legend.position="top")
          # +
            #scale_y_continuous(trans="log2")
        }
        
        gg1
        
        
      })
      famBarplot<-reactive({
        library(reshape)
        library(ggplot2)
       
        library(plyr)
        input_df<-getDatafam()
        if(input$selectfamscale=="Expression"){
          
          working_df<-melt(cbind(celltype = input_df$Cell.type, input_df[5:length(input_df)]))
          gg_df<-ddply(working_df, c("celltype","variable"), summarise,
                       mean_val=mean(value))
          gg1<-ggplot(gg_df, aes(variable, mean_val, fill=celltype, dodge=celltype)) +
            geom_bar(stat="identity",position=position_dodge())+ylab("Expression (TPM)")+
            xlab("Genes")+theme(axis.text=element_text(size=20),
                                axis.title=element_text(size=20))+theme(panel.margin = unit(0.5, "in"),
                                                                        legend.position="left",axis.text=element_text(size=20),
                                                                        axis.title=element_text(size=20),
                                                                        plot.title = element_text(size=20,lineheight=0.8, face="bold"))+
            ggtitle(paste0(input$selectfam, " family (Barplots of gene expression)")) 
          
          
        }
        else{
          working_df<-melt(cbind(celltype = input_df$Cell.type, input_df[5:length(input_df)]))
          working_df$logval<-log2(working_df$value)
          gg_df<-ddply(working_df, c("celltype","variable"), summarise,
                       mean_val=mean(value),
                       mean_log=mean(logval))
          gg1<-ggplot(gg_df, aes(variable, mean_log, fill=celltype, dodge=celltype)) +
            geom_bar(stat="identity",position=position_dodge())+ylab("Expression log2(TPM+1)")+xlab("Genes")+
            xlab("Genes")+theme(axis.text=element_text(size=20),
                                axis.title=element_text(size=20))+theme(panel.margin = unit(0.5, "in"),
                                                                        legend.position="left",axis.text=element_text(size=20),
                                                                        axis.title=element_text(size=20))+
            ggtitle(paste0(input$selectfam, " family (Barplots of gene expression)")) + 
            theme(plot.title = element_text(size=20,lineheight=0.8, face="bold"))
          
        }
      })
      pieplot<-function(cellt){
        m = list(
          l = 50,
          r = 50,
          b = 50,
          t = 50,
          pad = 4
        )
        dataf<<-famPieData()
        dataf<-subset(dataf, dataf$celltype==cellt)
        if(input$selectfamscale=="Expression"){
          cols <- getPalette(length(dataf$variable))
          plot_ly(data=dataf, labels = variable, values=mean_val , type = "pie",sort=FALSE,
                   name = "CHC1", showlegend = T )%>%
            layout(title = paste0(cellt," ",input$selectfam, " relative expression (TPM)"), margin = m)
        }else{
          plot_ly(data=dataf, labels = variable, values=mean_log , type = "pie", sort=FALSE, name = "CHC1", showlegend = T)%>%
            layout(title = paste0(cellt," ",input$selectfam," relative expression (log2(TPM+1))"), margin = m)
        }
      }
      
    }
    
    #Checklist_family
    if (TRUE) #Family lists
    {genes_per_family<-list()
    
    genes_per_family[["GABAa Receptors (ionotropic)"]]<-c(1:18)
    genes_per_family[["GABAb Receptors (metabotropic)"]]<-c(19)
    genes_per_family[["Not-NMDA Glutamate Rec. (ionotropic)"]]<-c(20:25,33:37)
    genes_per_family[["NMDA Glutamate Rec. (ionotropic)"]]<-c(26:32)
    genes_per_family[["Metabotropic Glutamate Receptors"]]<-c(38:44)
    #genes_per_family[["Receptors"]]<-c(1:44)
    genes_per_family[["Potassium-Channels"]]<-c(158:230)
    genes_per_family[["Sodium-Channels"]]<-c(263:276)
    genes_per_family[["Calcium-Channels"]]<-c(141:149)
    genes_per_family[["CNG-Channels"]]<-c(150:157)
    genes_per_family[["Synaptic release machinery genes"]]<-c(281:518)
    genes_per_family[["Ion channels and receptors"]]<-c(1:280)
    }
    
      reac_fam_cell <-  reactive({
      
      (dim(subset(raw_table, (raw_table$Cell.type %in% input$checklist_family))))
      
    })
      getDatafam<-reactive({
        current_indexes<-genes_per_family[[input$selectfam]]
        ##AQUI
        gene_list<-colnames(raw_table[5:length(raw_table)])
        current_genes<-gene_list[current_indexes]
        #Subset cells
        output_fam_df<-subset(raw_table, (raw_table$Cell.type %in% input$checklist_family))
        #Subset cols
        output_fam_df<-cbind(output_fam_df[1:4],output_fam_df[current_genes])
        #dim(output_fam_df)
        
        
      })
      
      famPieData<-reactive({
        library(reshape)
        library(ggplot2)
        
        library(plyr)
        input_df<-getDatafam()
        if(input$selectfamscale=="Expression"){
          
          working_df<-melt(cbind(celltype = input_df$Cell.type, input_df[5:length(input_df)]))
          return_df<-ddply(working_df, c("celltype","variable"), summarise,
                       mean_val=mean(value))
          }else{
          working_df<-melt(cbind(celltype = input_df$Cell.type, input_df[5:length(input_df)]))
          working_df$logval<-log2(working_df$value)
          return_df<-ddply(working_df, c("celltype","variable"), summarise,
                       mean_val=mean(value),
                       mean_log=mean(logval))
          }
          
        
      })
      
      #Boxplot renderPlot
    output$fam_boxplot<-renderPlot(print(famBoxplot()),height = 685,width=function(){
      w<-length(colnames(getDatafam()))*100
      })
    ##Barplot render
    output$fam_barplot<-renderPlot(print(famBarplot()),width=function(){
      w<-length(colnames(getDatafam()))*100
    })
    #output$fam_barplot2<-renderPlotly(ggplotly(famBarplot()))
    output$fam_barplot_width<-reactive(paste0(as.character(dim(getDatafam())[1]*10),"px"))
    output$fam_barplot2<-renderPlotly(
      {
        ggp_build <- plotly_build(famBarplot())
        
        ggp_build$layout$height = 685
        ggp_build$layout$width = length(colnames(getDatafam()))*100
        #ggp_build$layout$width = 20000
        #ggp_build$layout$legend$xanchor = "left"
        ggp_build$layout$legend$x <- -0.125
        #ggp_build$layout$annotations[["text"]]<- "kele"
        ggp_build # %>% layout(legend = list(x = 0.5, y = 0))
        #ggp <-(ggplotly(famBoxplot(), width = 20000)%>%layout(autosize=TRUE))
        #ggp$x$layout$width = 20000
        #ggp
        
      })
                                  #width=paste0(as.character(dim(getDatafam())[2])*100),"px")
    #output$trialtext<-renderText({famBoxplot()})
    if(TRUE) ##Pie-charts
      {output$fam_pie_CHC1<-renderPlotly(pieplot("CHC1"))
      output$fam_pie_CHC2<-renderPlotly(pieplot("CHC2"))
      output$fam_pie_Pv<-renderPlotly(pieplot("Pv"))
      output$fam_pie_Sst_CR<-renderPlotly(pieplot("Sst-CR"))
      output$fam_pie_Sst_Nos1<-renderPlotly(pieplot("Sst-Nos1"))
      output$fam_pie_VIP_CCK<-renderPlotly(pieplot("VIP-CCK"))
      output$fam_pie_VIP_CR<-renderPlotly(pieplot("VIP-CR"))
      
        
    }
    
  
  }
  
  #Third tab
  if(TRUE){
    Ind_sel_n<-reactive({
      length(input$IndSelectGenes)
    })
    getInData<-reactive({
      
      current_genes<-input$IndSelectGenes
      #Subset cells
      output_fam_df<-subset(raw_table, (raw_table$Cell.type %in% input$checklist_ind))
      #Subset cols
      output_ind_df<-cbind(output_fam_df[1:4],output_fam_df[current_genes])
      #dim(output_fam_df)
    }) 
    if (TRUE) #Plotting functions
    {
      percell_plot<-function(gene){
      library(reshape)
      library(ggplot2)
       
      
      gg_df<-getInData()
      gg_df$cellnumber<-gg_df$Sl..No
      if(input$selectIndscale=="Expression"){
        gg_df$expression<-gg_df[gene]
        gg_df$cellnumber<-gg_df$Sl..No
        g1<-ggplot(data=gg_df, aes(x=cellnumber, y=expression, fill=Cell.type) ) +
        geom_bar( stat="identity") + labs(fill = "Cell_type")+ xlab("") + ylab(paste0(gene," expression (TPM)"))+
        theme(panel.margin = unit(1, "in"),legend.position="left",axis.text=element_text(size=10),
        axis.title=element_text(size=10))+ ggtitle(paste0(gene, " single-cell expression")) + 
          theme(plot.title = element_text(size=20,lineheight=0.8, face="bold"))
        
      }else{
        gg_df_log<-cbind(gg_df[1:4],log2(gg_df[5:length(gg_df)]) )
        gg_df_log$expression<-gg_df_log[gene]
        gg_df_log$cellnumber<-gg_df$Sl..No
        g1<-ggplot(data=gg_df_log, aes(x=cellnumber, y=expression, fill=Cell.type) ) +
          geom_bar( stat="identity") + labs(fill = "Cell_type")+ xlab("") + ylab(paste0(gene," expression (log2(TPM+1))"))+
          theme(panel.margin = unit(1, "in"),legend.position="left",axis.text=element_text(size=10),
                axis.title=element_text(size=10))+ ggtitle(paste0(gene, " single-cell expression")) + 
          theme(plot.title = element_text(size=20,lineheight=0.8, face="bold"))
        
       }
      
      
      
      
      }
      indBoxplot<-reactive({
        library(reshape)
        library(ggplot2)
        #dev.off()
        input_df<-getInData()
        if(input$selectIndscale=="Expression"){
          
          gg_df<-melt(cbind(celltype = input_df$Cell.type, input_df[5:length(input_df)]))
          gg1<-ggplot(gg_df, aes(variable, value, fill=celltype, dodge=celltype)) +
            geom_boxplot(outlier.size=0, position=position_dodge())+ylab("Expression (TPM)")+
            xlab("Genes")+theme(panel.margin = unit(0.5, "in"),
                                legend.position="left",axis.text=element_text(size=20),
                                axis.title=element_text(size=20))+
            ggtitle(paste0(paste(input$IndSelectGenes, collapse=", "), " (Boxplots of gene expression)")) + 
            theme(plot.title = element_text(size=20,lineheight=0.8, face="bold"))
          
          }else
            {
          
          gg_df<-melt(cbind(celltype = input_df$Cell.type, input_df[5:length(input_df)]))
          gg_df$logval<-log2(gg_df$value)
          gg1<<-ggplot(gg_df, environment = environment(),aes(variable, logval, fill=celltype, dodge=celltype)) +
            geom_boxplot(outlier.size=0, position=position_dodge())+ylab("Expression log2(TPM+1)")+xlab("Genes")+
            xlab("Genes")+theme(panel.margin = unit(0.5, "in"),
                                legend.position="left",axis.text=element_text(size=20),
                                axis.title=element_text(size=20))+
            ggtitle(paste0(paste(input$IndSelectGenes, collapse=", "), " (Boxplots of gene expression)")) + 
            theme(plot.title = element_text(size=20,lineheight=0.8, face="bold"))
          #gg1<-gg1+theme(legend.position="left")
          #theme(legend.position="top")
          # +
          #scale_y_continuous(trans="log2")
        }
      })
      indBarplot<-reactive({
        library(reshape)
        library(ggplot2)
        #dev.off()
        library(plyr)
        input_df<-getInData()
        if(input$selectIndscale=="Expression"){
          
          working_df<-melt(cbind(celltype = input_df$Cell.type, input_df[5:length(input_df)]))
          gg_df<-ddply(working_df, c("celltype","variable"), summarise,
                       mean_val=mean(value))
          gg1<-ggplot(gg_df, aes(variable, mean_val, fill=celltype, dodge=celltype)) +
            geom_bar(stat="identity",position=position_dodge())+ylab("Expression (TPM)")+
            xlab("Genes")+theme(axis.text=element_text(size=20),
                                axis.title=element_text(size=20))+theme(panel.margin = unit(0.5, "in"),
                                                                        legend.position="left",axis.text=element_text(size=20),
                                                                        axis.title=element_text(size=20),
                                                                        plot.title = element_text(size=20,lineheight=0.8, face="bold"))+
            ggtitle(paste0(paste(input$IndSelectGenes, collapse=", "), " (Barplots of gene expression)")) 
          
          
        }
        else{
          working_df<-melt(cbind(celltype = input_df$Cell.type, input_df[5:length(input_df)]))
          working_df$logval<-log2(working_df$value)
          gg_df<-ddply(working_df, c("celltype","variable"), summarise,
                       mean_val=mean(value),
                       mean_log=mean(logval))
          gg1<-ggplot(gg_df, aes(variable, mean_log, fill=celltype, dodge=celltype)) +
            geom_bar(stat="identity",position=position_dodge())+ylab("Expression log2(TPM+1)")+xlab("Genes")+
            xlab("Genes")+theme(axis.text=element_text(size=20),
                                axis.title=element_text(size=20))+theme(panel.margin = unit(0.5, "in"),
                                                                        legend.position="left",axis.text=element_text(size=20),
                                                                        axis.title=element_text(size=20))+
            ggtitle(paste0(paste(input$IndSelectGenes, collapse=", "), " (Barplots of gene expression)")) + 
            theme(plot.title = element_text(size=20,lineheight=0.8, face="bold"))
          
        }
      })
      
      indPieData<-reactive({
        library(reshape)
        library(ggplot2)
        #dev.off()
        library(plyr)
        input_df<-getInData()
        if(input$selectIndscale=="Expression"){
          
          working_df<-melt(cbind(celltype = input_df$Cell.type, input_df[5:length(input_df)]))
          return_df<-ddply(working_df, c("celltype","variable"), summarise,
                           mean_val=mean(value))
        }else{
          working_df<-melt(cbind(celltype = input_df$Cell.type, input_df[5:length(input_df)]))
          working_df$logval<-log2(working_df$value)
          return_df<-ddply(working_df, c("celltype","variable"), summarise,
                           mean_val=mean(value),
                           mean_log=mean(logval))
        }
        
        
      })
      ind_pieplot<-function(cellt){
        m = list(
          l = 50,
          r = 50,
          b = 50,
          t = 50,
          pad = 4
        )
        dataf<-indPieData()
        dataf<-subset(dataf, dataf$celltype==cellt)
        if(input$selectIndscale=="Expression"){
          cols <- getPalette(length(dataf$variable))
          plot_ly(data=dataf, labels = variable, values=mean_val , type = "pie",sort=FALSE,
                  name = "CHC1", showlegend = T )%>%
            layout(title = paste0(cellt," relative expression of ", paste(input$IndSelectGenes, collapse=", "),  "(TPM)"), margin = m)
        }else{
          plot_ly(data=dataf, labels = variable, values=mean_log , type = "pie", sort=FALSE, name = "CHC1", showlegend = T)%>%
            layout(title = paste0(cellt," relative expression of ", paste(input$IndSelectGenes, collapse=", "),  "(log2(TPM+1))"), margin = m)
        }
      }
    }
    output$IndSelectnumber<-reactive(length(input$IndSelectGenes))
    #output$test<-renderText({ "que si"})
    outputOptions(output, "IndSelectnumber", suspendWhenHidden=FALSE)
    
    if(TRUE) #Per cell plots
    {
      output$indiv_gene1<-renderPlotly({
        
        ggplotly(percell_plot(input$IndSelectGenes[1]))
        
        # gg_df<-raw_table
        # gg_df$expression<-raw_table[input$IndSelectGenes[1]]
        # gg_df$cellnumber<-raw_table$Sl..No
        # 
        # g1<-ggplot(data=gg_df, aes(x=cellnumber, y=expression, fill=Cell.type) ) +
        #   geom_bar( stat="identity") + labs(fill = "Cell_type")+ xlab("") + ylab(paste0(input$IndSelectGenes[1]," expression (TPM)"))
        # plotly_build(g1)
      })
      output$indiv_gene2<-renderPlotly({
        
        ggplotly(percell_plot(input$IndSelectGenes[2]))
        
      })
      output$indiv_gene3<-renderPlotly({
        
        ggplotly(percell_plot(input$IndSelectGenes[3]))
        
      })
      output$indiv_gene4<-renderPlotly({
        
        ggplotly(percell_plot(input$IndSelectGenes[4]))
        
      })
      output$indiv_gene5<-renderPlotly({
        
        ggplotly(percell_plot(input$IndSelectGenes[5]))
        
      })
    }
    output$IndBox<-renderPlot(print(indBoxplot()),height = 685,width=function(){
      w<-length(colnames(getDatafam()))*100
    })
    output$IndBar<-renderPlotly({
        ggp_build <- plotly_build(indBarplot())
        ggp_build$layout$height = 685
        ggp_build$layout$width = length(colnames(getDatafam()))*100
        #ggp_build$layout$width = 20000
        #ggp_build$layout$legend$xanchor = "left"
        ggp_build$layout$legend$x <- -0.125
        #ggp_build$layout$annotations[["text"]]<- "kele"
        ggp_build # %>% layout(legend = list(x = 0.5, y = 0))
        #ggp <-(ggplotly(famBoxplot(), width = 20000)%>%layout(autosize=TRUE))
        #ggp$x$layout$width = 20000
        #ggp
        
      })
    output$checklist_ind<-reactive(input$checklist_ind)
    outputOptions(output, "checklist_ind", suspendWhenHidden=FALSE)
    
    if(TRUE) ##Pie-charts
    {output$ind_pie_CHC1<-renderPlotly(ind_pieplot("CHC1"))
    output$ind_pie_CHC2<-renderPlotly(ind_pieplot("CHC2"))
    output$ind_pie_Pv<-renderPlotly(ind_pieplot("Pv"))
    output$ind_pie_Sst_CR<-renderPlotly(ind_pieplot("Sst-CR"))
    output$ind_pie_Sst_Nos1<-renderPlotly(ind_pieplot("Sst-Nos1"))
    output$ind_pie_VIP_CCK<-renderPlotly(ind_pieplot("VIP-CCK"))
    output$ind_pie_VIP_CR<-renderPlotly(ind_pieplot("VIP-CR"))
    
    
    }

  }
  
  #Fourth tab
  if(TRUE){
    getPCAData1<-reactive({
      #current_genes<-input$IndSelectGenes
      #Subset cells
      output_fam_df<-subset(raw_table, (raw_table$Cell.type %in% input$checklist_PCA))
      #Subset cols
      #output_ind_df<-cbind(output_fam_df[1:4],output_fam_df[current_genes])
      #dim(output_fam_df)
    })
    
    getPCAData2<-reactive({
      current_indexes<-genes_per_family[[input$PCAselectfam]]
      gene_list<-colnames(raw_table[5:length(raw_table)])
      current_genes<-gene_list[current_indexes]
      #Subset cells
      output_fam_df<-subset(raw_table, (raw_table$Cell.type %in% input$checklist_PCA))
      #Subset cols
      output_ind_df<-cbind(output_fam_df[1:4],output_fam_df[current_genes])
      #dim(output_fam_df)
    })
    
    getPCAData3<-reactive({
      current_genes<-input$PCASelectGenes
      #Subset cells
      output_fam_df<-subset(raw_table, (raw_table$Cell.type %in% input$checklist_PCA))
      #Subset cols
      output_ind_df<-cbind(output_fam_df[1:4],output_fam_df[current_genes])
      #dim(output_fam_df)
    })
    
    #Plotting functions
    if (TRUE){
    
      calc_PCA<-function(data){
        
        if(input$selectPCAscale=="Expression"){
          input_data<-data
        }else{
          input_data<-cbind(data[1:4],log2(data[5:length(data)]) )
          }
        #Perform PCA
        meta_data<-input_data[1:4]
        
        gene_data<-input_data[5:length(input_data)]
        gene_data<-Filter(function(x) sd(x) != 0, gene_data)
        res.pca<-PCA(gene_data,graph = FALSE)
        reac_pca<<-reactive(res.pca)
        PC1 <- res.pca$ind$coord[,1]
        PC2 <- res.pca$ind$coord[,2]
        PC3 <- res.pca$ind$coord[,3]
        labs <- rownames(res.pca$ind$coord)
        #PCs <- data.frame(cbind(PC1,PC2,PC3))
        PCs <- data.frame(cbind(PC1,PC2,PC3, meta_data))
        #PCs$identifier<-paste(meta_data$Batch, meta_data$Cell.type, rownames(meta_data), sep = "_")
        #batch_celltype_celln<-paste(meta_data$Batch, meta_data$Cell.type, rownames(meta_data), sep = "_")
        #rownames(PCs) <-
        return(PCs)
      } 
      
      plot2DPCA<-function(PCs){
        library(rwantshue)
        
        scheme <- iwanthue()
        col_hcl <- list(c(0, 360),   # hue range [0,360]
                        c(0, 3),     # chroma range [0,3]
                        c(0, 1.5))   # lightness range [0,1.5]
        cols <- scheme$hex(n = length(unique(PCs$Batch)), force_mode = FALSE, quality = 50, color_space = col_hcl)
        cell_number<-rownames(PCs)
        Celltype<-factor(PCs$Cell.type)
        Batch<-factor(PCs$Batch)
        ggpc1<-ggplot(PCs, aes(PC1,PC2, label=cell_number, shape= Celltype, color= factor(Batch)))  +
          #scale_shape_manual(values=1:nlevels(Celltype))+ geom_point(size=5)+scale_color_manual(values=cols)
          scale_color_manual(values=cols) +scale_shape_manual(values=1:nlevels(PCs$Cell.type))+ geom_point(size=5) 
        
      }
      
      plot3DPCA<-function(PCs){
        
        
        
        Celltype<-factor(PCs$Cell.type)
        #paste("Clarity: ", clarity)
        
        plot_ly(PCs, x = PC1, y = PC2, z= PC3, text = paste("Cell type: ",as.character(Cell.type)), 
                color = as.character(Batch), symbol = factor(Cell.type) ,type="scatter3d", 
                symbols= c(  "cross","x","circle-open" , "square-open" ,"circle", "diamond" , "diamond-open"))
                #color = Cell.type, type="scatter3d")
        
        #plot_ly(df, x = x, y = y, z = z, type = "scatter3d", mode = "markers")
        #ggpc1<-ggplot(PCs, aes(PC1,PC2, label=cell_number, shape= Celltype, color= factor(Batch)))  +
          #scale_shape_manual(values=1:nlevels(Celltype))+ geom_point(size=5)+scale_color_manual(values=cols)
         # scale_color_manual(values=cols) +scale_shape_manual(values=1:nlevels(PCs$Cell.type))+ geom_point(size=5) 
        
      }
      
      calc_PC_contribution<-function(data){
        calc_PCA(data)
        library(factoextra)
        fviz_screeplot(reac_pca(), ncp=9)
      }
      
      calc_PC_variables<-function(data){
        calc_PCA(data)
        library(factoextra)
        fviz_pca_contrib(reac_pca(), choice = "var", axes = 1:3, top = 20)+
          #theme(axis.text.x = element_text(size = 5))
          theme(axis.text.x = element_text(size = 5))+ggtitle("Contribution of gene expression levels to Dimensions 1,2 and 3")
        
      }
      
    }
    
    output$PCA1<-renderPlotly({
      if (input$selectPCAmethod=="2DPCA"){
        ggplotly(plot2DPCA(calc_PCA(getPCAData1())))
      }else{
        plot3DPCA(calc_PCA(getPCAData1()))
      }
      
      
    })
    output$scree1<-renderPlotly({
      input$selectPCAscale
      ggplotly(calc_PC_contribution(getPCAData1()))
    })
    output$cont1<-renderPlotly({
       ggplotly(calc_PC_variables(getPCAData1()))
       
     })
    
    output$PCA2<-renderPlotly({
      if (input$selectPCAmethod=="2DPCA"){
        ggplotly(plot2DPCA(calc_PCA(getPCAData2())))
      }else{
        plot3DPCA(calc_PCA(getPCAData2()))
      }
      
      
    })
    output$scree2<-renderPlotly({
      input$selectPCAscale
      ggplotly(calc_PC_contribution(getPCAData2()))
    })
    output$cont2<-renderPlotly({
      ggplotly(calc_PC_variables(getPCAData2()))
      
    })
    
    output$PCA3<-renderPlotly({
      if (input$selectPCAmethod=="2DPCA"){
        ggplotly(plot2DPCA(calc_PCA(getPCAData3())))
      }else{
        plot3DPCA(calc_PCA(getPCAData3()))
      }
      
      
    })
    output$scree3<-renderPlotly({
      input$selectPCAscale
      ggplotly(calc_PC_contribution(getPCAData3()))
    })
    output$cont3<-renderPlotly({
      ggplotly(calc_PC_variables(getPCAData3()))
      
    })
    
  }
  
  #Fifth tab
  if (TRUE){
    #Differential Expression function
    adaptedDE <- function(set1, set2, disag_n=0 ,threshold=0.01){
      # 2 dataframes with the cells to be analyzed. Column 1 must be "Batch"
      # and contain batch number or ID. Rest of the columns must be genenames and contain expression values
      #  expressed in normalized units (TPM). Each row is a cell. It returns list of DE genes
      gene_list<-colnames(set1[2:length(set1)])
      nDEbatch_1<-list()
      nDEbatch_2<-list()
      set1_blist<-unique(set1[,1])
      set2_blist<-unique(set2[,1])
      genes<-colnames(set1[2:length(set1)])
      ####set1
      
      for (gene in genes){
        p_values<-rep(1,length(unique(set1[,1])))
        index=1
        nDE<-0
        for (batch in set1_blist){
          current<-subset(set1, set1[,1] == batch)
          test_out<-wilcox.test(current[,gene], set2[,gene], exact = FALSE)
          current_p<-test_out$p.value
          if(is.nan(current_p)){
            current_p<-1
          }
          p_values[index]<-current_p
          index=index+1
        }
        adjusted_ps<-p.adjust(p_values, method = "BH")
        nDE<-sum(adjusted_ps< threshold)
        nDEbatch_1[[gene]]<-nDE
      }
      
      ##decide if consensus
      consensus_n<-length(unique(set1[,1]))-disag_n
      for (gene in gene_list){
        if (nDEbatch_1[[gene]]<consensus_n){
          nDEbatch_1[[gene]]<-NULL
        }
      }
      
      ###Set2
      
      for (gene in genes){
        p_values<-rep(1,length(unique(set2[,1])))
        index=1
        nDE<-0
        for (batch in set2_blist){
          current<-subset(set2, set2[,1] == batch)
          
          test_out<-wilcox.test(current[,gene], set1[,gene], exact = FALSE)
          current_p<-test_out$p.value
          if(is.nan(current_p)){
            current_p<-1
          }
          p_values[index]<-current_p
          index=index+1
        }
        adjusted_ps<-p.adjust(p_values, method = "BH")
        nDE<-sum(adjusted_ps< threshold)
        nDEbatch_2[[gene]]<-nDE
      }
      ##decide if consensus
      
      consensus_n<-length(unique(set2[,1]))-disag_n
      for (gene in gene_list){
        if (nDEbatch_2[[gene]]<consensus_n){
          nDEbatch_2[[gene]]<-NULL
        }
      }
      
      DEgenes_1<-names(nDEbatch_1)
      DEgenes_2<-names(nDEbatch_2)
      DE_list<-intersect(DEgenes_1,DEgenes_2)
      
      #return(nDEbatch_2)
      return(DE_list)
    }
    group1Data<-eventReactive(input$Calc_DE,{
      group1_df<-subset(raw_table, (raw_table$Cell.type %in% input$check_Group1))
      if(input$DEselectfam =="All genes"){
        group1_df<-group1_df
      }else if(input$DEselectfam =="Specific genes"){
        current_genes<-input$DEspecific
        group1_df<-cbind(group1_df[1:4],group1_df[current_genes])
        
      }else{
        current_indexes<-genes_per_family[[input$DEselectfam]]
        gene_list<-colnames(raw_table[5:length(raw_table)])
        current_genes<-gene_list[current_indexes]
        group1_df<-cbind(group1_df[1:4],group1_df[current_genes])
      }
    })
    group2Data<-eventReactive(input$Calc_DE,{
      group2_df<-subset(raw_table, (raw_table$Cell.type %in% input$check_Group2))
      if(input$DEselectfam =="All genes"){
        group2_df<-group2_df
      }else if(input$DEselectfam =="Specific genes"){
        current_genes<-input$DEspecific
        group2_df<-cbind(group2_df[1:4],group2_df[current_genes])
        
      }else{
        current_indexes<-genes_per_family[[input$DEselectfam]]
        gene_list<-colnames(raw_table[5:length(raw_table)])
        current_genes<-gene_list[current_indexes]
        group2_df<-cbind(group2_df[1:4],group2_df[current_genes])
      }
    })
    
    DEpyramid<-eventReactive(input$Calc_DE,{
      library(reshape)
      library(ggplot2)
      library(plyr)
      #dev.off()
      
      
      input_df<-getInData()
      
        
        working_df<-melt(cbind(celltype = input_df$Cell.type, input_df[5:length(input_df)]))
        gg_df<-ddply(working_df, c("celltype","variable"), summarise,
                     mean_val=mean(value))
        gg1<-ggplot(gg_df, aes(variable, mean_val, fill=celltype, dodge=celltype)) +
          geom_bar(stat="identity",position=position_dodge())+ylab("Expression (TPM)")+
          xlab("Genes")+theme(axis.text=element_text(size=20),
                              axis.title=element_text(size=20))+theme(panel.margin = unit(0.5, "in"),
                                                                      legend.position="left",axis.text=element_text(size=20),
                                                                      axis.title=element_text(size=20),
                                                                      plot.title = element_text(size=20,lineheight=0.8, face="bold"))+
          ggtitle(paste0(paste(input$IndSelectGenes, collapse=", "), " (Barplots of gene expression)")) 
        
        set.seed(1)
        df0 <- data.frame(Age = factor(rep(x = 1:10, times = 2)), 
                          Gender = rep(x = c("Female", "Male"), each = 10),
                          Population = sample(x = 1:100, size = 20))
        
        head(df0)
        #   Age Gender Population
        # 1   1 Female         27
        # 2   2 Female         37
        # 3   3 Female         57
        # 4   4 Female         89
        # 5   5 Female         20
        # 6   6 Female         86
        
        library("ggplot2")
        ggplot(data = df0, aes(x = Age, y = Population, fill = Gender)) +
          geom_bar(data = subset(df0, Gender=="Female"),
                   stat = "identity") +
          geom_bar(data = subset(df0, Gender=="Male"),
                   stat = "identity",
                   position = "identity",
                   mapping = aes(y = -Population)) +
          scale_y_continuous(labels = abs) +
          coord_flip()
      
      
    })
        
    ntext <- eventReactive(input$Calc_DE, {
      out_list<-adaptedDE(group1Data()[4:length(group1Data())],
                          group2Data()[4:length(group2Data())],
                          disag_n = input$n_value, threshold = input$pvalue)
      })
    
    output$nText <- renderText({
      if (length(ntext())>0){
        paste(ntext(), collapse = ", ")
      }else{
        out<- "No differentially expressed genes"
      }
      
      
    })
    
    output$venn<-renderPlot({#dev.off()
      venn.plot <- draw.pairwise.venn((length(group1Data())-4), 
    (length(group1Data())-4),(length(group1Data())-4)- length(ntext()), c(paste(unique(group1Data()$Cell.type), collapse = ", "), 
    paste(unique(group2Data()$Cell.type), collapse = ", ")), cat.just = list(c(-1, -1), c(1, 1)),
    fill = c("blue", "red"))})
    output$DE_table<- DT::renderDataTable({
      library(plyr)
      set1<-cbind( celltype=group1Data()$Cell.type,group1Data()[ntext()])
      set2<-cbind(celltype= group2Data()$Cell.type,group2Data()[ntext()])
      dset1<- ddply(set1, "celltype", numcolwise(mean))
      dset1<-cbind(Group=rep(1,nrow(dset1)),dset1)
      dset2<-ddply(set2, "celltype", numcolwise(mean))
      dset2<-cbind(Group=rep(2,nrow(dset2)), dset2)
      joint_set<-rbind(dset1,dset2)
      DT::datatable(joint_set, options = list(dom = 't',scrollX = TRUE, searching= FALSE))
                    
    })
    
    
    
   
     output$DEpyramid2<-renderPlotly({
      library(reshape2)
       
       
      Group1lab<-paste(unique(group1Data()$Cell.type), collapse = ", ")
      Group2lab<-paste(unique(group2Data()$Cell.type), collapse = ", ")
       
      set1<-cbind(Group=rep(Group1lab,nrow(group1Data())),group1Data()[ntext()])
      set2<-cbind(Group=rep(Group2lab,nrow(group2Data())),group2Data()[ntext()])
      dset1<-ddply(set1, "Group", numcolwise(mean))
      dset2<-ddply(set2, "Group", numcolwise(mean))
      dset1<-cbind(dset1[1],round(dset1[2:length(dset1)], digits = 2))
      dset2<-cbind(dset2[1],round(dset2[2:length(dset2)], digits = 2))
      
      
      mset1<-melt(dset1)
      
      s1max<-max(mset1$value)
      
      mset1$Gene<-mset1$variable
      mset1$value<-(-1)*mset1$value
      mset2<-melt(dset2)
      
      mset2$Gene<-mset2$variable
      s2max<-max(mset2$value)
      
      maxval<-round(max(s1max,s2max),-2)
      minval<-round(min(s1max,s2max),-1)
      
      
      tickvalues = sort(c(-maxval, (-0.8)*maxval,(-0.6*maxval),(-0.4)*maxval,-0.2*maxval, 0, 0.2*maxval, 
                          (0.4)*maxval,(0.6*maxval),(0.8)*maxval,maxval, minval))
      ticklabs = as.character(abs(tickvalues))
      
      sets<-rbind(mset1,mset2)
      
      sets$abs_value<-abs(sets$value)
      plot_ly(sets, x = value, y = Gene, group = Group, type = 'bar', orientation = 'h',
              hoverinfo = 'y+text+name', text = abs_value)%>%
        layout(bargap = 0.1, barmode = 'overlay',
               xaxis = list(tickmode = 'array', tickvals = tickvalues,
                            ticktext = ticklabs))
      #if (length(ntext())<1){
       # return(NULL)      }
    })
    output$renderpyr<-renderUI({
      if(length(ntext())>0){
      x= 100+ 25*length(ntext())
      list(
        plotlyOutput("DEpyramid2", height = x)
      )}
      
    })
   }
}

#shinyApp(ui, server)
shinyApp(ui, server, options=list(shiny.port = 2020, host="0.0.0.0"))
#options=list(port = 7777

####################################################
