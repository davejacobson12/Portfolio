### This app has epi curve and map on the same plot

library(shiny) #yes
library(shinydashboard)
library(countrycode)
library(ggplot2) #yes
library(dplyr) #yes
# library(usmap) #yes
# library(readxl) #yes
# library(gridExtra) #yes
library(tidyr) #yes
library(plotly) #yes
# library(usdata) # may be useful for abbreviation conversion
# library(maditr)
library(stringr)
library(RColorBrewer)
library(shinythemes)
library(DT)
library(rsconnect)
library(bslib)
############# Epi Curve Processing Start ###############################

# library(plotly)
df <- read.csv("https://raw.githubusercontent.com/plotly/datasets/master/2011_us_ag_exports.csv")
dfMal <- read.csv("data/malaria_initial_meta.csv")

plotlyFunTest <- function(df){
    # df$hover <- with(df, paste(Submitter_State, '<br>', "Count", n))
    stateCount <- df %>% count(
        Submitter_State
    )
    # stateCount$hover <- with(stateCount, paste(Submitter_State, '<br>', "Count: ", n, sep = ""))
    stateCount$hover <- with(stateCount, paste(Submitter_State, '<br>', "Count: ", n, sep = ""))

    # give state boundaries a white border
    l <- list(color = toRGB("white"), width = 2)
    # specify some map projection/options
    g <- list(
    scope = 'usa',
    projection = list(type = 'albers usa'),
    showlakes = TRUE,
    lakecolor = toRGB('white')
    )

    fig <- plot_geo(stateCount, locationmode = 'USA-states')
    fig <- fig %>% add_trace(
        z = ~n, hovertext = ~hover, locations = ~Submitter_State,
        color = ~n, 
        colors = 'Purples',
        hoverinfo = "text"
    )
    fig <- fig %>% colorbar(title = "Count")
    fig <- fig %>% layout(
        title = 'Specimens Submitted for Malaria Genomic Surveillance',
        geo = g
    )

    fig 
}


plotlyWorld <- function(df){
    df <- df[!df$Reported_Travel_History == "",]
    countryCount <- df %>% count(
        Reported_Travel_History
    )
    countryCount <- as.data.frame(countryCount)
    # print(countryCount)
   countryCount$code <- countrycode(sourcevar = countryCount$Reported_Travel_History, destination = "iso3c", origin = "country.name.en")
    # print(countryCount)
    countryCount$n = as.numeric(countryCount$n)
    # countryCount <- coun
    # light grey boundaries
    l <- list(color = toRGB("grey"), width = 0.5)

    # specify map projection/options
    g <- list(
    showframe = TRUE,
    showcoastlines = TRUE,
    projection = list(type = 'natural earth')
    )

    fig <- plot_geo(countryCount)
    fig <- fig %>% add_trace(
        z = ~n, color = ~n, colors = 'Blues',
        text = ~Reported_Travel_History, locations = ~code, marker = list(line = l)
    )
    fig <- fig %>% colorbar(title = 'Countries', tickprefix = '$')
    fig <- fig %>% layout(
        # title = '2014 Global GDP<br>Source:<a href="https://www.cia.gov/library/publications/the-world-factbook/fields/2195.html">CIA World Factbook</a>',
        geo = g
    )

    fig
    }

countSequences <- function(df){
    allVec <- unique(df$LSDB_Sequence_ID)
    allVec <- allVec[!allVec == ""] 
    val <- length(allVec)
    # as.character(val)

    fig <- plot_ly(
        mode = "number",
        type = "indicator",
        value = val,
        title = list(text = "Specimens Sequenced", font = list(size = 26)),
        height = 300
    )
    fig <- fig %>% 
        layout(paper_bgcolor = "#bfbfbf")
    fig


}
countSubmissions <- function(df){
    val <- length(unique(df$LSDB_Specimen_ID))
    # as.character(val)

    fig <- plot_ly(
        mode = "number",
        type = "indicator",
        value = val,
        title = list(text = "Specimens Submitted", font = list(size = 26)),
        height = 300
    )
    fig <- fig %>% 
        layout(paper_bgcolor = "#e6e6e6")
    fig
    
}
countStates<- function(df){
    val <- length(unique(df$Submitter_State))
    # as.character(val)
    fig <- plot_ly(
        mode = "number",
        type = "indicator",
        value = val,
        title = list(text = "Submitter States", font = list(size = 26)),
        height = 300
    )
    fig <- fig %>% 
        layout(paper_bgcolor = "#bfbfbf")
    fig
    
    
    fig
}




vars <- setdiff(names(iris), "Species")

ui <- pageWithSidebar(
  headerPanel('Iris k-means clustering'),
  sidebarPanel(
    selectInput('xcol', 'XTEST Variable', vars),
    selectInput('ycol', 'Y Variable', vars, selected = vars[[2]]),
    numericInput('clusters', 'Cluster count', 3, min = 1, max = 9)
  ),
  mainPanel(
    plotOutput('plot1'),
    plotlyOutput('plotlyMap'),

  )
)


ui2 <- fluidPage(theme = shinytheme("cerulean"),
  navbarPage("Malaria Surveillance Dashboard",
	tabPanel("Specimen Submission Summary",

	sidebarLayout(
		sidebarPanel(style = "position:fixed;width:30.5%;",
			tabsetPanel(
				tabPanel( 
					helpText("Overview"),
					    strong(h4("Overview")),
                            p("The visualizations in this document are meant to supplement epidemiological investigations of ", em("Cyclospora cayetanensis."), "The data represented here 
                                have been collated from various outputs of the ",  em("Cyclospora cayetanensis"), " bioinformatic workflow. This means the data represented here only come from 
                                specimens, or sequences, that have been submitted to CDC for genotyping. In most rows, users can choose from different options under the", 
                                strong('Choose Variables to Display'), "tab to display different maps/graphs for the given variable. For example, in Row C, users can choose different Temporal Genetic Clusters (TGCs)
                                from the dropdown menu, and the map and plot in Row C will be updated to reflect abundance/collection date* information for the selected TGC. The final row of the page 
                                (Row G) allows users to select any of the plots from Rows A - D and view the plot as a larger image. Please refer to the", strong('Instructions (and tips) for Use'), "section below
                                for more information on the how to get the most out of this page. Pleae refer to the", strong('Legends'), "tab for more detailed information on each plot."),
                            )
                        )
                    ),
                     mainPanel(
                        fluidRow(h1("Specimen Submission Summary"),splitLayout(plotlyOutput('statesBox'), plotlyOutput('submissionsBox'), plotlyOutput('sequencedBox'))),
                        fluidRow(h1("Specimen Submission Map"), plotlyOutput('plotlyMap')),
                        fluidRow(h1("Travel History Map"), plotlyOutput('plotlyWorldMap'))
                        
                        )
        
        ) 
    ),

    tabPanel("Specimen Lookup"

    ),
    tabPanel("P. vivax Genotyping Results"
    ),
    tabPanel("P. falciparum Genotyping Results"
    )
  )
)

server <- function(input, output, session) {

  # Combine the selected variables into a new data frame
#   selectedData <- reactive({
#     iris[, c(input$xcol, input$ycol)]
#   })

#   clusters <- reactive({
#     kmeans(selectedData(), input$clusters)
#   })

#   output$plot1 <- renderPlot({
#     palette(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
#       "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"))

#     par(mar = c(5.1, 4.1, 0, 1))
#     plot(selectedData(),
#          col = clusters()$cluster,
#          pch = 20, cex = 3)
#     points(clusters()$centers, pch = 4, cex = 4, lwd = 4)
#   })
    output$plotlyMap <- renderPlotly({
		# newMap <- input$myTGC
		plotlyFunTest(dfMal)
		})

    output$plotlyWorldMap <- renderPlotly({
		# newMap <- input$myTGC
		plotlyWorld(dfMal)
		})
    output$statesBox <- renderPlotly({
        countStates(dfMal)
    })
    output$submissionsBox <- renderPlotly({
        countSubmissions(dfMal)
    })
    output$sequencedBox <- renderPlotly({
        countSequences(dfMal)
    })

}
# shinyApp(ui = ui, server = server)
shinyApp(ui = ui2, server = server)
# 					p("*All dates represented on the graph are dates of specimen collection provided by State Public Health Labs."),
# 					p("Please take special note of how the axis value change in different plots when different TGCs, Clusters, etc. are chosen. Each cluster has a different size and the plots automatically resize to 
# 						fit the data.")),
# 					tabPanel(
# 						helpText("Instructions"),
# 					strong(h4("Instructions (and tips) for Use")),
# 					p("Plots are organized by row, from Row A to Row F, and the plots in Rows C through F are reactive to user input, while the plots in Rows A and B are static. The options for display on each plot 
# 						can be seen in the various menus on the", strong('Choose Variables to Display'), "tab."),
# 					p("In each plot, users can hover over a portion of the plot to retrieve additional information about that segment of the graph."),
# 					p("On the epidemiological curves, users have two options for highlighting segments of interest. First, users can click on a single color bar in the plot itself, this will enhance all bar segments with this color (Genetic Cluster, TGC, Epi Cluster, or State depending on the plot)
# 						and bar segments from other colors will be greyed out (but still visible). Users can only highlight a single category this way and users can exit this view by double clicking anywhere in the white space on the plot.
# 						 Second, users can double click on the category of interest in the legend. This will hide all other categories from the plot; however, additional categories can be seen when single clicking on that category in the legend.
# 						 All other categories can be brough back by double clicking on the category that the user originally clicked on. Alternatively, users can single click on an element in the legend, and this will make that element disappear 
# 						 from the plot - the element can be brought back by clicking on that element in the legend."),
# 					p("Each plot has a control bar that users will see upon moving the cursor the top right of the plot. The save image option will save the plot exactly how the user sees the plot at the 
# 						moment the button was clicked - the downloaded PNG file will be called 'newplot.png'. The pan, zoom in, zoom out, and reset axes are best used in together. Users can zoom in or out by clicking the respective button, 
# 						then select the pan button to move around the plot, then finally press the home button to zoom back to the original view. The final option in the control bar is an option to control hover behavior; either to
# 						hover on a single element that the cursor is on, or to hover on multiple items near by (which will compare hover text). The default is to hover on a single item.")),
# 				tabPanel(
# 					helpText("Choose Variables to Display"),
# 					h5(strong("Rows A and B: Plots are Static")),
# 					p("No input choice"),
# 					selectInput("myTGC",label = "Row C: Choose a TGC to Display",choices = currentClus_tgc, selected = "TGC_2021_001"),
# 					selectInput("myEpi",label = "Rows D and E: Choose an Epi Cluster to Display",choices = currentClus_epi, selected = "Epi_Cluster_2021 August Butter 1"),
# 					# selectInput("graph", label = "Choose a plot to view in detail:", choices = c("Row A: Collection Date by TGC", "Row A: Collection Date by Epi Cluster", "Row B: TGC Membership by State", "Row B: TGC Membership Epi Curve")),
# 					checkboxGroupInput("myHaps",label = "Row F: Choose TGC Haplotypes to Compare (2 or more)", choices = choicesBox, inline = T),
# 					radioButtons("graph", label = "Row G: Choose a plot to view in detail:",
# 						choices = c("Row A: Epi Curve by Genetic Cluster", "Row B: Epi Curve by TGC", "Row B: Epi Curve by Epi Cluster", "Row C: TGC Membership by State", "Row C: TGC Membership Sunburst Plot",
# 							"Row D: Epi Cluster Membership by State", "Row D: Epi Cluster Epi Curve"))),
# 				tabPanel(
# 					helpText("Legends"),
# 					strong(h3("Legends for Each Row")),
# 					strong(h4("Row A")),
# 					p("Epidemiological curve for current season specimens and sequences sent for", em("Cyclospora cayetanensis."), 
# 						"The color categories in the legend represent genetic clusters where specimens from the current season clustered."),
# 					strong(h4("Row B")),
# 					p("Epidemiological curves for all specimens and sequences sent for", em("Cyclospora cayetanensis."),
# 						"First plot: The color categories in the legend represent temporal genetic clusters (TGCs) where specimens from the current season clustered. 
# 						Second plot: The color categories in the legend represent epidemiological cluster where specimens from the current season clustered."),
# 					strong(h4("Row C")),
# 					p("First plot: Map of U.S.A., colored by TGC abundance. Second plot: Epi linkages within selected TGC"),
# 					strong(h4("Row D")),
# 					p("First plot: Map of U.S.A., colored by Epi Cluster abundance. Second plot: Epi curve for each state with specimens in this Epi Cluster"),
# 					strong(h4("Row E")),
# 					p("Sunburst plot with proportional cluster and TGC membership of specimens in the Epi Cluster"),
# 					strong(h4("Row F")),
# 					p("Haplotype barcode for selected TGCs"),
# 					strong(h4("Row G")),
# 					p("Larger view of selected plot"))
				
# 				# 	helpText("Legends"),
# 				# 	strong(h3("Legends for Each Row")),
# 				# 	strong(h4("Row A")),
# 				# 	p("legend for Row A"),
# 				# 	strong(h4("Row B")),
# 				# 	p("legend for Row B"),
# 				# 	strong(h4("Row C")),
# 				# 	p("legend row C"),
# 				# 	strong(h4("Row D")),
# 				# 	p("legend row D"),
# 				# 	strong(h4("Row E")),
# 				# 	p("Legend Row E")))

			
# 		),
# 			tags$head(
# 			    tags$style(
# 			      HTML(
# 			        ".checkbox-inline { 
# 			                    margin-left: 0px;
# 			                    margin-right: 10px;
# 			          }
# 			         .checkbox-inline+.checkbox-inline {
# 			                    margin-left: 0px;
# 			                    margin-right: 10px;
# 			          }
# 			        "
# 			      )
# 			    ) 
# 			  )
			
# ),
# 		mainPanel( 
# 			fluidRow(h1("Row A"), plotlyOutput("allGenetics", height = "450px", width= "95%")),
# 			fluidRow(h1("Row B"), splitLayout(plotlyOutput("fullTGCs", width = "100%", height = "600px"),plotlyOutput("epiCurve", width = "98%", height = "600px"))),
# 			fluidRow(h1("Row C"), splitLayout(plotlyOutput("plotlyMap", width = "100%", height = "500px"),plotlyOutput("plotlySunburst_TGC", width = "98%", height = "500px"))),
# 			fluidRow(h1("Row D"), splitLayout(plotlyOutput("plotlyMap2", width = "100%", height = "500px"),plotlyOutput("epiCurve_ind", width = "98%", height = "500px"))),
# 			# fluidRow(h1("Row D"), splitLayout(plotlyOutput("plotlySunburst"), plotlyOutput("plotlyBarcode", height = "500px", width = "1200px"))),
# 			fluidRow(h1("Row E", plotlyOutput("plotlySunburst"))),
# 			fluidRow(h1("Row F", plotlyOutput("plotlyBarcode", height = "650px", width = "95%"))),
# 			fluidRow(h1("Row G", plotlyOutput("selected_graph", height = "750px", width = "95%"))),
# 		),

# 		)),
# #   tabPanel("Genetic Strains with Epi History",
# #            sidebarLayout(
# #              sidebarPanel("Test",
# # 			 	tabsetPanel(
# # 					 tabPanel(
# # 						 helpText("Genetic Strains Background"),
# # 						 p("This page displays how current current genetic strains map to epidemiological clusters from previous years. In the first table, each row is a genetic strain and each column is an epidemiological cluster from 2018-2021. The numbers represent the frequency of the strain in the epidemiological cluster from samples in our refernece population. For this reason, the values are not the total number of specimens in the epidemiological cluster, they are only meant to represen the primary strain for each epidemiological cluster"),

# # 						 checkboxGroupInput("myStrainHaps",label = "Row B: Choose Strain Haplotypes to Compare (2 or more)", choices = strainChoices, inline = T),
# # 						 downloadButton("downloadStrainEpi", "Download Strain + Epi History")

# # 					 )
# # 				 )),
# #            mainPanel(width = 8,
# # 		   h1("Strain Names Effective: ", strainEffectiveDate),
# #         #    fluidRow(h2("Row A"),DT::dataTableOutput("strainHistory")),
# # 		   fluidRow(h2("Row B", plotlyOutput("strainBarcode", height = "650px", width = "95%"))),
# #            )),
# # 	),
#   tabPanel("Genetic Strain Tables",
#   	sidebarLayout(
# 		  sidebarPanel("Test",
# 		  	tabsetPanel(
# 				  tabPanel(
# 					  helpText("Background"),
# 					#   p("Descript"),
# 					  p("This page displays how current current genetic strains map to epidemiological clusters from previous years. In the first table, each row is a genetic strain and each column is an epidemiological cluster from 2018-2021. The numbers represent the frequency of the strain in the epidemiological cluster from samples in our refernece population. For this reason, the values are not the total number of specimens in the epidemiological cluster, they are only meant to represen the primary strain for each epidemiological cluster")

# 					#   checkboxGroupInput("myStrainHaps",label = "Row B: Choose Strain Haplotypes to Compare (2 or more)", choices = strainChoices, inline = T)
# 					  ),
# 				  tabPanel(
# 					  helpText("Download CSV Files"),
# 					  p('Click on the buttons to download a csv file of each table'),
# 					  downloadButton("downloadStrainEpi", "Download Strain + Epi History"),
# 					  br(),
# 					  br(),
# 					  downloadButton("downloadStrainNames", "Download Strain Names"),
# 					  br(),
# 					  br(),
# 					  downloadButton("downloadMtKey", "Download Mitochondrial Key File"),
# 					  br(),
# 					  br(),
# 					  downloadButton("downloadNuKey", "Download Nuclear Key File")
# 				  ),
# 				  tabPanel(
# 					  helpText("Strain Barcode Choices"),
# 					  p("Choose Strain Haplotypes to Compare"),
# 					  checkboxGroupInput("myStrainHaps",label = "Select 2 or more", choices = strainChoices, inline = T)
# 					  )
# 			  )),
# 			mainPanel(width = 8,
# 				h1("Strain Names Effective: ", strainEffectiveDate),
# 				tabsetPanel(id = 'dataset',
# 				tabPanel("Strains + Historical Epi Links", DT::dataTableOutput("strainHistory")),
# 				tabPanel("Strain Names",DT::dataTableOutput("strainNames")),
# 				tabPanel("Mitochondrial Key",DT::dataTableOutput("mtKey")),
# 				tabPanel("Nuclear Key", DT::dataTableOutput("nuKey"))
				
				
# 				),
# 			fluidRow(h2("Row B", plotlyOutput("strainBarcode", height = "650px", width = "95%")))
# 	  )),
#   )
# #   tabPanel("Nuclear Key",
# #   	sidebarLayout(
# # 		  sidebarPanel("Test",
# # 		  	tabsetPanel(
# # 				  tabPanel(
# # 					  helpText("test1"),
# # 					  p('Desciption'),
# # 					  downloadButton("downloadNuKey", "Download Nuclear Key File")
# # 				  )
# # 			  )),
# # 			mainPanel(width = 8,
# # 				fluidRow(DT::dataTableOutput("nuKey"))
# # 	  )),
# #   )
# )
# )



### Epi curve processing start ### 

# inFile <- 

# masterPath <- list.files(path = "data", pattern = "*CYCLOSPORA_MASTER_SHEET.xlsx")
# print(masterPath)
# masterPath <- paste("data/",masterPath, sep = "")
# epiData <- read_excel(masterPath[1], sheet = "Sheet1")
# print("firstFile")

# file1 <- list.files(path = "data", pattern = "*RESULTING_CLUSTERS_clus*")
# file1 <- paste("data/",file1, sep  = "")
# clusFile <- read.table(file1[1], header = T, sep = "\t")
# print("secondFile")
# file2 <- list.files(path = "data", pattern = "*TGC_MEMBERSHIPS.txt*")
# file2 <- paste("data/",file2, sep = "")
# tgcFile <- read.table(file2[1], header = T, sep = "\t")
# print("thirdFile")
# barcodePath <- list.files(path = "data", pattern = "*_TGC_haplotypeSheet.txt")
# barcodePath <- paste("data/",barcodePath, sep = "")
# tgcBarcodeData <- read.table(barcodePath[1], header = T, sep = "\t")

# strainNamesPath <- list.files(path = "data", pattern = "*_Strain_haplotypeSheet.txt")
# strainNamesPath <- paste("data/", strainNamesPath, sep = "")
# strainBarcodeData <- read.table(strainNamesPath[1], header = T, sep ="\t")

# epiStrainsPath <- list.files(path = "data", pattern = "2022_startSeason_strainName_epiLinks.txt")
# epiStrainsPath <- paste("data/", epiStrainsPath, sep = "")
# epiStrains <- read.table(epiStrainsPath[1], header = T, sep ="\t")

# mtKeyPath <- list.files(path = "data", pattern = "2022_startingStrainName_97_stringency_Hmat_Mt_clean2.csv")
# mtKeyPath <- paste("data/", mtKeyPath, sep = "")
# max(count.fields(mtKeyPath[1])) -> maxFields_Mt
# mtKeyFile <- read.table(mtKeyPath[1],  sep = ",", col.names = paste0("V",seq_len(maxFields_Mt+1)), fill = T, header = F)

# nuKeyPath <- list.files(path = "data", pattern = "2022_startingStrainName_97_stringency_Hmat_Nu_clean2.csv")
# nuKeyPath <- paste("data/", nuKeyPath, sep = "")

# max(count.fields(nuKeyPath[1])) -> maxFields_nu
# nuKeyFile <- read.table(nuKeyPath[1],  sep = ",", col.names = paste0("V",seq_len(maxFields_nu+1)), fill = T, header = F)

# # epiStrains <- read.table("data/2022_inSeason_strainName_epiLinks.txtt", header = T, sep = "\t")
# # epiStrains <- read.table("data/2022_startSeason_strainName_epiLinks.txt", header = T, sep = "\t")
# names(epiStrains)[names(epiStrains) == "firstOption"] <- "Strain"

# #switch strain and epi cluster and hash out the reordering if you want to transpose the rows/columns for the epi/strain table
# maditr::dcast(epiStrains,Strain~Epi_Cluster_Name) -> epiStrains_count
# epiStrains_countReorder <- epiStrains_count[,c(1,5,9,13,19,20,22,2,3,4,6,15,18,16,17,10,11,12,14,21,7,8)]
# # epiStrains_count -> epiStrains_countReorder

# sort(unique(epiStrains$Strain)) -> latestStrainNames
# latestStrainNames <- epiStrains[match(latestStrainNames,epiStrains$Strain),]
# latestStrainNames <- latestStrainNames[,c(3,4)]


# # max(count.fields("data/mtStrains_clean.csv")) -> maxFields_Mt
# # mtKeyFile <- read.table("data/mtStrains_clean.csv", sep = ",", col.names = paste0("V",seq_len(maxFields_Mt+1)), fill = T, header = F)

# # max(count.fields("data/2022_startingStrainName_97_stringency_Hmat_Mt_clean2.csv")) -> maxFields_Mt
# # mtKeyFile <- read.table("data/2022_startingStrainName_97_stringency_Hmat_Mt_clean2.csv", sep = ",", col.names = paste0("V",seq_len(maxFields_Mt+1)), fill = T, header = F)

# # max(count.fields("data/nuStrains_clean.csv")) -> maxFields_nu
# # nuKeyFile <- read.table("data/nuStrains_clean.csv", sep = ",", col.names = paste0("V",seq_len(maxFields_nu+1)), fill = T, header = F)

# # max(count.fields("data/2022_startingStrainName_97_stringency_Hmat_Nu_clean2.csv")) -> maxFields_nu
# # nuKeyFile <- read.table("data/2022-03-31-1026_updated_NuKey_clean.csv", sep = ",", col.names = paste0("V",seq_len(maxFields_nu+1)), fill = T, header = F)

# # effectiveDate_path <- list.files(path = "data", pattern = "2022*inSeason_clusterNaming.csv")

# effectiveDate_path <- list.files(path = "data", pattern = "*start2022_clusterNaming.csv")
# print(effectiveDate_path)
# strainEffectiveDate <- strsplit(effectiveDate_path[1], split = "_")[[1]][1]

# epiCurve_data <- data.frame(matrix(nrow = length(epiData$Seq_ID), ncol = 4))
# epiCurve_data$Seq_ID <- epiData$Seq_ID 
# epiCurve_data$Collection_date <- epiData$Collection_date
# epiCurve_data$Current_Cluster_Code <- epiData$Current_Cluster_Code
# epiCurve_data$Epi_Cluster_Name <- epiData$Epi_Cluster_Name

# epiCurve_data <- merge(epiCurve_data, clusFile, by = "Seq_ID", all.x =T)

# epiData2 <- merge(epiCurve_data, tgcFile, by = "Seq_ID", all.x = T)

# epiData2$Most_Recent_Cluster_Code[is.na(epiData2$Most_Recent_Cluster_Code)] <- ""
# epiData2$Current_Cluster_Code[is.na(epiData2$Current_Cluster_Code)] <- ""
# epiData2$Most_Recent_Cluster_Code[epiData2$Most_Recent_Cluster_Code == "None"] <- ""
# epiData2$Current_Cluster_Code[epiData2$Current_Cluster_Code == "Dissolved"] <- ""
# epiData2 <- unite(epiData2, midTGC, c(Most_Recent_Cluster_Code, Current_Cluster_Code), sep ="")

# epiData2$finalTGC <- substr(epiData2$midTGC,1,8)
# epiData2$finalTGC[epiData2$finalTGC == ""] <- "none"
# # print(epiData2)
# epiCurve_data <- epiData2

# myVec <- vector()
# for(i in epiCurve_data$Cluster){
# 	# k <- as.character(i)
# 	# print(nchar(k))
# 	if(nchar(i) > 1 || is.na(i)){
# 		myVec <- c(myVec, i)
# 	}
# 	else {
# 		k <- paste("0", i, sep = "")
# 		myVec <- c(myVec, k)
# 		# print("short length")
# 	}
# }

# print("after clus")

# epiCurve_data$newCluster <- myVec
# epiCurve_data$Cluster <- paste("Cluster_", epiCurve_data$newCluster, sep = "")
# # print(epiCurve_data)
# #format state and date information.
# substr(epiCurve_data$Seq_ID, 2,3) -> epiCurve_data$State
# as.Date(epiCurve_data$Collection_date, format = "%m/%d/%y") -> epiCurve_data$dateFix

# #either keep all specimens regardless if they have a TGC ...
# # epiCurve_data$Current_Cluster_Code[epiCurve_data$Current_Cluster_Code == ""] <- "none"
# #Or remove specimens from collection curve if they don't have a TGC
# # epiCurve_data[!epiCurve_data$Current_Cluster_Code == "",] -> epiCurve_TGCs

# epiCurve_TGCs <- epiCurve_data
# paste("TGC_", epiCurve_TGCs$finalTGC, sep = "") -> epiCurve_TGCs$tgcName
# print("after name")
# epiCurve_TGCs[!grepl("DUP", epiCurve_TGCs$tgcName),] -> epiCurve_TGCs
# epiCurve_TGCs[!grepl("Dissolved", epiCurve_TGCs$tgcName),] -> epiCurve_TGCs
# epiCurve_TGCs[!epiCurve_TGCs$Epi_Cluster_Name == "",] -> epiCurve_epi
# epiCurve_epi[!epiCurve_epi$Epi_Cluster_Name == "n/a",] -> epiCurve_epi

# #to keep all specimens in epi curve plot
# epiCurve_data -> epiCurve_epi
# epiCurve_epi$Epi_Cluster_Name[epiCurve_epi$Epi_Cluster_Name== ""] <- "none"
# epiCurve_epi$Epi_Cluster_Name[epiCurve_epi$Epi_Cluster_Name == "n/a"] <- "none"
# epiCurve_epi$Epi_Cluster_Name[is.na(epiCurve_epi$Epi_Cluster_Name)] <- "none"

# paste("Epi_Cluster_", epiCurve_epi$Epi_Cluster_Name, sep = "") -> epiCurve_epi$epiRename
# epiCurve_epi$epiRename <- gsub("[()]", "_", epiCurve_epi$epiRename)


# ### Epi curve processing end ### 

# ### Epi curve function start ### 

# allTGCsPlot <- function(myData, var1){
# 	#get assortment of colors for this plot
# 	colLen <- length(unique(myData$tgcName))
# 	myCols <- sample(col_vector, colLen)
# 	#set light grey as color for specimens not in TGC
# 	myCols[colLen] <- "#d3d3d3"

# 	orderColorDF2 <- colorDF2[order(colorDF2$TGC),]
# 	orderColorDF2$Color -> colorVec
# 	colorVec[colLen] <- "#d3d3d3"

# 	#allow for highlighting by tgc name
# 	d <- highlight_key(myData, ~tgcName)

# 	#create histogram chart in ggplot with date along x axis and bars colored by tgc
# 	newGG <- ggplot(data = d, aes(x = dateFix, fill = tgcName, text = paste("Date:" ,dateFix))) + geom_bar() + scale_fill_manual(values = colorVec) + scale_x_date(date_breaks = "1 week", limits = as.Date(c("2021-05-01", "2021-09-30"))) +
# 	theme_bw() + theme(axis.text.x = element_text(angle = 25, size = 8), plot.title = element_text(hjust = 0.5)) + xlab ("Collection Date") + ggtitle("Epi Curve by TGC") 

# 	#make ggplot object into an interactive object with plotly
# 	gg <- ggplotly(newGG, tooltip = c("fill", "text")) %>% layout(legend = list(orientation = "h", y = -0.4, bordercolor = "Black", x= 0.075 ,font = list(size = 10)))
# 	highlight(gg, on = "plotly_click", off = "plotly_doubleclick")

# }

# allGeneticPlot <- function(myData, var1){
# 	# print(myData)
# 	#get assortment of colors for this plot
# 	colLen <- length(unique(myData$Cluster))
# 	# myCols <- sample(col_vector, colLen)
# 	myCols <- col_vector[1:colLen]
# 	#set light grey as color for specimens not in TGC
# 	# print(colLen)
# 	# print(length(myCols))
# 	myCols[colLen - 1] <- "#d3d3d3"
# 	# print(myCols)
# 	# print(myCols)
# 	# print(order(colorDF$Cluster))
# 	orderColorDF <- colorDF[order(colorDF$Cluster),]
# 	orderColorDF$Color -> colorVec
# 	colorVec[colLen - 1] <- "#d3d3d3"
# 	# print(orderColorDF)
# 	# print(orderColorDF$Color)
# 	#allow for highlighting by tgc name
# 	d <- highlight_key(myData, ~Cluster)

# 	#create histogram chart in ggplot with date along x axis and bars colored by tgc
# 	newGG <- ggplot(data = d, aes(x = dateFix, fill = Cluster, text = paste("Date:" ,dateFix))) + geom_bar() + scale_fill_manual(values = colorVec) + scale_x_date(date_breaks = "1 week", limits = as.Date(c("2021-05-01", "2021-09-30"))) +
# 	theme_bw() + theme(axis.text.x = element_text(angle = 25, size = 8), plot.title = element_text(hjust = 0.5)) + xlab ("Collection Date") + ggtitle("Epi Curve by Genetic Cluster") 

# 	#make ggplot object into an interactive object with plotly
# 	gg <- ggplotly(newGG, tooltip = c("fill", "text")) %>% layout(legend = list(orientation = "h", y = -0.50, bordercolor = "Black", x= 0.075 ,font = list(size = 10)))
# 	highlight(gg, on = "plotly_click", off = "plotly_doubleclick")

# } 

# allEpiPlot <- function(myData, var1){
# 	#get assortment of colors for this plot
# 	colLen <- length(unique(myData$Epi_Cluster_Name))
# 	myCols <- sample(col_vector, colLen)
# 	#set grey as color for specimens not in an epi cluster
# 	myCols[colLen] <- "#d3d3d3"


# 	#allow for highlighting by epi cluster name
# 	d <- highlight_key(myData, ~Epi_Cluster_Name)

# 	#create histogram chart in ggplot with date along x axis and bars colored by tgc
# 	newGG <- ggplot(data = d, aes(x = dateFix, fill = Epi_Cluster_Name, text = paste("Date:" ,dateFix))) + geom_bar() + scale_fill_manual(values = myCols) + scale_x_date(date_breaks = "1 week", limits = as.Date(c("2021-05-01", "2021-09-30"))) +
# 	theme_bw() + theme(axis.text.x = element_text(angle = 25, size = 8), plot.title = element_text(hjust = 0.5)) + xlab("Collection Date") + ggtitle("Epi Curve by Epi Cluster") 

# 	#make ggplot object into an interactive object with plotly
# 	gg <- ggplotly(newGG, tooltip = c("fill", "text")) %>% layout(legend = list(orientation = "h", y = -0.3, x= 0.075, bordercolor = "Black"))
# 	highlight(gg, on = "plotly_click", off = "plotly_doubleclick")
# }


# singleTGCPlot <- function(myData, var1){
# 	#get assortment of colors for this plot
# 	colLen <- length(unique(myData$State))
# 	myCols <- sample(col_vector, colLen)

# 	#subset data by a single TGC cluster, selected in rshiny
# 	singleTGC <- myData[myData$tgcName == var1,]

# 	#allow for highlighting by state
# 	d <- highlight_key(singleTGC, ~State)

# 	#create histogram chart in ggplot with date along x axis and bars colored by tgc
# 	newGG <- ggplot(data = d, aes(x = dateFix, fill = State, text = paste("Date:" ,dateFix))) + geom_bar() + scale_fill_manual(values = myCols) + scale_x_date(date_breaks = "1 week", limits = as.Date(c("2021-05-01", "2021-09-30"))) +
# 	theme_bw() + theme(axis.text.x = element_text(angle = 25, size = 8), plot.title = element_text(hjust = 0.5)) + xlab("Collection Date") + ggtitle(paste(var1, "TGC Membership Epi Curve")) 

# 	#make ggplot object into an interactive object with plotly
# 	gg <- ggplotly(newGG, tooltip = c("fill", "text")) %>% layout(legend = list(orientation = "h", bordercolor = "Black", y = -0.22, x = 0.075, xanchor = "middle"))
# 	highlight(gg, on = "plotly_selected", off = "plotly_doubleclick")
	
# }

# singleEpiPlot <- function(myData, var1){
# 	#get assortment of colors for this plot
# 	colLen <- length(unique(myData$State))
# 	myCols <- col_vector[1:colLen]
# 	# myCols <- sample(col_vector, colLen)

# 	#subset data by a single epi cluster, selected in rshiny
# 	singleEpi <- myData[myData$epiRename == var1,]

# 	#allow for highlighting by state
# 	d <- highlight_key(singleEpi, ~State)

# 	#create histogram chart in ggplot with date along x axis and bars colored by tgc
# 	newGG <- ggplot(data = d, aes(x = dateFix, fill = State, text = paste("Date:" ,dateFix))) + geom_bar() + scale_fill_manual(values = myCols) + scale_x_date(date_breaks = "1 week", limits = as.Date(c("2021-05-01", "2021-09-30"))) +
# 	theme_bw() + theme(axis.text.x = element_text(angle = 25, size = 8), plot.title = element_text(hjust = 0.5)) + xlab("Collection Date") +  ggtitle(paste(var1, "Epi Cluster Epi Curve")) 

# 	#make ggplot object into an interactive object with plotly
# 	gg <- ggplotly(newGG, tooltip = c("fill", "text")) %>% layout(legend = list(orientation = "h", y = -0.22, x = 0.075, bordercolor = "Black"))
# 	highlight(gg, on = "plotly_selected", off = "plotly_doubleclick")
# }


# ### Epi curve function end ### 

# ############# Epi Curve Processing End ###############################

# ############# Map Processing Start ###############################

# ### Map pre-processing begin ###

# #will need shell script to place the latest of each of these files into the correct folder
# # file1 <- list.files(path = "data", pattern = "*RESULTING_CLUSTERS_clus*")
# # file1 <- paste("data/",file1, sep  = "")
# # file2 <- list.files(path = "data", pattern = "*TGC_MEMBERSHIPS.txt*")
# # file2 <- paste("data/",file2, sep = "")
# # file3 <- list.files(path = "data", pattern = "*CYCLOSPORA_MASTER_SHEET.xlsx")
# # file3 <- paste("data/",file3, sep = "")

# #
# statesWithAKHI <- us_map(region = "states")

# # clusFile <- read.table(file1[1], header = T, sep = "\t")
# # tgcFile <- read.table(file2[1], header = T, sep = "\t")
# # # epiData <- read_excel(file3[1], sheet = "Sheet1")


# #manipulate file to get DF for cluster plot
# substr(clusFile$Seq_ID,2,3) -> clusFile$abbr
# substr(clusFile$Seq_ID,10,11) -> clusFile$year

# myVec <- vector()
# for(i in clusFile$Cluster){
# 	# k <- as.character(i)
# 	# print(nchar(k))
# 	if(nchar(i) > 1 || is.na(i)){
# 		myVec <- c(myVec, i)
# 	}
# 	else {
# 		k <- paste("0", i, sep = "")
# 		myVec <- c(myVec, k)
# 		# print("short length")
# 	}
# }

# clusFile$newCluster <- myVec
# clusFile$clusRename <- paste("Cluster_", clusFile$newCluster, sep = "")

# # paste("Cluster_", clusFile$Cluster, sep = "") -> clusFile$clusRename
# clusFile <- clusFile[clusFile$year == "21",]
# clusFile %>% maditr::dcast(abbr ~ clusRename) -> newTable
# newTable <- as.data.frame(newTable)
# # print(newTable)
# mergedStates <- left_join(statesWithAKHI, newTable, by = "abbr")
# colnames(mergedStates)[grep("Cluster",colnames(mergedStates))] -> currentClus

# #manipulate file to get DF for TGC plot
# substr(tgcFile$Seq_ID,2,3) -> tgcFile$abbr
# # print(tgcFile)
# paste("TGC_", tgcFile$Most_Recent_Cluster_Code, sep = "") -> tgcFile$tgcRename
# # paste("TGC_", tgcFile$Most_Recent_Cluster_Code, sep = "") -> tgcFile$tgcRename
# tgcFile %>% maditr::dcast(abbr ~ tgcRename) -> newTable_tgc
# # print(head(newTable_tgc))
# newTable_tgc <- as.data.frame(newTable_tgc)

# #check for this one
# # newTable_tgc <- newTable_tgc[, -which(names(newTable_tgc) %in% "TGC_NA")]


# # print(newTable_tgc)
# # newTable_tgc <- newTable_tgc[, -which(names(newTable_tgc) %in% "DUP")]
# # newTable_tgc <- newTable_tgc[, -which(names(newTable_tgc) %in% "TGC_Dissolved")]
# mergedStates_tgc <- left_join(statesWithAKHI, newTable_tgc, by = "abbr")
# colnames(mergedStates_tgc)[grep("TGC",colnames(mergedStates_tgc))] -> currentClus_tgc



# #manipulate file to get DF for epi plot
# epiData$Epi_Cluster_Name[epiData$Epi_Cluster_Name == "n/a"] <- NA
# testing <- merge(clusFile, epiData, by = "Seq_ID", all.x = T)
# fullMerge <- merge(testing, tgcFile, by = "Seq_ID", all.x = T)
# paste("Epi_Cluster_", fullMerge$Epi_Cluster_Name, sep = "") -> fullMerge$epiRename
# fullMerge$epiRename <- gsub("[()]", "_", fullMerge$epiRename)

# # fullMerge$epiRename <- str_replace(fullMerge$epiRename, "(", "_")
# # fullMerge$epiRename <- str_replace(fullMerge$epiRename, ")", "_")
# paste(fullMerge$clusRename, fullMerge$tgcRename, sep = "") -> fullMerge$clusTGC_rename
# fullMerge %>% maditr::dcast(abbr.y ~ epiRename) -> newTable_epi
# # print(newTable_epi)
# newTable_epi <- as.data.frame(newTable_epi)
# newTable_epi <- na.omit(newTable_epi)
# newTable_epi <- newTable_epi[, -which(names(newTable_epi) %in% "Epi_Cluster_NA")]
# colnames(newTable_epi)[colnames(newTable_epi) == "abbr.y"] <- "abbr"
# mergedStates_epi <- left_join(statesWithAKHI, newTable_epi, by = "abbr")
# colnames(mergedStates_epi)[grep("Epi_Cluster",colnames(mergedStates_epi))] -> currentClus_epi

# # print(fullMerge[grepl("Epi_Cluster_2021 KS Local Cluster 1", fullMerge$epiRename),])

# ### Map pre-processing End ###

# ### Map Function Start ###

# plotlyFun1 <- function(myData, var1, color, myData2){
# 	newCol <- length(colnames(myData))

# 	for(i in 1:length(myData$abbr)){
# 		myData[i,newCol + 1] <- sum(myData[i,2:newCol])
# 	}
# 	myData$abbr -> myData$code
# 	abbr2state(myData$code) -> myData$state
# 	colnames(myData)[newCol + 1] <- "stateTGCSum"
# 	sum(myData[,var1]) -> myData$completeTGCsum

# 	newCol2 <- length(colnames(myData2))
# 	for(i in 1:length(myData2$abbr)){
# 		myData2[i,newCol2 + 1] <- sum(myData2[i,2:newCol2])

# 	}
# 	colnames(myData2)[newCol2 + 1] <- "stateClusterSum"
# 	aTab <- merge(myData, myData2, by = "abbr", all.x = T)
# 	myData <- aTab
# 	myData$hover <- with(myData, paste(state, "<br>", state, " ", "Total Specimens Genotyped", ":", "\t",stateClusterSum, "<br>", state, " ", "Specimens in ", var1, ":", "\t", myData[,var1] , "<br>", "Total Specimens in ", var1,":", "\t",  completeTGCsum, sep = ""))

# 	# myData$hover <- with(myData, paste(state, '<br>', "Specimens in ", var1, ":", "\t", myData[,var1], "<br>", "Total Specimens in a TGC:", "\t",stateTGCSum, "<br>", "Total Specimens in this TGC:", "\t",completeTGCsum, "<br>", "state total specimens genotyped", "\t",stateClusterSum ,sep = ""))
# 	# var1 <- sym(var1)
# 	l <- list(color = toRGB("white"), width = 2)
# 	g <- list(
# 	  scope = 'usa',
# 	  projection = list(type = 'albers usa'),
# 	  showlakes = TRUE,
# 	  lakecolor = toRGB('white')
# 	)
# 	myData$category <- "state"
	
# 	fig <- plot_geo(myData, locationmode = 'USA-states')
# 	fig <- fig %>% add_trace(
# 	    z = ~myData[,var1], text = ~hover, locations = ~code,
# 	    color = ~myData[,var1], colors = color)
# 	fig <- fig %>% colorbar(title = "Cluster Abundance")
# 	fig <- fig %>% layout(
# 	    title = paste('TGC Membership', '<br>', var1),
# 	    geo = g
# 	  )
	
# }

# plotlyFun2 <- function(myData, var1, color, mergeFile){
# 	limData <- mergeFile[grepl(var1, mergeFile$epiRename),]
# 	# print(limData)

# 	#if we want to have the genetic cluster abundance or the TGC abundance in the hover window
# 	#dcast only prints character when there is a single value in a each state, which means that needs converting in loops below
# 	maditr::dcast(limData, abbr.y ~ clusRename) -> limTGC
# 	# print(is.character(limTGC))
# 	stateList <- colnames(limTGC)[2:length(colnames(limTGC))]
# 	myLen <- length(stateList) +1 
# 	for(j in 1:length(limTGC$abbr.y)){
# 		# print(j)
# 		for(k in 2:myLen){
# 				value <- as.numeric(limTGC[j,..k])
# 				#this if statement is necessary if cluster only has 1 occurrance in each state (honestly not sure why it is needed, but it works)
# 				if(is.na(value) == TRUE){
# 					limTGC[j,k] <- 1

# 				}
# 			}
# 			}
# 	# dcast(limData, abbr.y ~ Current_Cluster_Code.y) -> limTGC
# 	newCol2 <- length(colnames(limTGC))
# 	limTGC$total_sum <- "NA"
# 	#necessary for clusters with only a single value from a single state, or single value in multiple states
# 	for(i in 1:length(limTGC$abbr.y)){
# 		selected <- as.matrix(limTGC[i,2:myLen])
# 		selected <- sapply(selected, as.numeric)
# 		value <- sum(selected)
# 		limTGC[i,newCol2 + 1] <- value
# 	}
	
# 	limTGC$abbr.y -> limTGC$code
# 	abbr2state(limTGC$code) -> limTGC$state
# 	limTGC$category <- "state"
# 	as.integer(limTGC$total_sum) -> limTGC$total_sum
# 	colnames(limTGC)[2:newCol2] -> tgcVector
# 	as.data.frame(limTGC) -> limTGC
# 	pp <- vector()
# 	for(j in limTGC$state){
# 		stateData <- limTGC[limTGC$state  == j,]
# 		qq <- vector()
# 		for(i in 1:newCol2+1){
# 			newVar <- with(stateData, paste("<br>", colnames(stateData)[i], ":", "\t", stateData[,i], sep = ""))
# 			# print(newVar)
# 			qq <- paste(qq, newVar)
# 	}
# 	pp <- c(pp, qq)
# 	}
# 	limTGC$hover <- with(limTGC, paste(state, pp)) 
# 	# limTGC$hover <- with(limTGC, paste(state, "<br>", colnames(limTGC)[3], limTGC[,3], colnames(limTGC)[4], limTGC[,4], colnames(limTGC)[5], limTGC[,5]))
# 	l <- list(color = toRGB("white"), width = 2)
# 	g <- list(
# 	  scope = 'usa',
# 	  projection = list(type = 'albers usa'),
# 	  showlakes = TRUE,
# 	  lakecolor = toRGB('white')
# 	)
# 	fig <- plot_geo(limTGC, locationmode = 'USA-states')
# 	fig <- fig %>% add_trace(
# 	    z = ~limTGC$total_sum, text= ~hover, locations = ~code,
# 	    color = ~limTGC$total_sum, colors = color 
#       )
# 	fig <- fig %>% colorbar(title = "Cluster Abundance")
# 	fig <- fig %>% layout(
# 	    title = paste('Epi Cluster Membership', '<br>', var1),
# 	    geo = g
# 	  )
	
# }

# ### Map Function End ###

# ############# Map Processing Start ####################################

# ############# Sunburst Processing Start ###############################

# ### Sunburst plot uses fullMerge dataframe from a previous prepocessing step

# sunburstPlot <- function(myData, var1, color){
# 	#get the epi cluster to make a sunburst plot
# 	filtData <- myData[grepl(var1, myData$epiRename),]
# 	filtData <- filtData[!grepl("NA", rownames(filtData)),]
# 	numEpi <- length(filtData$epiRename)

# 	#get the cluster and tgc info
# 	newDat2 <- filtData[,c("clusRename", "tgcRename")]

# 	#count the number of times each genetic cluster / TGC appears in this epi cluster
# 	count(newDat2, clusRename) -> clusCount
# 	count(newDat2, tgcRename) -> tgcCount

# 	#the epi cluster is the parent of each genetic cluster (at least as it is organized now)
# 	clusCount$parent <- var1

# 	#merge tgc count and newDat2 to find parent of each tgc (which genetic cluster does the TGC come from)
# 	merge(tgcCount, newDat2, by = "tgcRename", all.x =T) -> parentTest
# 	parentTest <- parentTest[!duplicated(parentTest$tgcRename),]

# 	#rename the cluster and TGC dfs to follow sunburst terminology. Then put in single dg
# 	colnames(parentTest) <- c("label", "value", "parent")
# 	colnames(clusCount) <- c("label", "value", "parent")
# 	newDat3 <- rbind(parentTest, clusCount)

# 	#initialize column to hold text for hover
# 	newDat3$hover <- "NA"

# 	#loop through each cluster/TGC and add hover text (whole loop is to add hover text)
# 	for(i in newDat3$label){
# 		#one loop for genetic clusters ...
# 		matcher <- filtData[filtData$clusRename == i,]

# 		#second loop for TGCs
# 		matcher2 <- filtData[filtData$tgcRename == i,]

# 		if(length(rownames(matcher)) > 0){
# 			#work with clusters, get the count of the cluster in each state
# 			maditr::dcast(matcher, clusRename ~ abbr.y) -> stateData
# 			stateList <- colnames(stateData)[2:length(colnames(stateData))]

# 			## qq is vector that will hold hover text
# 			qq <- vector()
# 			qq <- c(i)
# 			#loop through states 
# 			for(k in stateList){
# 				value <- as.numeric(stateData[1,..k])
# 				#this if statement is necessary if cluster only has 1 occurrance in each state (honestly not sure why it is needed, but it works)
# 				if(is.na(value) == TRUE){
# 					value <- 1
# 				}
# 				#passte togehter hover text
# 				newVar <- paste(k, " Spec. in Epi and ", i,":", "\t", value, sep = "")
# 				qq <- paste(qq, "<br>", newVar)
# 			}
# 		}
# 		else {
# 			#work with tgcs
# 			qq <- vector()
# 			qq <- c(i)
# 			maditr::dcast(matcher2, tgcRename ~ abbr.y) -> stateData2
# 			stateList2 <- colnames(stateData2)[2:length(colnames(stateData2))]
# 			for(k in stateList2){
# 				value2 <- as.numeric(stateData2[1,..k])
# 				#same loop needed as for the genetic clusters
# 				if(is.na(value2) == TRUE){
# 					value2 <- 1
# 				}
# 				newVar2 <- paste(k, " Spec. in Epi and ", i,":", "\t", value2, sep = "")
# 				qq <- paste(qq, "<br>", newVar2)
# 			} 
# 		} 
# 		#add in hover text in correct location for each label
# 		newDat3[newDat3$label == i,4] <- qq
# 	}
	
# 	colorTry <- merge(newDat3, colorDF, by.x = "label", by.y = "Cluster", all.x = T)
# 	colorTry2 <- merge(colorTry, colorDF, by.x = "parent", by.y = "Cluster", all.x = T)
# 	colorTry2[is.na(colorTry2)] <- ""
# 	colorTry3 <- unite(colorTry2, newColor, c(Color.x, Color.y), sep ="")
# 	# newDat3 <- colorTry3
# 	newDat3 <- colorTry3

# 	newDat4 <- newDat3[rev(order(colorTry3$value)),]
# 	finalCols <- unique(newDat4$newColor)
# 	# print(newDat4)
# 	# print(newDat4)


# 	#make the sunburst figure and format the hover text
# 	fig <- plot_ly(labels = newDat4$label, parents = newDat4$parent, values = newDat4$value, type = "sunburst", branchvalues = "total")

# 	fig <- fig %>% add_trace(hovertemplate = paste(newDat4$hover,"<br>", "Total Specimens in Epi Cluster:", "\t", numEpi,'<extra></extra>', sep ="")) %>% layout(title = paste(var1, " Cluster Membership ", "(n = ", numEpi,")", sep = ""), colorway = finalCols)

# 	fig
# }

# sunburstPlot_TGC <- function(myData, var1){
#   filtData <- myData[grepl(var1, myData$tgcRename),]
#   filtData <- filtData[!grepl("NA", rownames(filtData)),]
#   filtData$epiRename[filtData$epiRename == "Epi_Cluster_NA"] <- "Unknown Epi Linkage"
#   numTGC <- length(filtData$tgcRename)
  
#   #get the cluster and tgc info
#   newDat2 <- filtData[,c("epiRename")]
#   newDat2 <- as.data.frame(newDat2)
  
#   #count the number of times each genetic cluster / TGC appears in this epi cluster
#   count(newDat2, newDat2) -> epiCount
#   #count(newDat2, tgcRename) -> tgcCount
  
#   #the epi cluster is the parent of each genetic cluster (at least as it is organized now)
#   epiCount$parent <- var1
#   colnames(epiCount) <- c("label", "value", "parent")
#   epiCount$hover <- "NA"
  
#   if(length(epiCount$label) > 1){
#     for(i in 1:length(epiCount$label)){
#       writer <- paste(epiCount[i,1], ":", "\t", epiCount[i,2], sep = "\t")
#       epiCount[i,4] <- writer
#     }
#   } else{
#     writer <- paste(epiCount[1,1], ":", "\t", epiCount[1,2], sep = "\t")
#     epiCount[1,4] <- writer
#   }
#   fig <- plot_ly(labels = epiCount$label, parents = epiCount$parent, values = epiCount$value, type = "sunburst", branchvalues = "total")
  
#   fig <- fig %>% add_trace(hovertemplate = paste(epiCount$hover,"<br>", "Total Specimens in TGC:", "\t", numTGC,'<extra></extra>', sep ="")) %>% layout(title = paste(var1, " Cluster Membership ", "(n = ", numTGC,")", sep = ""))
  
#   fig
# }
# ### Sunburst function code end ###

# ############# Sunburst Processing End ###############################


# ############# Barcode Processing Start ###############################

# ### Barcode pre-processing code start ###

# ## The files read in below are outputs from pre-processing done by other scripts
# # the haplotype sheet is a haplotype sheet of consensus haplotypes for each cluster (from clusterStrainNaming.py script)
# # the tgcLink file is a file of current year specimens with seq ID in first column, TGC in second column, and genetic cluster in third column

# #The code commented out below is for genetic clusters and linking to TGCs, the current running code uses consensus barcodes from the TGCs themselves, not from the genetic clusters

# # barcodeData <- read.table("data/testHapSheet2.txt", header = T, sep = "\t")
# # tgcLink <- read.table("data/newestLink.txt", header =T, sep = "\t")

# #rename columns in link file for processing. Then process to put in order, remove specimens without TGCs, generate matrix for plotly heatmap
# # colnames(tgcLink) <- c("linkSeqID", "TGC", "Cluster")
# # tgcLink <- tgcLink[!grepl("None", tgcLink$TGC),]
# # tgcLink <- tgcLink[order(tgcLink$TGC),]
# # link2 <- merge(barcodeData, tgcLink, by.x = "Seq_ID",by.y = "Cluster", all.y = T)
# # link3 <- distinct(link2, TGC, .keep_all = T)
# # link3 <- link3 %>% drop_na(TGC)
# # link3 <- link3[rev(order(link3$TGC)),]
# # rownames(link3) <- link3$TGC
# # link3 <- link3[, -which(names(link3) %in% c("Seq_ID", "TGC", "Seq_ID.y"))]
# # as.matrix(link3) -> barcodeMat

# # the haplotype sheet is a haplotype sheet of consensus haplotypes for each cluster (from tgc_clusterHapSheet.py script)
# # tgcBarcodeData <- read.table("data/testHapSheet_tgc.txt", header = T, sep = "\t")

# # file1 <- list.files(path = "data", pattern = "*RESULTING_CLUSTERS_clus*")
# # file1 <- paste("data/",file1, sep  = "")
# # clusFile <- read.table(file1[1], header = T, sep = "\t")



# # print(tgcBarcodeData)

# tgcBarcodeData <- tgcBarcodeData[rev(order(tgcBarcodeData$Seq_ID)),]
# rownames(tgcBarcodeData) <- tgcBarcodeData$Seq_ID
# tgcBarcodeData <- tgcBarcodeData[, -which(names(tgcBarcodeData) %in% c("Seq_ID", "Most_Recent_Cluster_Code", "Status_of_specimen"))]
# as.matrix(tgcBarcodeData) -> barcodeMat
# # print(barcodeMat)

# #extract index position of each marker, used for making outline boxes and labelling the barcode output in function
# which(grepl("Mt_Cmt", colnames(barcodeMat)) == TRUE) -> firstMarker
# which(grepl("Mt_MSR", colnames(barcodeMat)) == TRUE) -> secondMarker
# which(grepl("Nu_360i2", colnames(barcodeMat)) == TRUE) -> thirdMarker
# which(grepl("Nu_378", colnames(barcodeMat)) == TRUE) -> fourthMarker
# which(grepl("Nu_CDS1", colnames(barcodeMat)) == TRUE) -> fifthMarker
# which(grepl("Nu_CDS2", colnames(barcodeMat)) == TRUE) -> sixthMarker
# which(grepl("Nu_CDS3", colnames(barcodeMat)) == TRUE) -> seventhMarker
# which(grepl("Nu_CDS4", colnames(barcodeMat)) == TRUE) -> eigthMarker

# firstX <- -0.5
# junctionEnd <- tail(firstMarker, n = 1) - 0.5
# msrEnd <- tail(secondMarker, n = 1) - 0.5
# nu360i2End <- tail(thirdMarker, n = 1) - 0.5
# nu378End <- tail(fourthMarker, n = 1) - 0.5
# cds1End <- tail(fifthMarker, n = 1) - 0.5
# cds2End <- tail(sixthMarker, n = 1) - 0.5
# cds3End <- tail(seventhMarker, n = 1) - 0.5
# cds4End <- tail(eigthMarker, n = 1) - 0.5

# #sort all TGCs to put as options in the checkbox menu
# choicesBox <- sort(rownames(barcodeMat))

# ### Barcode pre-processing code end ###

# ### Barcode function code start ###

# hapBarcodePlot <- function(myData, var1){
# 	#get dataframe of the selected TGCs to make barcode plot
# 	selectedClus <- myData[which(row.names(myData) %in% var1),] 

# 	#trying to figure out best way to place outline boxes and text. Still not perfect
# 	#plotly heatmap has no gridelines at least in R (to my knowledge)
# 	height <- length(rownames(selectedClus)) - 0.5
# 	textHeight <-- length(rownames(selectedClus))
# 	textHeight <- -0.75

# 	#make plotly heatmap adding outline boxes and text annotations as labels
# 	fig <- plot_ly(z = selectedClus, x = colnames(selectedClus), y = rownames(selectedClus), type = "heatmap", colorscale = "Greys", reversescale = T, showscale = F)
# 	fig <- fig %>% layout(shapes = list(list(type = "rect", fillcolor = "white", bordercolor = "black", opacity = 0.4, x0 = junctionEnd, x1 = msrEnd, xref = "x", y0 = -0.5, y1 = height, yref = "y"), 
# 		list(type = "rect", fillcolor = "blue", opacity = 0.4, x0 = firstX, x1 = junctionEnd, xref = "x", y0 = -0.5, y1 = height, yref = "y"), 
# 		list(type = "rect", fillcolor = "green", opacity = 0.4, x0 = msrEnd, x1 = nu360i2End, xref = "x", y0 = -0.5, y1 = height, yref = "y"),
# 		list(type = "rect", fillcolor = "yellow", opacity = 0.4, x0 = nu360i2End, x1 = nu378End, xref = "x", y0 = -0.5, y1 = height, yref = "y"),
# 		list(type = "rect", fillcolor = "pink", opacity = 0.4, x0 = nu378End, x1 = cds1End, xref = "x", y0 = -0.5, y1 = height, yref = "y"),
# 		list(type = "rect", fillcolor = "blue", opacity = 0.4, x0 = cds1End, x1 = cds2End, xref = "x", y0 = -0.5, y1 = height, yref = "y"),
# 		list(type = "rect", fillcolor = "white", opacity = 0.4, x0 = cds2End, x1 = cds3End, xref = "x", y0 = -0.5, y1 = height, yref = "y"),
# 		list(type = "rect", fillcolor = "red", opacity = 0.4, x0 = cds3End, x1 = cds4End, xref = "x", y0 = -0.5, y1 = height, yref = "y")), 
# 		annotations = list(list(text = "Mt_Junction", x = median(firstMarker), y =textHeight, showarrow =F, font = list(size = 25)),
# 		list(text = "Mt_MSR", x = median(secondMarker), y = textHeight, showarrow =F, font = list(size = 25)),
# 		list(text = "Nu_360i2", x = median(thirdMarker), y = textHeight, showarrow =F, font = list(size = 25)),
# 		list(text = "Nu_378", x = median(fourthMarker), y = textHeight, showarrow =F, font = list(size = 25)),
# 		list(text = "Nu <br> CDS1", x = median(fifthMarker)-1, y = textHeight, showarrow =F, font = list(size = 13)),
# 		list(text = "Nu <br> CDS2", x = median(sixthMarker)-1, y =textHeight, showarrow =F, font = list(size = 13)), 
# 		list(text = "Nu <br> CDS3", x = median(seventhMarker)-1, y = textHeight, showarrow =F, font = list(size = 13)),
# 		list(text = "Nu <br> CDS4", x = median(eigthMarker)-1, y = textHeight, showarrow =F, font = list(size = 13))), 
# 		xaxis = list(showticklabels = F, showgrid = T),
# 		showlegend = F, title = "Haplotype Barcode")
# 	fig
# }

# ### Barcode function code end ###

# ### Strain barcode processing start ###

# strainBarcodeData <- strainBarcodeData[rev(order(strainBarcodeData$Seq_ID)),]
# rownames(strainBarcodeData) <- strainBarcodeData$Seq_ID
# # tgcBarcodeData <- tgcBarcodeData[, -which(names(tgcBarcodeData) %in% c("Seq_ID", "Most_Recent_Cluster_Code", "Status_of_specimen"))]
# as.matrix(strainBarcodeData) -> strainBarcodeMat
# # print(barcodeMat)

# #extract index position of each marker, used for making outline boxes and labelling the barcode output in function
# which(grepl("Mt_Cmt", colnames(strainBarcodeMat)) == TRUE) -> strainfirstMarker
# which(grepl("Mt_MSR", colnames(strainBarcodeMat)) == TRUE) -> strainsecondMarker
# which(grepl("Nu_360i2", colnames(strainBarcodeMat)) == TRUE) -> strainthirdMarker
# which(grepl("Nu_378", colnames(strainBarcodeMat)) == TRUE) -> strainfourthMarker
# which(grepl("Nu_CDS1", colnames(strainBarcodeMat)) == TRUE) -> strainfifthMarker
# which(grepl("Nu_CDS2", colnames(strainBarcodeMat)) == TRUE) -> strainsixthMarker
# which(grepl("Nu_CDS3", colnames(strainBarcodeMat)) == TRUE) -> strainseventhMarker
# which(grepl("Nu_CDS4", colnames(strainBarcodeMat)) == TRUE) -> straineigthMarker

# strainfirstX <- -0.5
# strainjunctionEnd <- tail(strainfirstMarker, n = 1) - 0.5
# strainmsrEnd <- tail(strainsecondMarker, n = 1) - 0.5
# strainnu360i2End <- tail(strainthirdMarker, n = 1) - 0.5
# strainnu378End <- tail(strainfourthMarker, n = 1) - 0.5
# straincds1End <- tail(strainfifthMarker, n = 1) - 0.5
# straincds2End <- tail(strainsixthMarker, n = 1) - 0.5
# straincds3End <- tail(strainseventhMarker, n = 1) - 0.5
# straincds4End <- tail(straineigthMarker, n = 1) - 0.5

# strainChoices <- sort(rownames(strainBarcodeMat))
# print(strainChoices)

# strainBarcodePlot <- function(myData, var1){
# 	print(var1)
# 	#get dataframe of the selected Strain to make barcode plot
# 	selectedStrains <- myData[which(row.names(myData) %in% var1),] 

# 	#trying to figure out best way to place outline boxes and text. Still not perfect
# 	#plotly heatmap has no gridelines at least in R (to my knowledge)
# 	height <- length(rownames(selectedStrains)) - 0.5
# 	textHeight <-- length(rownames(selectedStrains))
# 	textHeight <- -0.75

# 	#make plotly heatmap adding outline boxes and text annotations as labels
# 	fig <- plot_ly(z = selectedStrains, x = colnames(selectedStrains), y = rownames(selectedStrains), type = "heatmap", colorscale = "Greys", reversescale = T, showscale = F)
# 	fig <- fig %>% layout(shapes = list(list(type = "rect", fillcolor = "white", bordercolor = "black", opacity = 0.2, x0 = strainjunctionEnd, x1 = strainmsrEnd, xref = "x", y0 = -0.5, y1 = height, yref = "y"), 
# 		list(type = "rect", fillcolor = "blue", opacity = 0.2, x0 = strainfirstX, x1 = strainjunctionEnd, xref = "x", y0 = -0.5, y1 = height, yref = "y"), 
# 		list(type = "rect", fillcolor = "green", opacity = 0.2, x0 = strainmsrEnd, x1 = strainnu360i2End, xref = "x", y0 = -0.5, y1 = height, yref = "y"),
# 		list(type = "rect", fillcolor = "yellow", opacity = 0.2, x0 = strainnu360i2End, x1 = strainnu378End, xref = "x", y0 = -0.5, y1 = height, yref = "y"),
# 		list(type = "rect", fillcolor = "pink", opacity = 0.2, x0 = strainnu378End, x1 = straincds1End, xref = "x", y0 = -0.5, y1 = height, yref = "y"),
# 		list(type = "rect", fillcolor = "blue", opacity = 0.2, x0 = straincds1End, x1 = straincds2End, xref = "x", y0 = -0.5, y1 = height, yref = "y"),
# 		list(type = "rect", fillcolor = "white", opacity = 0.2, x0 = straincds2End, x1 = straincds3End, xref = "x", y0 = -0.5, y1 = height, yref = "y"),
# 		list(type = "rect", fillcolor = "red", opacity = 0.2, x0 = straincds3End, x1 = straincds4End, xref = "x", y0 = -0.5, y1 = height, yref = "y")), 
# 		annotations = list(list(text = "Mt_Junction", x = median(strainfirstMarker), y =textHeight, showarrow =F, font = list(size = 25)),
# 		list(text = "Mt_MSR", x = median(strainsecondMarker), y = textHeight, showarrow =F, font = list(size = 25)),
# 		list(text = "Nu_360i2", x = median(strainthirdMarker), y = textHeight, showarrow =F, font = list(size = 25)),
# 		list(text = "Nu_378", x = median(strainfourthMarker), y = textHeight, showarrow =F, font = list(size = 25)),
# 		list(text = "Nu <br> CDS1", x = median(strainfifthMarker)-1, y = textHeight, showarrow =F, font = list(size = 13)),
# 		list(text = "Nu <br> CDS2", x = median(strainsixthMarker)-1, y =textHeight, showarrow =F, font = list(size = 13)), 
# 		list(text = "Nu <br> CDS3", x = median(strainseventhMarker)-1, y = textHeight, showarrow =F, font = list(size = 13)),
# 		list(text = "Nu <br> CDS4", x = median(straineigthMarker)-1, y = textHeight, showarrow =F, font = list(size = 13))), 
# 		xaxis = list(showticklabels = F, showgrid = T),
# 		showlegend = F, title = "Haplotype Barcode")
# 	fig
# }
# ############### Barcode Processing End ##################################

# #color vector for epi curve plots from I Want Hue
# # col_vector <- c("#e67e70","#40c85d","#c4239a","#75c046","#8d5fde","#a3ba37","#4e6fe5","#c7b03b","#9b3eb3","#439b38","#cd6ee4","#52c485","#e75dc1","#49883c","#604fb7","#e19a35","#598de7","#e95a36","#45c7c8","#de334d","#329981","#e6478e","#428f5c","#aa3d94","#557018","#a981e3","#889234","#884b9f","#9fb76e","#655eaf","#a78331","#476caa","#bb3326","#56a6d9","#af511b","#a79ade","#71661a","#d37bcb","#306a3c","#d94166","#70bd95","#b33a73","#557638","#e586bf","#277257","#e8814b","#7d5b99","#b26e34","#ce92c0","#515b24","#a35e90","#807a45","#864069","#d6a46f","#974863","#8c5e2f","#e38197","#aa5444","#ab3e4f","#9d565b")
# col_vector <- c("#a366a4", 
# "#4dc754",
# "#9159de",
# "#6eb62c",
# "#a037b4",
# "#b4c231",
# "#3666dc",
# "#dfb338",
# "#3c4dbb",
# "#42972e",
# "#c468e0",
# "#77c45f",
# "#6d4bba",
# "#8ea935",
# "#7b74ef",
# "#e08d24",
# "#618bed",
# "#e8682d",
# "#4067be",
# "#bdb64a",
# "#9549b6",
# "#3fc982",
# "#bd309a",
# "#3d9e59",
# "#e457bb",
# "#387725",
# "#e67ee3",
# "#6a821c",
# "#4e45a3",
# "#a28a23",
# "#695dbf",
# "#d0a145",
# "#9f82e9",
# "#516615",
# "#9c6ccd",
# "#9ec174",
# "#e43588",
# "#5dccad",
# "#c8361b",
# "#43ccd7",
# "#e03960",
# "#73c086",
# "#a72c72",
# "#428c5a",
# "#a94899",
# "#457636",
# "#e1629f",
# "#235e31",
# "#d94542",
# "#18a7c8",
# "#a93627",
# "#63baea",
# "#9c4813",
# "#509be0",
# "#c05e2c",
# "#3161a2",
# "#df8e4f",
# "#5e4393",
# "#799b51",
# "#785bae",
# "#6a6914",
# "#c690e1",
# "#415a1f",
# "#7e4895",
# "#bdb26f",
# "#505099",
# "#a36b25",
# "#a6a5e8",
# "#746015",
# "#7b81c3",
# "#f37b5f",
# "#2ba198",
# "#b63465",
# "#419d7a",
# "#af3c4c",
# "#1a6447",
# "#dd708c",
# "#367042",
# "#dd90c6",
# "#597236",
# "#6f64a6",
# "#908644",
# "#505b8f",
# "#dda573",
# "#4384b4",
# "#815519",
# "#8b4773",
# "#2f7b63",
# "#d96762",
# "#585c26",
# "#b66a92",
# "#766f3a",
# "#ed97a6",
# "#865630",
# "#8b4a5f",
# "#ad7e4e",
# "#964154",
# "#e68e7b",
# "#a3583c",
# "#ae6160")


# colLen <- length(unique(epiCurve_data$Cluster))
# myCols <- col_vector[1:colLen]
# colorDF <- data.frame(matrix(ncol = 2, nrow =colLen))
# # print(unique(epiCurve_data$Cluster[order(epiCurve_data$Cluster)]))
# print(choicesBox)
# colorDF[,1] <- unique(epiCurve_data$Cluster[order(epiCurve_data$Cluster)])
# colorDF[,2] <- myCols
# colnames(colorDF) <- c("Cluster", "Color")
# colorDF$Color[colorDF$Cluster == "Cluster_NA"] <- "#d3d3d3"
# # print(colorDF$Color[colorDF$Cluster == "Cluster_NA"])
# # print("original colors")
# # print(colorDF)

# colLen2 <- length(unique(epiCurve_TGCs$tgcName))
# myCols2 <- col_vector[1:colLen2]
# colorDF2 <- data.frame(matrix(ncol = 2, nrow =colLen2))
# colorDF2[,1] <- unique(epiCurve_TGCs$tgcName)
# colorDF2[,2] <- myCols2
# colnames(colorDF2) <- c("TGC", "Color")
# colorDF2$TGC[colorDF2$TGC == "TGC_NA"] <- "#d3d3d3"

# ui <- fluidPage(theme = shinytheme("cosmo"),
#   navbarPage("Cyclo Analysis Portal",
# 	tabPanel("Cyclospora cayetanensis U.S. Map Visualization",

# 	sidebarLayout(
# 		sidebarPanel(style = "position:fixed;width:30.5%;",
# 			tabsetPanel(
# 				tabPanel( 
# 					helpText("Overview"),
# 					strong(h4("Overview")),
# 					p("The visualizations in this document are meant to supplement epidemiological investigations of ", em("Cyclospora cayetanensis."), "The data represented here 
# 						have been collated from various outputs of the ",  em("Cyclospora cayetanensis"), " bioinformatic workflow. This means the data represented here only come from 
# 						specimens, or sequences, that have been submitted to CDC for genotyping. In most rows, users can choose from different options under the", 
# 						strong('Choose Variables to Display'), "tab to display different maps/graphs for the given variable. For example, in Row C, users can choose different Temporal Genetic Clusters (TGCs)
# 						from the dropdown menu, and the map and plot in Row C will be updated to reflect abundance/collection date* information for the selected TGC. The final row of the page 
# 						(Row G) allows users to select any of the plots from Rows A - D and view the plot as a larger image. Please refer to the", strong('Instructions (and tips) for Use'), "section below
# 						for more information on the how to get the most out of this page. Pleae refer to the", strong('Legends'), "tab for more detailed information on each plot."),
# 					p("*All dates represented on the graph are dates of specimen collection provided by State Public Health Labs."),
# 					p("Please take special note of how the axis value change in different plots when different TGCs, Clusters, etc. are chosen. Each cluster has a different size and the plots automatically resize to 
# 						fit the data.")),
# 					tabPanel(
# 						helpText("Instructions"),
# 					strong(h4("Instructions (and tips) for Use")),
# 					p("Plots are organized by row, from Row A to Row F, and the plots in Rows C through F are reactive to user input, while the plots in Rows A and B are static. The options for display on each plot 
# 						can be seen in the various menus on the", strong('Choose Variables to Display'), "tab."),
# 					p("In each plot, users can hover over a portion of the plot to retrieve additional information about that segment of the graph."),
# 					p("On the epidemiological curves, users have two options for highlighting segments of interest. First, users can click on a single color bar in the plot itself, this will enhance all bar segments with this color (Genetic Cluster, TGC, Epi Cluster, or State depending on the plot)
# 						and bar segments from other colors will be greyed out (but still visible). Users can only highlight a single category this way and users can exit this view by double clicking anywhere in the white space on the plot.
# 						 Second, users can double click on the category of interest in the legend. This will hide all other categories from the plot; however, additional categories can be seen when single clicking on that category in the legend.
# 						 All other categories can be brough back by double clicking on the category that the user originally clicked on. Alternatively, users can single click on an element in the legend, and this will make that element disappear 
# 						 from the plot - the element can be brought back by clicking on that element in the legend."),
# 					p("Each plot has a control bar that users will see upon moving the cursor the top right of the plot. The save image option will save the plot exactly how the user sees the plot at the 
# 						moment the button was clicked - the downloaded PNG file will be called 'newplot.png'. The pan, zoom in, zoom out, and reset axes are best used in together. Users can zoom in or out by clicking the respective button, 
# 						then select the pan button to move around the plot, then finally press the home button to zoom back to the original view. The final option in the control bar is an option to control hover behavior; either to
# 						hover on a single element that the cursor is on, or to hover on multiple items near by (which will compare hover text). The default is to hover on a single item.")),
# 				tabPanel(
# 					helpText("Choose Variables to Display"),
# 					h5(strong("Rows A and B: Plots are Static")),
# 					p("No input choice"),
# 					selectInput("myTGC",label = "Row C: Choose a TGC to Display",choices = currentClus_tgc, selected = "TGC_2021_001"),
# 					selectInput("myEpi",label = "Rows D and E: Choose an Epi Cluster to Display",choices = currentClus_epi, selected = "Epi_Cluster_2021 August Butter 1"),
# 					# selectInput("graph", label = "Choose a plot to view in detail:", choices = c("Row A: Collection Date by TGC", "Row A: Collection Date by Epi Cluster", "Row B: TGC Membership by State", "Row B: TGC Membership Epi Curve")),
# 					checkboxGroupInput("myHaps",label = "Row F: Choose TGC Haplotypes to Compare (2 or more)", choices = choicesBox, inline = T),
# 					radioButtons("graph", label = "Row G: Choose a plot to view in detail:",
# 						choices = c("Row A: Epi Curve by Genetic Cluster", "Row B: Epi Curve by TGC", "Row B: Epi Curve by Epi Cluster", "Row C: TGC Membership by State", "Row C: TGC Membership Sunburst Plot",
# 							"Row D: Epi Cluster Membership by State", "Row D: Epi Cluster Epi Curve"))),
# 				tabPanel(
# 					helpText("Legends"),
# 					strong(h3("Legends for Each Row")),
# 					strong(h4("Row A")),
# 					p("Epidemiological curve for current season specimens and sequences sent for", em("Cyclospora cayetanensis."), 
# 						"The color categories in the legend represent genetic clusters where specimens from the current season clustered."),
# 					strong(h4("Row B")),
# 					p("Epidemiological curves for all specimens and sequences sent for", em("Cyclospora cayetanensis."),
# 						"First plot: The color categories in the legend represent temporal genetic clusters (TGCs) where specimens from the current season clustered. 
# 						Second plot: The color categories in the legend represent epidemiological cluster where specimens from the current season clustered."),
# 					strong(h4("Row C")),
# 					p("First plot: Map of U.S.A., colored by TGC abundance. Second plot: Epi linkages within selected TGC"),
# 					strong(h4("Row D")),
# 					p("First plot: Map of U.S.A., colored by Epi Cluster abundance. Second plot: Epi curve for each state with specimens in this Epi Cluster"),
# 					strong(h4("Row E")),
# 					p("Sunburst plot with proportional cluster and TGC membership of specimens in the Epi Cluster"),
# 					strong(h4("Row F")),
# 					p("Haplotype barcode for selected TGCs"),
# 					strong(h4("Row G")),
# 					p("Larger view of selected plot"))
				
# 				# 	helpText("Legends"),
# 				# 	strong(h3("Legends for Each Row")),
# 				# 	strong(h4("Row A")),
# 				# 	p("legend for Row A"),
# 				# 	strong(h4("Row B")),
# 				# 	p("legend for Row B"),
# 				# 	strong(h4("Row C")),
# 				# 	p("legend row C"),
# 				# 	strong(h4("Row D")),
# 				# 	p("legend row D"),
# 				# 	strong(h4("Row E")),
# 				# 	p("Legend Row E")))

			
# 		),
# 			tags$head(
# 			    tags$style(
# 			      HTML(
# 			        ".checkbox-inline { 
# 			                    margin-left: 0px;
# 			                    margin-right: 10px;
# 			          }
# 			         .checkbox-inline+.checkbox-inline {
# 			                    margin-left: 0px;
# 			                    margin-right: 10px;
# 			          }
# 			        "
# 			      )
# 			    ) 
# 			  )
			
# ),
# 		mainPanel( 
# 			fluidRow(h1("Row A"), plotlyOutput("allGenetics", height = "450px", width= "95%")),
# 			fluidRow(h1("Row B"), splitLayout(plotlyOutput("fullTGCs", width = "100%", height = "600px"),plotlyOutput("epiCurve", width = "98%", height = "600px"))),
# 			fluidRow(h1("Row C"), splitLayout(plotlyOutput("plotlyMap", width = "100%", height = "500px"),plotlyOutput("plotlySunburst_TGC", width = "98%", height = "500px"))),
# 			fluidRow(h1("Row D"), splitLayout(plotlyOutput("plotlyMap2", width = "100%", height = "500px"),plotlyOutput("epiCurve_ind", width = "98%", height = "500px"))),
# 			# fluidRow(h1("Row D"), splitLayout(plotlyOutput("plotlySunburst"), plotlyOutput("plotlyBarcode", height = "500px", width = "1200px"))),
# 			fluidRow(h1("Row E", plotlyOutput("plotlySunburst"))),
# 			fluidRow(h1("Row F", plotlyOutput("plotlyBarcode", height = "650px", width = "95%"))),
# 			fluidRow(h1("Row G", plotlyOutput("selected_graph", height = "750px", width = "95%"))),
# 		),

# 		)),
# #   tabPanel("Genetic Strains with Epi History",
# #            sidebarLayout(
# #              sidebarPanel("Test",
# # 			 	tabsetPanel(
# # 					 tabPanel(
# # 						 helpText("Genetic Strains Background"),
# # 						 p("This page displays how current current genetic strains map to epidemiological clusters from previous years. In the first table, each row is a genetic strain and each column is an epidemiological cluster from 2018-2021. The numbers represent the frequency of the strain in the epidemiological cluster from samples in our refernece population. For this reason, the values are not the total number of specimens in the epidemiological cluster, they are only meant to represen the primary strain for each epidemiological cluster"),

# # 						 checkboxGroupInput("myStrainHaps",label = "Row B: Choose Strain Haplotypes to Compare (2 or more)", choices = strainChoices, inline = T),
# # 						 downloadButton("downloadStrainEpi", "Download Strain + Epi History")

# # 					 )
# # 				 )),
# #            mainPanel(width = 8,
# # 		   h1("Strain Names Effective: ", strainEffectiveDate),
# #         #    fluidRow(h2("Row A"),DT::dataTableOutput("strainHistory")),
# # 		   fluidRow(h2("Row B", plotlyOutput("strainBarcode", height = "650px", width = "95%"))),
# #            )),
# # 	),
#   tabPanel("Genetic Strain Tables",
#   	sidebarLayout(
# 		  sidebarPanel("Test",
# 		  	tabsetPanel(
# 				  tabPanel(
# 					  helpText("Background"),
# 					#   p("Descript"),
# 					  p("This page displays how current current genetic strains map to epidemiological clusters from previous years. In the first table, each row is a genetic strain and each column is an epidemiological cluster from 2018-2021. The numbers represent the frequency of the strain in the epidemiological cluster from samples in our refernece population. For this reason, the values are not the total number of specimens in the epidemiological cluster, they are only meant to represen the primary strain for each epidemiological cluster")

# 					#   checkboxGroupInput("myStrainHaps",label = "Row B: Choose Strain Haplotypes to Compare (2 or more)", choices = strainChoices, inline = T)
# 					  ),
# 				  tabPanel(
# 					  helpText("Download CSV Files"),
# 					  p('Click on the buttons to download a csv file of each table'),
# 					  downloadButton("downloadStrainEpi", "Download Strain + Epi History"),
# 					  br(),
# 					  br(),
# 					  downloadButton("downloadStrainNames", "Download Strain Names"),
# 					  br(),
# 					  br(),
# 					  downloadButton("downloadMtKey", "Download Mitochondrial Key File"),
# 					  br(),
# 					  br(),
# 					  downloadButton("downloadNuKey", "Download Nuclear Key File")
# 				  ),
# 				  tabPanel(
# 					  helpText("Strain Barcode Choices"),
# 					  p("Choose Strain Haplotypes to Compare"),
# 					  checkboxGroupInput("myStrainHaps",label = "Select 2 or more", choices = strainChoices, inline = T)
# 					  )
# 			  )),
# 			mainPanel(width = 8,
# 				h1("Strain Names Effective: ", strainEffectiveDate),
# 				tabsetPanel(id = 'dataset',
# 				tabPanel("Strains + Historical Epi Links", DT::dataTableOutput("strainHistory")),
# 				tabPanel("Strain Names",DT::dataTableOutput("strainNames")),
# 				tabPanel("Mitochondrial Key",DT::dataTableOutput("mtKey")),
# 				tabPanel("Nuclear Key", DT::dataTableOutput("nuKey"))
				
				
# 				),
# 			fluidRow(h2("Row B", plotlyOutput("strainBarcode", height = "650px", width = "95%")))
# 	  )),
#   )
# #   tabPanel("Nuclear Key",
# #   	sidebarLayout(
# # 		  sidebarPanel("Test",
# # 		  	tabsetPanel(
# # 				  tabPanel(
# # 					  helpText("test1"),
# # 					  p('Desciption'),
# # 					  downloadButton("downloadNuKey", "Download Nuclear Key File")
# # 				  )
# # 			  )),
# # 			mainPanel(width = 8,
# # 				fluidRow(DT::dataTableOutput("nuKey"))
# # 	  )),
# #   )
# )
# )


# server <- function(input, output){
#   output$strainHistory <- DT::renderDataTable({
# 	epiStrains_countReorder
#    })
# 	output$mtKey <- DT::renderDataTable({
# 	datatable(mtKeyFile, rownames = FALSE)

#    })
#    	output$nuKey <- DT::renderDataTable({
# 	datatable(nuKeyFile, rownames = FALSE)
#    })
#     output$strainNames <- renderDataTable({
# 	datatable(latestStrainNames, rownames = FALSE)
#    })
#    output$downloadNuKey <- downloadHandler(
# 	   filename = function() {
# 		   paste0("cyclospora_nuKey", strainEffectiveDate, ".csv", sep = "")
# 		    },
# 	   content = function(file){
# 		   write.table(nuKeyFile, file, row.names = F, sep =",", quote = F)
# 	   } 
#    )
#     output$downloadMtKey <- downloadHandler(
# 	   filename = function() {
# 		   paste0("cyclospora_mtKey_effective_", strainEffectiveDate,".csv", sep = "")
# 		    },
# 	   content = function(file){
# 		   write.table(mtKeyFile, file, row.names = F, sep = ",", quote = F)
# 	   } 
#    )
#     output$downloadStrainNames <- downloadHandler(
# 	   filename = function() {
# 		   paste0("cyclospora_strain_names_effective_", strainEffectiveDate, ".csv", sep = "")
# 		    },
# 	   content = function(file){
# 		   write.table(latestStrainNames, file, row.names = F, sep =",", quote = F)
# 	   } 
#    )
#       output$downloadStrainEpi <- downloadHandler(
# 	   filename = function() {
# 		   paste0("cyclospora_strain_epi_history_effective_", strainEffectiveDate, ".csv", sep = "")
# 		    },
# 	   content = function(file){
# 		   write.table(epiStrains_countReorder, file, row.names = F, sep = ",", quote = F)
# 	   } 
#    )
# 	plot0 <- reactive({
# 		aVar <- input$currentTGC
# 		allGeneticPlot(epiCurve_data, aVar)
# 		})
# 	plot1 <- reactive({
# 		aVar <- input$currentTGC
# 		allTGCsPlot(epiCurve_TGCs, aVar)
# 		})
# 	plot2 <- reactive({
# 		aVar <- input$myEpi
# 		allEpiPlot(epiCurve_epi, aVar)
# 		})
# 	plot3 <- reactive({
# 		newMap <- input$myTGC
# 		plotlyFun1(newTable_tgc, newMap, "Blues", newTable)
# 		})
# 	plot4 <- reactive({
# 	  newSun <- input$myTGC
# 	  sunburstPlot_TGC(fullMerge, newSun)
# 		})
# 	plot5 <- reactive({
# 		newMap2 <- input$myEpi
# 		plotlyFun2(newTable_epi, newMap2, "Purples", fullMerge)
# 		})
# 	plot6 <- reactive({
# 		aVar <- input$myEpi
# 		singleEpiPlot(epiCurve_epi, aVar)
# 		})
# 	graphInput <- reactive({
# 		switch(input$graph, 
# 			"Row A: Epi Curve by Genetic Cluster" = plot0(),
# 			"Row B: Epi Curve by TGC" = plot1(),
# 			"Row B: Epi Curve by Epi Cluster" = plot2(),
# 			"Row C: TGC Membership by State" = plot3(),
# 			"Row C: TGC Membership Sunburst Plot" = plot4(),
# 			"Row D: Epi Cluster Membership by State" = plot5(),
# 			"Row D: Epi Cluster Epi Curve" = plot6()
# 			)
# 		})


# 	output$selected_graph <- renderPlotly({
# 		graphInput()
# 		})

# 	output$fullTGCs <- renderPlotly({
# 		aVar <- input$currentTGC
# 		allTGCsPlot(epiCurve_TGCs, aVar)
# 		})
# 	output$individualTGCs <- renderPlotly({
# 		aVar <- input$myTGC
# 		singleTGCPlot(epiCurve_TGCs, aVar)
# 		})
# 	output$epiCurve<- renderPlotly({
# 		aVar <- input$myEpi
# 		allEpiPlot(epiCurve_epi, aVar)
# 		})
# 	output$epiCurve_ind <- renderPlotly({
# 		aVar <- input$myEpi
# 		singleEpiPlot(epiCurve_epi, aVar)
# 		})

# 	output$allGenetics <- renderPlotly({
# 		aVar <- input$myEpi
# 		allGeneticPlot(epiCurve_data, aVar)
# 		})

# 	output$plotlyMap <- renderPlotly({
# 		newMap <- input$myTGC
# 		plotlyFun1(newTable_tgc, newMap, "Blues", newTable)
# 		})
# 	output$plotlyMap2 <- renderPlotly({
# 		newMap2 <- input$myEpi
# 		plotlyFun2(newTable_epi, newMap2, "Purples", fullMerge)
# 		})
# 	output$plotlySunburst <- renderPlotly({
# 		newPie <- input$myEpi
# 		sunburstPlot(fullMerge, newPie)

# 		})
# 	output$plotlySunburst_TGC <- renderPlotly({
# 	  newSun <- input$myTGC
# 	  sunburstPlot_TGC(fullMerge, newSun)
	  
# 	})
# 	output$plotlyBarcode <- renderPlotly({
# 		aVar <- input$myHaps
# 		validate(
# 			need(length(aVar) >= 2, "Please Select at least Two Clusters"))
		
# 		hapBarcodePlot(barcodeMat, aVar)
# 		})
# 	output$strainBarcode <- renderPlotly({
# 		bVar <- input$myStrainHaps
# 		validate(
# 			need(length(bVar) >= 2, "Please Select at least Two Clusters"))
		
# 		strainBarcodePlot(strainBarcodeMat, bVar)
# 		})

# }
# shinyApp(ui = ui, server = server)