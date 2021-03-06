

library("readxl")
library(shiny)
library("shinydashboard")
library("dplyr")
library(reshape2)
library(xlsx)
library(tibble)
library(stringr)

# Define UI for application that draws a histogram
ui <- dashboardPage(
    dashboardHeader(title = 'TrioDx Test'),
    
    dashboardSidebar(
        sidebarMenu(
            menuItem("Merge Files",
                     tabName="cnv_merger",
                     icon=icon("cloud-upload-alt", lib="font-awesome"))
        )
    ),
    
    dashboardBody(
        tabItems(
            tabItem(tabName = "cnv_merger",
                    fluidRow(
                        box(
                            title="File Upload",
                            fileInput("files", "Choose Files to be Merged", width='80%', multiple=TRUE, accept=c(".xlsx", 'xls')),
                            downloadButton("download_merged","Download Result"),
                            verbatimTextOutput("default"),
                            #verbatimTextOutput("debug"),
                            #tableOutput('contents'),
                            align='center',
                            width=6
                        )),
                    fluidRow(
                        
                        box(
                            title="QC Table Lookup",
                            tableOutput('qcTable'),
                            align='center',
                            width=12
                        )
                        
                    ),
                    fluidRow(
                        box(
                            title="Merged Data Table Lookup",
                            tableOutput('contents'),
                            align="center",
                            width=12
                        )
                    )
                    )
        )
        
    )
)

