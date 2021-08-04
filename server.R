

library("readxl")
library(shiny)
library(reshape2)
library(xlsx)
library(dplyr)
library(tibble)
library(stringr)
library(openxlsx)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    
    getData <- reactive({
        if(is.null(input$files)){
            return(NULL)
        } else {
            
            #result_header <- c("#", "Barcodes", "Lab ID", "E", "N", "RdRp", "GAPDH", "Interpretation", "Assay", "Date", "Used Extraction Instrument", "Used PCR Instrument")
            
            #result_df <- setNames(data.frame(matrix(ncol=12, nrow=0)), result_header)
            
            nfiles <- nrow(input$files)
            pcr_data = list()
            
            # Declare an empty data frame for the result
            res = data.frame()
            qc = data.frame()
            ntcExtQC <- c()
            ntcPCRQC <- c()
            pcQC <- c()
            
            for (i in 1:nfiles) {
                # Find the Starting Index of the data using the Row Index of 'Sample Name'
                f = read_excel(input$files$datapath[i],col_names = FALSE)
                #f = read_excel(input$files[[i, 'datapath']],col_names = FALSE)
                skipIndex <- which(f == 'Sample Name')
                
                #pcr_data <- read_excel(input$files$datapath[i], skip=skipIndex-1)
                pcr_data[[i]] <- read_excel(input$files[[i, 'datapath']], skip=skipIndex-1)
                
                # Replace NA with 0 for CT mean 
                ct_mean <- names(pcr_data[[i]])[4]
                pcr_data[[i]][[ct_mean]] <- ifelse(is.na(pcr_data[[i]][[ct_mean]]), 0, pcr_data[[i]][[ct_mean]])
                
                # Add E gene CT mean to result_df
                eData <- pcr_data[[i]] %>% filter(`Target Name` == 'E', ignore.case = TRUE) %>% select(-c(3))
                eData$`Ct Mean` <- round(as.numeric(eData$`Ct Mean`),1)
                colnames(eData)[1] <- "Lab ID"
                
                # Add N gene CT mean to result_df
                nData <- pcr_data[[i]] %>% filter(`Target Name` == 'N',ignore.case = TRUE) %>% select(-c(3))
                nData$`Ct Mean` <- round(as.numeric(nData$`Ct Mean`),1)
                
                # Add RdRp gene CT mean to result_df
                rdData <- pcr_data[[i]] %>% filter(`Target Name` == 'RdRp',ignore.case = TRUE) %>% select(-c(3))
                rdData$`Ct Mean` <- round(as.numeric(rdData$`Ct Mean`),1)
                
                # Add GAPDH CT value to result_df
                gapData <- pcr_data[[i]] %>% filter(`Target Name` == 'GAPDH',ignore.case = TRUE)
                gapData$`Ct Mean` <- round(as.numeric(gapData$`Ct Mean`),1)
                gapData$`Ct Mean`[gapData$`Ct Mean` == 0] <- 45.0
                #gapData <- gapData %>% select(-c('CT'))
                #gapData$`CT` <- round(as.numeric(gapData$`CT`),1)
                #gapData$`CT`[is.na(gapData$`CT`)] <- 'Undetermined'
                
                
                # Paste into Dataframe (E, N, RdRp, GAPDH in an order)
                combine <- data.frame('Lab Id'=eData$`Lab ID`, 'E'=eData$`Ct Mean`, 'N'=nData$`Ct Mean`, 'RdRp'=rdData$`Ct Mean`, 'GAPDH'=gapData$`Ct Mean`)
                names(combine)[1] <- 'Lab Id'
                
                
                
                ### Generate QC Data ### 
                # Extract Lab ID from the filename 
                # Remove Date format in front
                labID <- sub(".xls$", "", basename(input$files$name[i]))

                labID <- labID %>%
                            str_replace(pattern="\\d{4}-\\d{1,2}-\\d{1,2}", "") %>%
                            str_trim() %>%
                            str_replace_all(" ", repl="-") %>%
                            str_replace_all("-", repl="_")
                
                # Get Date
                testDate <- sub(".xls$", "", basename(input$files$name[i]))
                testDate <- testDate %>% str_extract(pattern="\\d{4}-\\d{1,2}-\\d{1,2}")
                
                # Extract NTC-Extraction & NTC-PCR & Positive Control 
                #ntcExtract <- combine %>% filter(`Lab Id` == 'NTC-Extraction',ignore.case = TRUE)
                ntcExtract <- combine %>% filter(grepl('NTC-Extraction', `Lab Id`,ignore.case = TRUE))
                ntcExtract$GAPDH[ntcExtract$GAPDH == 'Undetermined'] <- 0.0
                ntcExtract$E[ntcExtract$E == 0] <- 45.0
                ntcExtract$N[ntcExtract$N == 0] <- 45.0
                ntcExtract$RdRp[ntcExtract$RdRp == 0] <- 45.0
                ntcExtract$GAPDH[ntcExtract$GAPDH == 0] <- 45.0
                ntcExtract$GAPDH <- as.numeric(ntcExtract$GAPDH)
                #combine <- combine %>% filter(`Lab Id` != 'NTC-Extraction',ignore.case = TRUE)
                combine <- combine %>% filter(!grepl('NTC-Extraction', `Lab Id`,ignore.case = TRUE))
        
                
                #ntcPCR <- combine %>% filter(`Lab Id` == 'NTC-PCR',ignore.case = TRUE)
                ntcPCR <- combine %>% filter(grepl('NTC-PCR', `Lab Id`,ignore.case = TRUE))
                ntcPCR$GAPDH[ntcPCR$GAPDH == 'Undetermined'] <- 0.0
                ntcPCR$E[ntcPCR$E == 0] <- 45.0
                ntcPCR$N[ntcPCR$N == 0] <- 45.0
                ntcPCR$RdRp[ntcPCR$RdRp == 0] <- 45.0
                ntcPCR$GAPDH[ntcPCR$GAPDH == 0] <- 45.0
                ntcPCR$GAPDH <- as.numeric(ntcPCR$GAPDH)
                #combine <- combine %>% filter(`Lab Id` != 'NTC-PCR',ignore.case = TRUE)
                combine <- combine %>% filter(!grepl('NTC-PCR', `Lab Id`,ignore.case = TRUE))
                
                pc <- combine %>% filter(grepl('Positive Control', `Lab Id`,ignore.case = TRUE),ignore.case = TRUE)
                combine <- combine %>% filter(!grepl('Positive Control', `Lab Id`,ignore.case = TRUE),ignore.case = TRUE)

                # Quality PASS/FAIL Test
                if(ntcExtract$E >= 38.0 & ntcExtract$N >= 38.0 & ntcExtract$RdRp >= 38.0 & ntcExtract$GAPDH >= 38.0) {
                    ntcExtQC <- c(ntcExtQC, "Pass")
                } else {
                    ntcExtQC <- c(ntcExtQC, "Fail")
                }
                
                if(ntcPCR$E >= 38.0 & ntcPCR$N >= 38.0 & ntcPCR$RdRp >= 38.0 & ntcPCR$GAPDH >= 38.0) {
                    ntcPCRQC <- c(ntcPCRQC, "Pass")
                } else {
                    ntcPCRQC <- c(ntcPCRQC, "Fail")
                }
                
                if(pc$E < 38.0 & pc$N < 38.0 & pc$RdRp < 38.0 & pc$GAPDH < 38.0) {
                    pcQC <- c(pcQC, "Pass")
                } else {
                    pcQC <- c(pcQC, "Fail")
                }
                
                
                
                qcData <- data.frame('Plate#'=NA, 'Date'=testDate, 'Plate lot #'=NA, 'PCR Instrument #'=NA, 'Lab ID'=labID, 'Non-Template Control-Extraction Lot#'=NA, 
                    'NTC-Ext_E'=ntcExtract$E, 'NTC-Ext_N'=ntcExtract$N, 'NTC-Ext_RdRp'=ntcExtract$RdRp,'NTC-Ext_GAPDH'=ntcExtract$GAPDH, 'NTC Ext Quality'=NA,
                    'Non-Template Control-RT-qPCR Lot#'=NA, 'NTC-PCR_E'=ntcPCR$E, 'NTC-PCR_N'=ntcPCR$N, 'NTC-PCR_RdRp'=ntcPCR$RdRp,'NTC-PCR_GAPDH'=ntcPCR$GAPDH,'NTC PCR Quality'=NA,
                    'Positive Control Lot#'=NA, 'PC_E'=pc$E, 'PC_N'=pc$N, 'PC_RdRp'=pc$RdRp,'PC_GAPDH'=pc$GAPDH, 'PC Quality'=NA)
                
                names(qcData) <- gsub(".", " ", names(qcData), fixed = TRUE)
                
                # Finally, do the rBIND 
                res <- rbind(res, combine)
                qc <- rbind(qc, qcData)

            }
            qc$`NTC Ext Quality` <- ntcExtQC
            qc$`NTC PCR Quality` <- ntcPCRQC
            qc$`PC Quality` <- pcQC
            
            res <- cbind('#'=NA, 'Barcodes'=NA,res,"Interpretation"=NA, "Assay"='TrioDx', "Date"=testDate, "Used Extraction Instrument"='KingFisher Flex', "Used PCR Instrument"='QuantStudio6' )
            
            return(list(result = res, qc = qc))
            
            # GAPDH CT value -> to Numeric

        }
    })
    
    
    #output$default <- renderText({     getData() })

    output$contents <- renderTable( 
        getData()$result 
    )
    
    output$qcTable <- renderTable( 
        getData()$qc
    )
    
    
    output$download_merged <- downloadHandler(
        filename=function(){
            paste("Result_", Sys.Date(), ".xlsx", sep="")
        },
        content=function(file){
            
            retest <- data.frame('Index'=NA,'Barcode'=NA, 'Lab ID'=NA, 'E'=NA, 'N'=NA, 'RdRp'=NA, 'GAPDH'=NA, 'Interpretation'=NA, 'RNA Plate#'=NA, 'RNA Well#'=NA, 'Re-PCR Plate Well#'=NA, 'E'=NA, 'N'=NA, 'RdRp'=NA, 'GAPDH'=NA, 'Interpretation'=NA)
            
            wb <- createWorkbook()
            addWorksheet(wb,sheetName='Result')
            addWorksheet(wb,sheetName='DAILY QC')
            addWorksheet(wb,sheetName='RETEST')

            redStyle <- createStyle(fontColour="#9C0006", bgFill="#FFC7CE")
            yellowStyle <- createStyle(fontColour="#9C5700", bgFill="#FFEB9C")
            mergedStyle <- createStyle(halign="center", fontSize = 18, textDecoration="bold")
            mergedStyle2 <- createStyle(halign="center", textDecoration="bold")
            lightBorder <- createStyle(border='TopBottomLeftRight', borderStyle='thin')
            qcHeaderStyle <- createStyle(halign = "center",wrapText = TRUE,valign = "center", border='TopBottomLeftRight', borderStyle='thin', textDecoration="bold")
            #headerBorder <- createStyle(border=c('top','left','right'), borderStyle='medium')
            #headerBorder2 <- createStyle(border=c('bottom','left','right'), borderStyle='medium')
            
            #rightBorder <- createStyle(border='right', borderStyle = 'medium')
            
            writeDataTable(wb,sheet='Result', x=getData()$result, keepNA=FALSE, rowNames=FALSE)
            
            writeData(wb, sheet='DAILY QC', 'Daily Quality Control', startCol = 1, startRow = 1, colNames = FALSE, rowNames = TRUE)
            writeData(wb, sheet='DAILY QC','Non Template Control -Extraction', startCol = 6, startRow = 2, colNames = FALSE, rowNames = TRUE)
            writeData(wb, sheet='DAILY QC','Non Template Control - RT-qPCR', startCol = 12, startRow = 2, colNames = FALSE, rowNames = TRUE)
            writeData(wb, sheet='DAILY QC','Positive Control', startCol = 18, startRow = 2, colNames = FALSE, rowNames = TRUE)
            
            writeData(wb,sheet='DAILY QC', x=getData()$qc, headerStyle=qcHeaderStyle, borders='all', borderStyle = 'thin',  startRow=3, keepNA=FALSE, rowNames=FALSE)
            writeData(wb, sheet='RETEST', x=retest, keepNA=FALSE, rowNames=FALSE)
            
            mergeCells(wb, sheet='DAILY QC', cols=1:length(getData()$qc), rows=1)
            mergeCells(wb, sheet='DAILY QC', cols=1:5, rows=2)
            mergeCells(wb, sheet='DAILY QC', cols=6:11, rows=2)
            mergeCells(wb, sheet='DAILY QC', cols=12:17, rows=2)
            mergeCells(wb, sheet='DAILY QC', cols=18:length(getData()$qc), rows=2)
            addStyle(wb, sheet='DAILY QC', mergedStyle, rows=1, cols=1)
            addStyle(wb, sheet='DAILY QC', mergedStyle2, rows=2, cols=6)
            addStyle(wb, sheet='DAILY QC', mergedStyle2, rows=2, cols=12)
            addStyle(wb, sheet='DAILY QC', mergedStyle2, rows=2, cols=18)
            
            #addStyle(wb, sheet='DAILY QC', rightBorder, rows=4:(nrow(getData()$qc)+3), cols=5)
            
            setColWidths(wb, sheet='DAILY QC', ignoreMergedCells = TRUE,cols = 5, widths = "auto")

            conditionalFormatting(wb, sheet='Result', style=yellowStyle, cols=4:6, rows=1:(nrow(getData()$result)+1), rule="AND(D1<38, D1>0)" )
            conditionalFormatting(wb, sheet='Result', style=redStyle, cols=7, rows=1:(nrow(getData()$result)+1), rule="AND(G1>=38, G1<=45)" )

            saveWorkbook(wb, file=file, overwrite=TRUE)

            #write.xlsx(getData()$result,file, sheetName='Result', showNA = FALSE, rowNames=FALSE, append = FALSE)
            #write.xlsx(getData()$qc,file, sheetName='DAILY QC', showNA = FALSE, rowNames=FALSE,append = TRUE)

        }
    )

})
