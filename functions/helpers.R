resetInput_to_initialState <-function(output){
    shinyjs::reset("proteinOrder")
    shinyjs::reset("line_dot_size")
    shinyjs::reset("wordsize")
    shinyjs::reset("supportLimit")
    shinyjs::reset("kmerlength")
    shinyjs::reset("filetype")
    shinyjs::reset("host")
    shinyjs::reset("proteinNames")
    shinyjs::reset("hostname")
    shinyjs::reset("proteinNames_secondHost")
    shinyjs::reset("hostname_secondHost")
    shinyjs::reset("MSAfile")
    shinyjs::reset("MSAfile_secondHost")
    shinyjs::reset("Metafile")
    shinyjs::reset("Metafasta")
    shinyjs::reset("Meta")
    shinyjs::reset("inmetafilename")
    shinyjs::reset("inmetafasta")
    shinyjs::disable(id="downloadDiMA")

    #clear output
    output$alert <- renderUI({})
    output$alertSample <- renderUI({})
    #output$protein_selection <- renderUI({})
    output$plot_worldmap<- renderPlot({})
    output$countrytable<- renderDataTable({})
    output$plot_time<- renderPlot({})
    output$timetable<- renderDataTable({})
    output$metademoSee<- renderDataTable({})

    output$plot1<- renderPlot({})
    output$plotEntropy<- renderPlot({})
    output$plot2<- renderPlot({})
    output$plot3<- renderPlot({})
    output$plot4<- renderPlot({})
    output$plot7<- renderPlot({})
    output$entropyTable<-renderDataTable({})
    output$plot7_seqs<-renderDataTable({})
}

resetMetaDataInput <- function(output){
    shinyjs::reset("Metafile")
    shinyjs::reset("Metafasta")
    shinyjs::reset("Meta")
    shinyjs::reset("inmetafilename")
    shinyjs::reset("inmetafasta")
    shinyjs::disable(id="downloadDiMA")

    #clear output
    output$alert <- renderUI({})
    output$alertSample <- renderUI({})
    output$plot_worldmap<- renderPlot({})
    output$countrytable<- renderDataTable({})
    output$plot_time<- renderPlot({})
    output$timetable<- renderDataTable({})
    output$metademoSee<- renderDataTable({})

}

#download sample data
downloadSampleData<-function(output){
    output$downloadSampleData <- downloadHandler(
        filename <- function() {
            paste("vDiveR_sample_input_", Sys.Date(), ".zip", sep = "")
        },

        content <- function(file) {
            file.copy("www/vDiveR_sample_input.zip", file)
        },
        contentType = "application/zip"
    )
}

#plot generators
generate_worldmap <-function(input, output, plot_worldmap){
    output$plot_worldmap <- renderPlot({
        plot_worldmap()
    })

    output$plot_worldmap_download <- downloadHandler(
        filename = function() { paste("plot_world_map", '.jpg', sep='') },
        content = function(file) {
            ggsave(file, plot = plot_worldmap(), width=input$width_wm, height=input$height_wm,unit="in", device = "jpg", dpi=input$dpi_wm)
        }
    )

}

generate_timeplot <-function(input, output, plot_time){
    output$plot_time <- renderPlot({
        plot_time()
    })

    output$plot_time_download <- downloadHandler(
        filename = function() { paste("plot_time", '.jpg', sep='') },
        content = function(file) {
            ggsave(file, plot = plot_time(), width=input$width_tm, height=input$height_tm,unit="in", device = "jpg", dpi=input$dpi_tm)
        }
    )

}

generate_plot1<-function(input, output, plot1){
    output$plot1 <- renderPlot({
        plot1()
    })

    output$plot1_download <- downloadHandler(
        filename = function() { paste("plot_entropy_incidence", '.jpg', sep='') },
        content = function(file) {
            ggsave(file, plot = plot1(), width=input$width, height=input$height,unit="in", device = "jpg", dpi=input$dpi)
        }
    )

}

generate_plotEntropy<-function(input, output, plotEntropy){
    output$plotEntropy <- renderPlot({
        plotEntropy()
    })

    output$plotEntropy_download <- downloadHandler(
        filename = function() { paste("plot_entropy", '.jpg', sep='') },
        content = function(file) {
            ggsave(file, plot = plotEntropy(), width=input$width, height=input$height,unit="in", device = "jpg", dpi=input$dpi)
        }
    )

}

generate_plot2<-function(input, output, plot2){
    output$plot2<- renderPlot({
        plot2()
    })

    output$plot2_download <- downloadHandler(
        filename = function() { paste("plot_relationship_entropy_total_variants", '.jpg', sep='') },
        content = function(file) {
            ggsave(file, plot = plot2(), width=input$width2, height=input$height2,unit="in", device = "jpg", dpi=input$dpi2)
        }
    )

    output$info_plot2 <- renderText({
        HTML(paste0("Total variants (%) = ", input$plot2_click$x, "<br>",em("k"),"-mer entropy (bits) = ",input$plot2_click$y))
    })
}

generate_plot3<-function(input, output, plot3){
    output$plot3<- renderPlot({
        plot3()
    })

    output$plot3_download <- downloadHandler(
        filename = function() { paste("plot_dynamics_diversity_motifs_proteome", '.jpg', sep='') },
        content = function(file) {
            ggsave(file, plot = plot3(), width=input$width3, height=input$height3,unit="in", device = "jpg", dpi=input$dpi3,bg='white')
        }
    )
}

generate_plot4<-function(input, output, plot4){
    output$plot4<- renderPlot({
        plot4()
    })

    output$plot4_download <- downloadHandler(
        filename = function() { paste("plot_dynamics_diversity_motifs_proteins", '.jpg', sep='') },
        content = function(file) {
            ggsave(file, plot = plot4(), width=input$width4, height=input$height4,unit="in", device = "jpg", dpi=input$dpi4, bg='white')
        }
    )
}

generate_plot7<-function(input, output, plot7){
    output$plot7<- renderPlot({
        plot7()
    })

    output$plot7_download <- downloadHandler(
        filename = function() { paste("plot_conservation_levels_protein", '.jpg', sep='') },
        content = function(file) {
            ggsave(file, plot = plot7(),  width=input$width7, height=input$height7, unit="in", device = "jpg", dpi=input$dpi7)
        })
}


generate_entropyTable<-function(data, output, proteinName){
    #get position of min entropy, min, max of entropy and total variant
    entropyTable <- data %>%
        dplyr::group_by(proteinName) %>%
        dplyr::summarise(
            Position = gsub("((?:\\d+,){2}\\d+),", "\\1,\n", paste0(position[which(entropy == min(entropy))], collapse = ",")),
            minEntropy = format(round(min(entropy),digits=2),nsmall=2),
            maxEntropy = format(round(max(entropy), digits = 2),nsmall=2),
            minTotalVariants = format(round(min(totalVariants.incidence),digits=2),nsmall=2),
            maxTotalVariants = format(round(max(totalVariants.incidence), digits = 2),nsmall=2)
        )
    #rename table df
    names(entropyTable)<- c("Protein Name","Position (Minimum Entropy)","Minimum Entropy","Maximum Entropy","Minimum Total Variants (%)","Maximum Total Variants (%)")

    output$entropyTable <- renderDataTable(
        entropyTable,
        width = "100%",
        options = list(scrollX=TRUE, scrollCollapse=TRUE)
    )
}

generate_CCS_HCS_table<-function(input, output, data){
    # for now not splitted by hosts
    output$plot7_seqs <- renderDataTable({
        vDiveR::concat_conserved_kmer(data=data.frame(data), kmer=input$kmerlength,
                              threshold_pct = as.numeric(input$conserv_percent),
                              conservation_level=input$conserv_lvl)[[input$table_type]]
    })

    output$conservSeq_download <- downloadHandler(
        filename =  function() {paste0(input$conserv_lvl, ".", input$table_type)},
        content = function(fname) {
            df <- vDiveR::concat_conserved_kmer(data=data.frame(data), kmer=input$kmerlength,
                                        threshold_pct = as.numeric(input$conserv_percent),
                                        conservation_level=input$conserv_lvl)[[input$table_type]]
            write.table(df, file = fname, col.names = ifelse(input$table_type == "csv", TRUE, FALSE), sep = ",", row.names = FALSE, quote = FALSE)
        }
    )
}

detect_date_format <- function(dates) {
    formats <- c("%d/%m/%Y", "%Y-%m-%d", "%m/%d/%Y", "%Y/%m/%d", "%d-%b-%Y", "%b %d, %Y", "%Y%m%d", "%d-%m-%y", "%m-%d-%Y", "%B %d, %Y") # Add more as needed
    
    for (fmt in formats) {
        parsed <- as.Date(dates, format = fmt)
        if (all(!is.na(parsed))) {
            return(fmt)
        }
    }
    return("Unknown format")
}

get_plot_time_tick_interval <- function(date_column) {
    # Get the range of dates
    date_range <- range(date_column, na.rm = TRUE)
    start_date <- date_range[1]
    end_date <- date_range[2]
    # Calculate the number of months between the two dates
    months_diff <- interval(start_date, end_date) %/% months(1)
    # Calculate the tick interval for 10 ticks
    tick_interval <- ceiling(months_diff / 10)
    return(tick_interval)
}