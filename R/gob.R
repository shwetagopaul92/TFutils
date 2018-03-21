
#' gadget to help sort through tags naming TFs
#' @importFrom shiny selectInput dataTableOutput reactive renderDataTable observeEvent
#' @importFrom shiny runGadget stopApp
#' @import miniUI
#' @importFrom GSEABase geneIdType SymbolIdentifier "geneIdType<-" setName geneIds
#' @param gscoll a GSEABase GeneSetCollection
#' @param initTF character(1) initial TF string for app
#' @rawNamespace importClassesFrom(GenomicRanges, GRanges)
#' @note Will use gwascat_hg19 to look for 'MAPPED_GENE' field entries matching targets.
#' @export
TFtargs = function(gscoll=TFutils::tftColl, initTF="VDR") {
  ui <- miniPage(gadgetTitleBar("Search for a TF; its targets will be checked for mapped status in GWAS catalog"), 
                 miniContentPanel(
                    selectInput("tfsel", "TF:", names(gscoll),
                         selected = initTF),
                    dataTableOutput("tab"))
                ) # end page
  server <- function(input, output, session) {
    getTab = reactive({
      grabTab(input$tfsel, TFutils::tftColl, org.Hs.eg.db::org.Hs.eg.db, TFutils::gwascat_hg19)
      })
    output$tab <- renderDataTable({
      getTab()
    })
    observeEvent(input$done, {
      stopApp(getTab())
    })
  }
  runGadget(ui, server)
}
