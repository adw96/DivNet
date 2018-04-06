#' savePDF
#'
#' This function saves a pdf image to a given filepath, with LaTeX fonts
#'
#' @param filepath computer specific filepath for saving, defaults to current working directory
#' @param img command used to generate a plot
#' @param latex whether to use Computer Modern font, defaults to true. Note if there are issues may need font_import(), font_install('fontcm'), loadfonts()
#'
#' @export
savePDF <- function(filepath, img, latex = TRUE) {
    if (latex == TRUE) {
        pdf(filepath, family = "CM Roman")
        img
        dev.off()
        embed_fonts(filepath, outfile = filepath)
    }
    if (latex == FALSE) {
        pdf(filepath)
        img
        dev.off()
    }
}
