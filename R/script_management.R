#' Install All Scripts from benchmarkQFeatures
#'
#' This function installs all R scripts located in the `inst/scripts` directory
#' of the benchmarkQFeatures package to a specified destination folder.
#' If no destination folder is provided, it defaults to the current working directory.
#' All installed scripts will be made executable.
#'
#' @param destination_folder A character string specifying the folder where the
#' scripts should be installed. Defaults to the current working directory.
#'
#' @examples
#' # Install all scripts to the current working directory
#' install_all_scripts()
#'
#' @export
installScripts <- function(destinationFolder = getwd()) {
    scriptsDir <- system.file("scripts", package = "benchmarkQFeatures")

    if (scriptsDir == "") {
        stop("Scripts directory not found in the package")
    }

    destinationFolder <- path.expand(destinationFolder)
    if (!dir.exists(destinationFolder)) {
        dir.create(destinationFolder, recursive = TRUE)
        message("Created directory: ", destinationFolder)
    }

    scriptFiles <- list.files(scriptsDir, pattern = "\\.R$",
                               full.names = TRUE)

    if (length(scriptFiles) == 0) {
        stop("No R scripts found in the package scripts directory.")
    }

    for (script in scriptFiles) {
        scriptName <- gsub(".R$", "", basename(script))
        dest <- file.path(destinationFolder, scriptName)
        file.copy(script, dest, overwrite = TRUE)
        Sys.chmod(dest, mode = "0755")
        message("Installed and made executable: ", dest)
    }

    message("All scripts have been installed to: ", destinationFolder)
}
