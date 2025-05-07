
# define BIO - function to check if BIOUNIT date is file
BIOextract <- function(x){
    tryCatch(
        expr = {
            BIO<-suppressWarnings(biounit(x))
            message("Biounit successfully extracted.")
        },
        error = function(e){
            message('NO Bio unit data in file!')
            print(e)
        }
    )
    return(BIO)
}
#

