# function to import Rda file but with a custom name
import_rda <- function(file, object_name = NA, rm_rda = TRUE) {
  # load the Rda file in the function's environment (not the global one)
  load(file)
  
  # extract the object's name from the file name
  if(is.na(object_name)) {
    object_name <- stringr::str_replace(string = file,
                                        pattern = 'data/tracking-data/',
                                        replacement = '')
    object_name <- stringr::str_replace(object_name, '.Rda', replacement = '')
  }
  
  # return the imported object
  return(get(object_name))
}

# test the function
if(FALSE) {
  head(import_rda('data/tracking-data/Oreamnos_americanus.Rda'))
  head(import_rda('data/tracking-data//Rangifer_tarandus.Rda',
                  object_name = 'data'))
}
