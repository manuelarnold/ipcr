#' @noRd
# Function to add significance stars based on p-values
add_stars <- function(p) {
  if (is.na(p)) {return(NA)}
  else if (p < 0.001) return("***")
  else if (p < 0.01) return("**")
  else if (p < 0.05) return("*")
  else if (p < 0.1) return(".")
  else return("")
}
