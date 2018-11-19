# Glaze Formulator: Computes a glaze recipe given a Unity Molecular Formula
# and a list of ingredients.
#     Copyright (C) 2018  Daniel J. Lizotte
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# This function generates the client-side HTML for a basic range input
rangeInput <- function(inputId, label = "", value = 0.5, width = "4em", step = "0.125") {
  tagList(
    # This makes web page load the JS file in the HTML head.
    # The call to singleton ensures it's only included once
    # in a page.
    shiny::singleton(
      shiny::tags$head(
        shiny::tags$script(src = "range-input-binding.js")
      )
    ),
    #shiny::tags$label(label, `for` = inputId, style = "display: inline-block;"),
    shiny::tags$input(id = inputId, type = "range", value = value, min = "0", max = "1", step = step,
                      style = paste0("display: inline-block; vertical-align: middle; width: ", width, ";"))
  )
}

# Send an update message to a range input on the client.
# This update message can change the value and/or label.
updateRangeInput <- function(session, inputId,
                           label = NULL, value = NULL) {
  
  message <- dropNulls(list(label = label, value = value))
  session$sendInputMessage(inputId, message)
}

# Given a vector or list, drop all the NULL items in it
dropNulls <- function(x) {
  x[!vapply(x, is.null, FUN.VALUE=logical(1))]
}