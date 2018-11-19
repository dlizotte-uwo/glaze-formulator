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

library(dplyr)
library(shiny)
library(shinyjs)
library(quadprog)
#library(nnls)
#library(numDeriv)
source("range-input.R")

# Style information: Interface colours for oxides, get form-groups to display as inline-block.
oxide.colour.styles <- HTML("\
                            .form-group {
                            display: inline-block;
                            }

                            .oxide-colors-R2O {
                            color: #d32f2fx
                            }
                            .oxide-colors-RO {
                            color: #1976d2
                            }
                            .oxide-colors-Al2O3,.oxide-colors-PbO,.oxide-colors-SnO2 {
                            }
                            .oxide-colors-B2O3,.oxide-colors-SiO2 {
                            color: #7b1fa2
                            }
                            .oxide-colors-K2O,.oxide-colors-Li2O,.oxide-colors-Na2O {
                            color: #d32f2f
                            }
                            .oxide-colors-BaO,.oxide-colors-BeO,.oxide-colors-CaO,.oxide-colors-MgO,.oxide-colors-SrO {
                            color: #1976d2
                            }
                            .oxide-colors-CdO,.oxide-colors-CoO,.oxide-colors-Cr2O3,.oxide-colors-Cu2O,.oxide-colors-CuO,
                            .oxide-colors-Fe2O3,.oxide-colors-FeO,.oxide-colors-MnO,.oxide-colors-MnO2,.oxide-colors-NiO,
                            .oxide-colors-TiO2,.oxide-colors-V2O5,.oxide-colors-ZnO,.oxide-colors-ZrO,.oxide-colors-ZrO2 {
                            color: #f57c00
                            }
                            .oxide-colors-F {
                            color: #757575
                            }
                            .oxide-colors-P2O5 {
                            color: #0097a7
                            }
                            .oxide-colors-fill-R2O {
                            fill: #d32f2f
                            }
                            .oxide-colors-fill-RO {
                            fill: #1976d2
                            }
                            .oxide-colors-fill-Al2O3,.oxide-colors-fill-PbO,.oxide-colors-fill-SnO2 {
                            fill: #388e3c
                            }
                            .oxide-colors-fill-B2O3,.oxide-colors-fill-SiO2 {
                            fill: #7b1fa2
                            }
                            .oxide-colors-fill-K2O,.oxide-colors-fill-Li2O,.oxide-colors-fill-Na2O {
                            fill: #d32f2f
                            }
                            .oxide-colors-fill-BaO,.oxide-colors-fill-BeO,.oxide-colors-fill-CaO,.oxide-colors-fill-MgO,.oxide-colors-fill-SrO {
                            fill: #1976d2
                            }
                            .oxide-colors-fill-CdO,.oxide-colors-fill-CoO,.oxide-colors-fill-Cr2O3,
                            .oxide-colors-fill-Cu2O,.oxide-colors-fill-CuO,.oxide-colors-fill-Fe2O3,
                            .oxide-colors-fill-FeO,.oxide-colors-fill-MnO,.oxide-colors-fill-MnO2,
                            .oxide-colors-fill-NiO,.oxide-colors-fill-TiO2,.oxide-colors-fill-V2O5,
                            .oxide-colors-fill-ZnO,.oxide-colors-fill-ZrO,.oxide-colors-fill-ZrO2 {
                            fill: #f57c00
                            }
                            .oxide-colors-fill-F {
                            fill: #757575
                            }
                            .oxide-colors-fill-P2O5 {
                            fill: #0097a7
                            }
                            
                            h2 {
                            margin-top: 5px;
                            }
                            ")

#Read in excel sheet from glazy
#Columns are oxides + LOI (percentages by mass)
#Rows are materials
materials <- readxl::read_excel("materials.xlsx")
#Missing values are actually zero
materials[is.na(materials)] <- 0

#Create matrix by removing non-numeric columns and transposing
#Remove headers, transpose
materials_wt_matrix <- t(materials[,c(-1,-2,-3)])
#Columns are now ingredients
colnames(materials_wt_matrix) <- materials$Ingredient

#Convert resulting matrix to data frame. Each row is an oxide.
materials_wt <- as.data.frame(materials_wt_matrix)

#Molecular weights for oxides etc. in g/mol
#Order is important (obviously) must match row order of materials_wt_matrix
molwts <- c(SiO2 = 60.08, Al2O3 = 101.96, B2O3 = 69.6182, MgO = 40.304, CaO = 56.0774,
            SrO = 103.619, BaO = 153.326, ZnO = 81.38, Li2O = 29.88, K2O = 94.20,
            Na2O = 61.98, P2O5 = 283.886/2, Fe2O3 = 159.69, TiO2 = 79.866,
            MnO = 70.9374, MnO2 = 86.9368, `F` = 18.998, PbO = 223.20)

# Identify which oxides should be considered fluxes for purpose of unity formula
fluxes <- c("K2O","Na2O","Li2O","CaO","MgO","ZnO","SrO","BaO")

#Gives number of moles of each oxide in 100g of ingredient
#Omit LOI in mol matrix
materials_mol_matrix <- materials_wt_matrix[-19,] / molwts
#materials_mol_matrix <- rbind(materials_mol_matrix, materials_mol_matrix['Na2O',] + materials_mol_matrix['K2O',])
#rownames(materials_mol_matrix)[19] <- "KNaO"
#Remove materials that have no relevant oxides
materials_mol <- as.data.frame(materials_mol_matrix) %>% 
  dplyr::select(which(sapply(.,function(c) sum(c) > 0)))

#' Returns a list of glaze recipes given a target and a formula.
#' E.g.
#' form = Target ~ -1 + `Frit 3134` + `Vansil W-10 Wollastonite` + 
#'            `EPK, EP Kaolin, Edgar Plastic` + `SIL-CO-SIL 45`
#' Note that -1 ensures no intercept, which makes no sense in this case.
#' Resulting recipes minimize the squared error between Target and Actual 
#' Unity Molecular Formula. Recipes are computed using constrained quadratic programming.
#'
#' @param targ The target UMF oxide values
#' @param form The formula, of the form ~ Ingredient1 + Ingredient2 + ... + Ingredientp
#' @param infl Vector of influences ("importances") for each oxide
#' @param oxideNames List of oxide names corresponding to the target values
#' @param combineKNa Boolean specifying whether to consider moles of K and Na as interchangeable
#'
#' @return A list containing one or more recipes
#'
glazeRecipes <- function(targ, form, infl = setNames(rep(1, length(targ)), oxideNames = names(targ)),
                        combineKNa = FALSE) {
  #"targ" must be a named vector with all appropriate oxides etc.; 
  # these may be in different order from oxides in the model data
  #"form" must be a formula object
  #"infl" must be a named vector with all appropriate oxides etc.; 
  # these may be in different order from oxides in the model data
  
  #Set target in the materials_mol data frame to the given target molecular formula
  #Set influnces as well
  #This ensures that oxides "line up" properly if they're given out of order

  materials_mol$Target = 0
  materials_mol[names(targ),"Target"] <- targ
  materials_mol$Influence = 0
  materials_mol[names(infl),"Influence"] <- infl
  
  #print("Influence")
  #print(infl)
  
  #Create mole matrix for this optimization
  mole_matrix_orig <- model.matrix(form, data = materials_mol)
  target_orig <- materials_mol$Target
  names(target_orig) <- rownames(mole_matrix_orig)
  influence_orig <- materials_mol$Influence
  names(influence_orig) <- rownames(mole_matrix_orig)

  #Make a copy in case we need to change it for K and Na
  mole_matrix <- mole_matrix_orig
  target <- target_orig
  influence <- influence_orig
  opt_fluxes <- fluxes
  
  #Possibly combine K and Na accounting
  if(combineKNa) {
    #print("Combining K and Na")
    mole_matrix <- rbind(mole_matrix, KNaO = mole_matrix["K2O",] + mole_matrix["Na2O",])
    target["KNaO"] <- target["K2O"] + target["Na2O"]
    #If errors for K and Na were equal, setting the KNaO influence
    #to 0.25 times sum of influences gives same total impact to the error function
    influence["KNaO"] <- 0.25*(influence["K2O"] + influence["Na2O"])
    
    mole_matrix <- mole_matrix[!rownames(mole_matrix) %in% c("K2O","Na2O"), , drop = FALSE]
    target <- target[!names(target) %in% c("K2O","Na2O")]
    influence <- influence[!names(influence) %in% c("K2O","Na2O")]
    
    # If combining KNa, must re-define what "fluxes" are.
    opt_fluxes <- c(opt_fluxes[!opt_fluxes %in% c("K2O","Na2O")], "KNaO")
  } else {
    #print("Not Combining K and Na because")
    #print(combineKNa)
  }

  #Computes a recipe or list of recipes for a given matrix of ingredients.
  #mole_matrix_in may have 1 fewer row than mole_matrix_orig_in
  #if K and Na were combined. Results are always in terms of the original oxide list, however.
  #mask indicates which columns (ingredients) of the two matrices should be used.
  #this is in case the problem is under-determined, in which case the function is called recursively with
  #one or more ingredients removed.
  glazeRecipe <- function(mole_matrix_in, mole_matrix_orig_in, mask = rep(TRUE, ncol(mole_matrix_in))) {

    #Remove masked ingredients/columns
    mole_matrix <- mole_matrix_in[, mask, drop = FALSE]
    mole_matrix_orig <- mole_matrix_orig_in[, mask, drop = FALSE]
    
    #Prepare list to hold recipes
    recipe <- list()
    
    #Compute QR factorization to determine whether the problem is underdetermined.
    qrmm <- qr(t(mole_matrix * sqrt(influence)))
    
    #If rank is less than number of columns, then problem is underdetermined.
    if(qrmm$rank < ncol(mole_matrix)) {
      # Find columns in the reduced matrix that are in the null space, hence are collinear.
      nullSet <- if(qrmm$rank == 0L) seq_len(ncol(mole_matrix)) else -seq_len(qrmm$rank)
      nullmm <- qr.Q(qrmm, complete = TRUE)[, nullSet, drop = FALSE]
      # cSet containts all columns that are candidate for removal.
      cSet <- apply(abs(nullmm) > .Machine$double.eps, 1, any)
      #Try removing each column that can be perfectly predicted.
      for(cc in which(cSet)) {
        maskCol <- seq(ncol(mole_matrix_in))[mask][cc]
        #If column is currently included and no "later" columns are currently excluded
        #then we mask that column and call this function recursively
        if(all(mask[maskCol:ncol(mole_matrix_in)])) {
          newMask <- mask
          newMask[maskCol] <- FALSE
          #Try finding a recipe without that column (ingredient)
          #append it to our current list of results
          recipe <- append(recipe,glazeRecipe(mole_matrix_in, mole_matrix_orig_in, newMask))
        }
      }
      #Return the list of recipes found.
      return(recipe)
    }
    
    #If we reach this point, problem is not underdetermined and we will find exactly 1 optimal
    #recipe using the given set of ingredients.
    
    #Quadratic objective measures total squared error between target and recipe UMF
    #weighted by influence
    Dmat <- t(mole_matrix) %*% diag(influence) %*% mole_matrix
    dvec <- t(target) %*% diag(influence) %*% mole_matrix
  
    #Equality constraints ensure that moles of fluxes sum to 1.0
    Amateq <- matrix((names(target) %in% opt_fluxes) %*% mole_matrix, nrow = 1)
    bveceq <- matrix(1)
    meq = 1
    
    #If the equality constraint is impossible to satisfy
    #because there are no fluxes to sum to 1.0, ignore it.
    #recipe will match the relative proportions of other oxides
    #as closely as possible.
    if(!any(Amateq > 0)) {
      Amateq <- NULL
      bveceq <- NULL
      meq = 0
    }
    
    #Inequality constraints ensure no negative recipe amounts
    Amatineq <- diag(ncol(mole_matrix))
    bvecineq <- matrix(rep(0,ncol(mole_matrix)), ncol = 1)
    
    #These matrices combine the equality and inequality constraints.
    Amat <- t(rbind(Amateq, Amatineq))
    bvec <- rbind(bveceq, bvecineq)
  
    #Note meq indicates how many equality constraints are in the problem.
    #(in our case, 1 or 0)
    qp_solution <- solve.QP(Dmat, dvec, Amat, bvec, meq)
  
    #Recover optimal mixture, ensuring no negative amounts caused
    #by rounding error
    optmix <- pmax(qp_solution$solution,0)
    names(optmix) <- colnames(mole_matrix_orig)
    
    #Recover optimal mixture in terms of original ingredient list 
    #(in case any were eliminated because problem was underdetermined)
    optmixorig <- rep(0,ncol(mole_matrix_orig_in))
    names(optmixorig) <- colnames(mole_matrix_orig_in)
    optmixorig[names(optmix)] <- optmix
    
    #Recover optimal unity formula, always in terms of original oxides (not combining K and Na)
    unity <- pmax(mole_matrix_orig_in %*% optmixorig, 0)
    names(unity) <- rownames(mole_matrix_orig_in)
    
    #Normalize recipe to 100%
    weights <- 100 * optmixorig / sum(optmixorig)
  
    #Also give back the original Target, for convenience
    outTarget <- round(target_orig, digits = 3)
    names(outTarget) <- rownames(mole_matrix_orig)

    #Return weights (ingredient amounts), unity (unity formula for computed recipe) 
    #and target (originally desired unity formula), and optimization result
    return(list
           (list(weights = round(weights, digits = 2),
                unity = round(unity, digits = 3),
                target = outTarget,
                optresult = qp_solution)))
  }

  #Computes a list containing one recipe if the optimal solution is unique,
  #or a list of recipes if the problem is underdetermined
  allRecipes <- glazeRecipe(mole_matrix, mole_matrix_orig)
  
  #Some of the recipes from above may not be optimal.
  #Remove any non-optimal recipes.
  recipeValues <- sapply(allRecipes, function(e) {e$optresult$value})
  optRecipes <- (recipeValues - min(recipeValues)) < 1e-8
  allRecipes <- allRecipes[optRecipes]
  
  #Identify any duplicate recipes
  allMixes <- sapply(allRecipes, function(e) {e$weights})
  #Each column is a mix
  duplicatedMixes <- duplicated(allMixes, MARGIN = 2)

  #Return list of recipes, only those that are not
  #duplicates. All are equally optimal.
  return(allRecipes[!duplicatedMixes])
}


# Re-order for interface - fluxes, stabilizers, glass-formers, other
oxides <- rownames(materials_mol)[c(9,11,10,4,5,6,7,8,2,3,1,14,12,13,15,16,17,18)]
# Indices corresponding to fluxes, for convenience
fluxIndices <- which(oxides %in% fluxes)
# Oxide values for initialization. (Just an example.)
initialvalues <- c(0,0.2,0.1,0,0.7,0,0,0,0.4,0.15,3.0,0,0,0,0,0,0,0)
oxidetypes <- c(rep("FluxR2O",3),rep("FluxRO",5),rep("Stabilizer",2),rep("GlassFormer",2),rep("Other",6))

# Colours to indicate unity formula values that are too low, perfect, and too high.
errorColor <- colorRamp(c("#00ffff","#ffffff","#ff0000"))

# Format a number without sci notation, and keep as many digits as possible (do
# we really need to go beyond 15 digits?)
formatNoSci <- function(x) {
  if (is.null(x)) return(NULL)
  format(x, scientific = FALSE, digits = 15)
}

# Convenience function for AND
`%AND%` <- function(x, y) {
  if (!is.null(x) && !is.na(x))
    if (!is.null(y) && !is.na(y))
      return(y)
  return(NULL)
}


# Function to create checkbox input but with no borrom margin.
myCheckboxInput <- function (inputId, label, value = FALSE, width = NULL) 
{
  value <- restoreInput(id = inputId, default = value)
  inputTag <- tags$input(id = inputId, type = "checkbox")
  if (!is.null(value) && value) 
    inputTag$attribs$checked <- "checked"
  div(class = "form-group shiny-input-container", style = if (!is.null(width)) 
    paste0("width: ", validateCssUnit(width), "; margin-bottom: 0px;"), div(class = "checkbox", 
                                                        tags$label(inputTag, tags$span(label))))
}

# Function that creates a UI for the given "oxide" in the unity formula.
oxideInput <- function (inputId, label, value, min = NA, max = NA, step = NA, type = "FluxR2O", width = NULL) 
{
  value <- restoreInput(id = inputId, default = value)

  inputTag <- tags$input(id = inputId, type = "number", class = "form-control", 
                         style = "display: inline-block; width: 6em; height: 2.1em; vertical-align: middle", 
                         value = formatNoSci(value))
  if (!is.na(min)) 
    inputTag$attribs$min = min
  if (!is.na(max)) 
    inputTag$attribs$max = max
  if (!is.na(step)) 
    inputTag$attribs$step = step

  resultTag <- textOutput(outputId = paste0(inputId, "Result"), inline = TRUE)
  resultTag$attribs$style <- "padding-left: 0.5em; display: inline-block; 
    width: 5em; color: black; border-width: 3px; border-style: none; vertical-align: middle;"

  influenceTag <- rangeInput(inputId = paste0(inputId, "Influence"), width = "8em")
  
  tags$div(style = "width: 26em;", 
        div(class = "form-group shiny-input-container",
        style = if (!is.null(width)) paste0("width: ", validateCssUnit(width), "; margin: 0;"),
        influenceTag,
        HTML(paste0(label %AND% tags$label(label, class = paste0("oxide-colors-",label), 
                               `for` = inputId, 
                               style = paste0("width: 4em; margin: 0 0 0 1em; text-align: right; padding-right: 5px; display: inline-block; vertical-align: middle;")), 
                               inputTag)), resultTag))
#  div(class = "form-group shiny-input-container", 
#      style = if (!is.null(width)) paste0("width: ", validateCssUnit(width), "; ", "margin-bottom: 5px;"), 
#     rangeInput(inputId = paste0(inputId, "Influence")),
#      label %AND% tags$label(label, class = paste0("oxide-colors-",label), 
#                             `for` = inputId, style = paste0("width: 4em; text-align: right; padding-right: 5px;")), 
#                             inputTag, resultTag)
}

# Main UI code
ui <- function(request) {
  fluidPage(
    useShinyjs(), # Need Shinyjs to change border colours on the fly to show errors
    tags$head(includeHTML("www/google-analytics.js")),
    tags$style(type = "text/css",oxide.colour.styles), # Set glazy colours for oxides
    tags$style(type = "text/css","@import url('https://fonts.googleapis.com/css?family=Montserrat:400,700');\n* {font-family: Montserrat,'Helvetica Neue',sans-serif;}"),
    title = "Glaze Formulator",
    tags$h1("The Glaze Formulator", style = "margin-top: 5px; margin-bottom: 5px;"), 
    HTML('<p>Welcome to the Glaze Formulator! This application will take a target <a href="https://www.google.com/search?q=unity+molecular+formula" target="_blank">unity molecular formula</a> 
           together with a list of ingredients to use, and it will produce a glaze recipe that approximates the target unity formula as closely as possible.
           Enter your desired target formula in the boxes under <b>Unity Formula</b> (there is a simple example formula to start from), then choose the materials you want 
           using the box under <b>Materials</b> (you can type to search). Your recipe will appear under <b>Recipe</b>, and its unity formula will appear under <b>Actual</b>. 
           Large discrepancies between the <b>Target</b> unity formula and the <b>Actual</b> unity formula will be highlighted. You can approximate the example unity formula pretty well 
           with EPK, Silica, Frit 3134, Neph Sy, and Whiting. Give it a try.</p>'),
    tags$hr(style = "margin-top: 0px; margin-bottom:3px;"),
    # Inputs (Unity formula) on the left, Materials/ingredient choices in the middle, recipe on the right.
    fluidRow(
      column(5, tags$h2("Unity Formula"),
                tags$hr(style = "margin-top: 0px; margin-bottom:3px;"),
                tags$div(style = "width: 100%; display: table;", 
                          tags$div(style = "display: table-cell; vertical-align: middle; width: 8em;", 
                          tags$div(style = "width: 8em;",
                                   HTML(paste0(
                          actionButton(inputId = "ZeroInfluences", label = HTML("<span class=\"glyphicon glyphicon-fast-backward\" aria-hidden=\"true\"></span>"), width = "2em", style = "padding-left: 0; padding-right: 0;"),
                          actionButton(inputId = "ResetInfluences", label = "Reset", width = "4em", style = "padding-left: 0; padding-right: 0;"),
                          actionButton(inputId = "MaxInfluences", label = HTML("<span class=\"glyphicon glyphicon-fast-forward\" aria-hidden=\"true\"></span>"), width = "2em", style = "padding-left: 0; padding-right: 0;"))))),
                          tags$div(style = "width: 18em; display: inline-block; text-align: center;",
                                  actionButton(inputId = "NormalizeFluxes", label = "Normalize", width = "10em", style = "margin: auto; display: block; padding: 6px 6px;"),
                                  myCheckboxInput(inputId = "KNaOBox", label = HTML('K <span class="glyphicon glyphicon-transfer" aria-hidden="true"></span> Na'), width = "6em"))),
                tags$div(style = "width: 26em; display: table;", 
                         span(HTML("<span class=\"glyphicon glyphicon-minus-sign\" aria-hidden=\"true\"></span> Influence <span class=\"glyphicon glyphicon-plus-sign\" aria-hidden=\"true\"></span>"),style = "display: inline-block; text-align: center; width: 8em; font-weight: 350;"),
                         HTML(paste0(
                                  span("Oxide",style = "display: inline-block; text-align: right; width: 4em; margin-left: 1em; padding-right: 5px; font-weight: 700;"),
                                  span("Target", style = "display: inline-block; text-align: center; width: 6em; font-weight: 700;"),
                                  span("Actual", style = "padding-left: 5px; font-weight: 700;")))),
                tags$hr(style = "margin-top: 0px; margin-bottom:3px;"),
                lapply(1:length(oxides), function(i) {
                oxideInput(inputId = oxides[i], label = oxides[i], type = oxidetypes[i], width="60em",
                value = initialvalues[i], min = 0, step = 0.01)})
      ),
      column(3, tags$h2("Materials"),
                tags$hr(style = "margin-top: 0px; margin-bottom:3px;"),
                selectInput(label = "Choose Materials", inputId = "Materials", choices = colnames(materials_mol), width="100%", multiple = TRUE)
             ),
      column(4, tags$h2(textOutput(outputId = "RecipeHeader", inline = TRUE)),
             tags$hr(style = "margin-top: 0px; margin-bottom:3px;"),
             textInput(inputId = "RecipeName", label="Recipe Name", width = "100%"),
             tags$h4(textOutput(outputId = "RecipeStatus", inline = TRUE), style = "color: #f57c00;"),
             tags$h4(textOutput(outputId = "TargetStatus", inline = TRUE), style = "color: #f57c00;"),
             uiOutput("Recipe"),
             tags$hr(style = "margin-top: 0px; margin-bottom:3px;"),
             tags$br(),
             bookmarkButton())),
      tags$hr(),
      HTML("<p><b>For advanced users:</b> Sometimes, a perfect approximation is not possible with a given collection of materials. 
           In that case, a compromise must be made. For example, if you try to approximate 
           the example formula using Frit 3195 instead of 3134, the sodium will be a little low. 
           However, you can tell the program to try harder to match the sodium value at the 
           expense of other oxides by adjusting its <b>Influence</b>.
           If you turn up the influence for sodium (move slider to the right) you will see your recipe change, 
           you will see that the sodium value matches better, but the potassium will become higher than the target. 
           If you move the slider to the left, that oxide will be less important. Move it all the way to the left
           and that oxide will be ignored when determining the recipe.
           The sliders let you determine which target values are most important to you and which ones don't matter. You can also click the 
           checkbox marked K <span class=\"glyphicon glyphicon-transfer\" aria-hidden=\"true\"></span> Na if you want the optimizer to ignore the
           distinction between those two oxides. If the box is checked, the optimizer will try to ensure that the sum of the oxides in the recipe UMF
           matches the sum of the targets for those oxides.</p><p><b>Multiple recipes:</b> for some settings, multiple different recipes may be found
           that have the same UMF. In this case, all recipes will be displayed. Any blend of these recipes will also have the same UMF.</p>")
  )
}

# Round a pre-normalized vector v to digits of precision, ensuring the vector sums exactly to given total
exactRound <- function(v, digits, total = 1) {
  roundFactor <- 10^digits
  newv <- floor(v * roundFactor)
  #Ensure sums to to roundFactor by adjusting largest value (causes minimum change of ratios)
  mi <- which.max(newv)
  newv[mi] <- newv[mi] + (total*roundFactor - sum(newv))
  newv <- newv / roundFactor
}

# Main server code
server <- function(input, output, session) {
  
  # Reactive variable that gets the desired influence of each oxide from the sliders
  influences <- debounce(reactive({
    t <- sapply(seq(length(oxides)), function(i) as.numeric(input[[paste0(oxides[i], "Influence")]]))
    names(t) <- oxides
    # "Raw" influences from sliders range from 0 to 1
    # These are transformed to a more useful range, with slow increase from 0 to 1 if slider,
    # is left of middle, then fast increase from 1 to 100 if slider is to right of middle. 
    transform <- function(i) { if(i <= 0.5) {2*i} else {1 + (2*i - 1)*99} }
    tt <- sapply(t, transform)
    tt
  }), 100) # Debounce of 100ms
  
  # Reactive variable that gives the target oxide values from the input boxes
  # Data validation sets non-numeric or negative inputs to 0 (internally and in UI)
  targetOxides <- debounce(reactive({
    t <- sapply(seq(length(oxides)), function(i) input[[oxides[i]]])
    names(t) <- oxides
    # Update invalid inputs to zero
    invalidInputs <- which(is.na(t) | t < 0)
    lapply(invalidInputs, function(i) updateNumericInput(session, oxides[i], value = 0))
    t[invalidInputs] <- 0
    t
  }), 100) # Debounce of 100ms
  
  # Reactive variable that is sum of target fluxes. Used to watch for non-normalized inputs.
  targetFluxSum <- reactive(sum(targetOxides()[fluxIndices]))
  
  # Watch the targetFluxSum, change label of the Normalize Fluxes button if inputs become unnormalized.
  observe(
    if(targetFluxSum() == 0) {
      updateActionButton(session, inputId = "NormalizeFluxes", label = "<span class=\"glyphicon glyphicon-alert\" aria-hidden=\"true\"></span> Normalize <span class=\"glyphicon glyphicon-alert\" aria-hidden=\"true\"></span>")
      output$TargetStatus <- renderText("Flux targets must sum to a positive number.")
    } else if(targetFluxSum() != 1.0) {
      updateActionButton(session, inputId = "NormalizeFluxes", label = "<span class=\"glyphicon glyphicon-alert\" aria-hidden=\"true\"></span> Normalize <span class=\"glyphicon glyphicon-alert\" aria-hidden=\"true\"></span>")
      output$TargetStatus <- renderText("")
    } else {
      updateActionButton(session, inputId = "NormalizeFluxes", label = "Normalize")
      output$TargetStatus <- renderText("")
    }
  )
  
  # Renormalize if the Normalize Fluxes button is clicked.
  # Updates the input boxes if necessary, which invalidates the targetOxides and targetFluxSum
  observeEvent(input$NormalizeFluxes,
               {
                 if(targetFluxSum() > 0) {
                   #Get inputs
                   oxideTargets <- as.numeric(sapply(seq(length(oxides)), function(i) input[[oxides[i]]]))
                   #Normalize
                   newTargets <- oxideTargets / targetFluxSum()
                   #Round to 3 decimal places
                   newTargets[fluxIndices] <- exactRound(newTargets[fluxIndices], 3)
                   #Update the inputs with normalized values
                   lapply(fluxIndices, function(i) updateNumericInput(session, oxides[i], value = newTargets[i]))
                   output$TargetStatus <- renderText("")
                 } else {
                   output$TargetStatus <- renderText("Flux targets must sum to a positive number.")
                 }
               })
  
  # Reset influences if Reset button is clicked.
  observeEvent(input$ResetInfluences,
               {
                   lapply(seq(length(oxides)), function(i) updateRangeInput(session, paste0(oxides[i], "Influence"), value = 0.5))
               })

  # Zero influences if zero button is clicked.
  observeEvent(input$ZeroInfluences,
               {
                 lapply(seq(length(oxides)), function(i) updateRangeInput(session, paste0(oxides[i], "Influence"), value = 0))
               })

  # Max influences if max button is clicked.
  observeEvent(input$MaxInfluences,
               {
                 lapply(seq(length(oxides)), function(i) updateRangeInput(session, paste0(oxides[i], "Influence"), value = 1))
               })
  
  # Recipe is a reactive variable, depends on the targetOxides() reactive variable and on input$Materials
  recipe <- reactive(
    if(!is.null(input$Materials) && any(influences() > 0)) {
      #Create formula string
      allMaterials <- paste(input$Materials, collapse = "` + `")
      #Create formula from formula string
      form <- as.formula(paste0("Target ~ -1 + `", allMaterials, "`"))
      #Use the glazeRecipes function to find optimal recipes given the
      #list of materials and the unity targets.
      r <- glazeRecipes(targetOxides(), form, infl = influences(), combineKNa = input$KNaOBox)
      if(is.null(r) || any(!is.finite(r[[1]]$unity))) {
        # Something has gone very wrong.
        output$RecipeStatus <- renderText("Unable to compute recipe.")
        print(r)
        NULL
      } else if (length(r) > 1) {
        # Underdetermined problem; multiple recipes found
        output$RecipeStatus <- renderText("Multiple recipes found with same unity formula.")
        r
      } else {
        # Seems fine
        output$RecipeStatus <- renderText("")
        r
      }
    } else {
      if(is.null(input$Materials)) {
        output$RecipeStatus <- renderText("Select materials to compute recipe.")
      } else if (all(influences() == 0)) {
        output$RecipeStatus <- renderText("At least one influence must be non-zero.")
      }
      NULL # If no Materials entered or influences all 0, return NULL
    }
  )
  
  # Render header for recipe
  output$RecipeHeader <- renderText("Recipe")
  
  # Update the actual unity formula of the recipe if it changes
  lapply(oxides, function(ox) { output[[paste0(ox, "Result")]] <- renderText(if(is.null(recipe())) "" else recipe()[[1]]$unity[[ox]]) })
  
  # Set the borders of the displayed unity formula to reflect errors.
  # If oxide is above target, colour is red. If below, colour is cyan.
  observe({
              jsCommand <- ""
              if(!is.null(recipe())) {
                #Compute errors and display, along with coloured borders
                for(o in oxides) {
                 oResult <- paste0(o, "Result")
                 u <- recipe()[[1]]$unity[o]
                 t <- recipe()[[1]]$target[o]
                 # Convert errors to colour. This is subjective.
                 err = (u - t)
                 err1 <- abs(err / (u + t + 0.01))*(abs(err) > 0.01)*sign(err)
                 err2 <- sqrt(abs(err1)) * sign(err1)
                 c <- rgb(errorColor(err2 / 2 + 0.5), max = 255)
                 #Adjust colours based on error.
                 jsCommand <- paste0(jsCommand, "$('#",oResult,"').css({'border-style':'solid','border-color':'",c,"'});")
                }
                runjs(jsCommand)
              } else {
                #Recipe is NULL; blank out unity formula and turn off borders.
                for(o in oxides) {
                  oResult <- paste0(o, "Result")
                  jsCommand <- paste0("$('#",oResult,"').css({'border-style':'none','border-color':'white'});")
                  runjs(jsCommand)
                }
              }
    })
  
  # Show the actual computed recipes in Column 3.
  output$Recipe <- renderUI({
    tables <- lapply(recipe(), function(r) {renderTable({
      recipeDF <- data.frame(Amount = exactRound(r$weights, 2, 100))
      # Trim off any initial or trailing back-ticks from column names.
      rownames(recipeDF) <- gsub("`$","",gsub("^`","",names(r$weights)))
      recipeDF
      }, rownames = TRUE, width = "100%")})
    list(tables)
  })
}

shinyApp(ui, server, enableBookmarking = "url")
