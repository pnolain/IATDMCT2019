library(shiny)
library(DT)
library(rhandsontable)
library(tidyverse)
library(mrgsolve)

# Helper functions
convert_time <- function(time, unit_in, unit_out){
  time_units <- c("minutes" = 1 / 60, "hours" = 1, "days" = 24, "weeks" = 24 * 7)

  conversion_factor <- as.numeric(time_units[unit_in] / time_units[unit_out])

  time * conversion_factor
}

# wrapper function around shiny::req, preventing conflits with mrgsolve::req
shinyreq <- shiny::req

#Define user-input dependent functions for output
shinyServer(function(input, output, session) {

  # User explicitly defined widgets
  custom_widgets <- list(
    WT = list(widget = "slider", label = "Weight", min = 35, max = 250, step = 1),
    AGE = list(widget = "slider", label = "Age", min = 18, max = 90, step = 1),
    SEX = list(widget = "select", label = "Gender", choices = c("Male" = 0, "Female" = 1)),
    STATIN = list(widget = "select", label = "Statin co-administration", choices = c("No" = 0, "Yes" = 1)),
    HDSTATIN = list(widget = "select", label = "High dose statin", choices = c("No" = 0, "Yes" = 1)),
    MONO = list(widget = "select", label = "Monotherapy", choices = c("No" = 0, "Yes" = 1)),
    BSLDLC = list(widget = "numeric", label = "Baseline LDL-C", min = 0, max = 400, step = 1),
    TBSPCSK = list(widget = "numeric", label = "Baseline Total PCSK9", min = 0, max = NA, step = 0.1),
    TPCSK = list(widget = "numeric", label = "Total PCSK9", min = 0, max = NA, step = 0.1),
    FBSPCSK = list(widget = "numeric", label = "Baseline Free PCSK9", min = 0, max = NA, step = 0.1),
    FPCSK = list(widget = "numeric", label = "Free PCSK9", min = 0, max = NA, step = 0.1),
    DISST = list(widget = "select", label = "Disease status", choices = c("Healthy volunteer" = 0, "Patient" = 1)))
  
  # mrgsolve model
  pkpd <- mrgsolve::mread_cache("PKPD", project = "www")
  tmdd <- mrgsolve::mread_cache("TMDD", project = "www")

  withProgress({
    # C++ function for AUC fast computation
    Rcpp::sourceCpp("www/auc.cpp")
  }, message = "Startup compilation", detail = "Please wait...")
  
  pkpd_annotations <- pkpd@annot$data %>%
    as_tibble() %>%
    mutate(unit = ifelse(unit == ".", NA, unit),
           options = ifelse(options == ".", NA, options))
  
  tmdd_annotations <- tmdd@annot$data %>%
    as_tibble() %>%
    mutate(unit = ifelse(unit == ".", NA, unit),
           options = ifelse(options == ".", NA, options))

  covariates_names <- bind_rows(list(pkpd = pkpd_annotations,
                                     tmdd = tmdd_annotations), .id = "model") %>%
    filter(options == "covariate")

  pkpd_cov_names <- covariates_names %>% filter(model == "pkpd")
  tmdd_cov_names <- covariates_names %>% filter(model == "tmdd")

  # model instances----
  studies <- list(
    PKPD = list(
      model = pkpd, 
      problem = pkpd@code[2],
      time_unit = "hours",
      routes = c(SC = 1, IV = 2),
      annotations = pkpd_annotations,
      outputs = pkpd_annotations %>% filter(options == "simshiny"),
      parameters = pkpd@param@data[which(!(names(pkpd@param@data) %in% c(pkpd_cov_names$name)))],
      covariates = pkpd@param@data[pkpd_cov_names$name]
      ),
    TMDD = list(
      model = tmdd,
      problem = tmdd@code[2],    time_unit = "days",
      routes = c(SC = 1, IV = 2),
      annotations = tmdd_annotations,
      outputs = tmdd_annotations %>% filter(options == "simshiny"),
      parameters = tmdd@param@data[which(!(names(tmdd@param@data) %in% c(tmdd_cov_names$name)))],
      covariates = tmdd@param@data[tmdd_cov_names$name]
  ))
  
  
  # Generate application widgets (covariates/parameters)
  param_covs_annotations <- map(studies, function(x){
    x$annotations %>%
      filter(block == "PARAM") %>%
      mutate(widget = pmap(list(name, descr, unit), function(name, descr, unit){
        custom <- custom_widgets[[name]]
        
        wid <- list(widget = "numeric", width = "200px")
        
        if(!is.na(descr)){
          wid$label <- ifelse(is.na(unit), descr, sprintf("%s (%s)", descr, unit))
        } else {
          wid$label <-  ifelse(is.na(unit), name, sprintf("%s (%s)", name, unit))
        }
        
        if(!is.null(custom))
          wid[names(custom)] <- custom[names(custom)]
        
        wid
      }))
    
    })
  
  app_widgets <- map(param_covs_annotations,  ~ set_names(.$widget, .$name))

  # Application specific variables----
  design_dosing_unit <- "mg"
  available_time_units <- c("days", "weeks")
  default_time_unit <- "days"
  default_precision <- 6/24
  default_dosing_interval <- 14
  default_last_time <- 180
  max_last_time <- 365

  default_design <- list(
    start = 0,
    end = default_last_time,
    precision = default_precision,
    time_unit = default_time_unit,
    administrations = tibble(
      time = numeric(),
      dose = numeric(),
      route = character()
    )
  )

  default_study <- studies$PKPD

  # Reactive Values----
  v <- reactiveValues(
    design = default_design,
    selected_study = default_study,
    history = tibble(sim = list()),
    reference_lines = tibble(output = character(), label = character(), value = numeric()),
    x = NULL,
    y = NULL,
    computation_duration = 0,
    computation_times = 0,
    reset_history_list = 0
  )

  # Hide the loading message when the rest of the server function has executed
  hide(id = "loading-content", anim = TRUE, animType = "fade", time = 2.5)

  # Functions----
  update_run_name <- function() {
    if (nrow(v$design$administrations) == 0) {
      updateTextInput(session, "run_name", value = "")
    } else {
      first_admin <-  v$design$administrations %>% slice(1)

      # add first administration
      regimens <-
        tibble(
          first_time = first_admin$time,
          dose = first_admin$dose,
          route = first_admin$route,
          interval = NA,
          N = 1
        )

      if (nrow(v$design$administrations) > 1)
      {
        current_regimen_id <- 1
        previous_admin <- first_admin

        for (i in 2:nrow(v$design$administrations))
        {
          current_regimen <- regimens[current_regimen_id, ]
          current_admin <- v$design$administrations %>% slice(i)

          new_regimen <- NULL

          if ((current_admin$dose != current_regimen$dose) |
              (current_admin$route != current_regimen$route))
          {
            new_regimen = tibble(
              first_time = current_admin$time,
              dose = current_admin$dose,
              route = current_admin$route,
              interval = NA,
              N = 1
            )

          } else if (is.na(current_regimen$interval)) {
            regimens[current_regimen_id, ]$interval = current_admin$time - previous_admin$time

          } else if (((current_admin$time - previous_admin$time) != current_regimen$interval)) {
            new_regimen = tibble(
              first_time = current_admin$time,
              dose = current_admin$dose,
              route = current_admin$route,
              interval = NA,
              N = 1
            )
          }

          if (!is.null(new_regimen))
          {
            regimens <- bind_rows(regimens, new_regimen)
            current_regimen <- new_regimen
            current_regimen_id <- current_regimen_id + 1
          } else {
            regimens[current_regimen_id, ]$N <-
              regimens[current_regimen_id, ]$N + 1
          }

          previous_admin = current_admin
        }
      }

      singular_dtu <-
        substr(v$design$time_unit, 1, nchar(v$design$time_unit) - 1)

      admin_groups <- regimens %>%
        mutate(regimen = paste(
          dose,
          design_dosing_unit,
          paste0("(", route,
                 if_else(
                   N > 1, paste(
                     " every",
                     round(interval, 2),
                     if_else(interval == 1, singular_dtu, v$design$time_unit)
                   ), ""
                 )),
          "-",
          paste(N, paste0("dose", ifelse(N > 1, "s)", ")")))
        ))

      design_name <- paste(admin_groups$regimen, collapse = ", ")

      updateTextInput(session, "run_name", value = design_name)
    }
  }

  convert_design_time_unit <- function(design, new_unit) {
    old_unit <- design$time_unit

    design$start <- convert_time(design$start, old_unit, new_unit)
    design$end <- convert_time(design$end, old_unit, new_unit)
    design$precision <- convert_time(design$precision, old_unit, new_unit)

    design$administrations$time <- convert_time(design$administrations$time, old_unit, new_unit)
    design$time_unit <- input$design_time_unit

    design
  }

  run_simulation <- function()
  {
    sim_design <- v$design

    sim_design$start <- shinyreq(input$design_time_range[1])
    sim_design$end <- shinyreq(input$design_time_range[2])
    sim_design$precision <- convert_time(shinyreq(input$precision), "hours", v$design$time_unit)

    admins <- v$design$administrations

    parameters_panel <- "parameters_ui_output"
    covariates_panel <- "covariates_ui_output"

    temp_parameters <- v$selected_study$parameters
    temp_covariates <- v$selected_study$covariates

    for (p in names(v$selected_study$parameters))
    {
      ui_elem <-
        input[[getUIElementId(parameters_panel, v$selected_study$model, p)]]

      if (!is.null(ui_elem))
        temp_parameters[[p]] <- as.numeric(ui_elem)
    }

    if (!is.null(v$selected_study$covariates)) {
      for (c in names(v$selected_study$covariates))
      {
        ui_elem <-
          input[[getUIElementId(covariates_panel, v$selected_study$model, c)]]

        if (!is.null(ui_elem))
          temp_covariates[[c]] <- as.numeric(ui_elem)
      }
    }

    run_number <- nrow(v$history) + 1

    sim <- NULL

    tryCatch({
      duration <- system.time({

        sim_model <- v$selected_study$model
        
        mw_factor <- ifelse(input$selected_model == "TMDD", 6.8493, 1) # to nM

        ds <- admins %>%
          mutate(cmt = plyr::revalue(route, v$selected_study$routes, warn_missing = FALSE) %>% as.numeric,
                 evid = 1,
                 time = time %>% convert_time(sim_design$time_unit, v$selected_study$time_unit)) %>%
          rename(amt = dose) %>%
          mutate(amt = mw_factor * amt) %>% 
          select(time, amt, evid, cmt)
        
        if(nrow(ds) == 0){
          ds <- ds %>% add_row(time = 0, amt = 0, evid = 1, cmt = 1)
        }

        if(input$simulate_iiv){
          mrg_sim <- sim_model %>%
            param(c(temp_parameters, temp_covariates)) %>%
            ev(as.ev(ds)) %>%
            mrgsim(start = sim_design$start %>% convert_time(sim_design$time_unit, v$selected_study$time_unit),
                   end = sim_design$end %>% convert_time(sim_design$time_unit, v$selected_study$time_unit),
                   delta = sim_design$precision %>% convert_time(sim_design$time_unit, v$selected_study$time_unit),
                   nid = input$iiv_n)
        } else {
          mrg_sim <- sim_model %>%
            param(c(temp_parameters, temp_covariates)) %>%
            ev(as.ev(ds)) %>%
            zero_re() %>% 
            mrgsim(start = sim_design$start %>% convert_time(sim_design$time_unit, v$selected_study$time_unit),
                   end = sim_design$end %>% convert_time(sim_design$time_unit, v$selected_study$time_unit),
                   delta = sim_design$precision %>% convert_time(sim_design$time_unit, v$selected_study$time_unit))
        }
      
        sim_outputs <- mrg_sim@data %>% as_tibble() %>%
          select(ID, time, one_of(v$selected_study$outputs$name)) %>%  as_tibble %>% unique

        if(sim_design$time_unit != v$selected_study$time_unit)
          sim_outputs$time <- sim_outputs$time %>% convert_time(v$selected_study$time_unit, sim_design$time_unit)

        run_name <- ifelse(input$run_name != "", paste(run_number, "-", input$run_name), as.character(run_number))
        
        if(input$simulate_iiv)
          run_name <- str_c(run_name, " [N=", input$iiv_n, "]")

        sim <- tibble(simulation = run_name,
                      model = list(sim_model),
                      design = list(tibble(start = sim_design$start, end = sim_design$end, precision = sim_design$precision,
                                           time_unit = sim_design$time_unit)),
                      administrations = list(admins),
                      outputs = list(sim_outputs),
                      iiv = input$simulate_iiv)

      })

      v$computation_duration <- duration["elapsed"]
      v$computation_times <- nrow(sim_outputs)

      print(duration)

    },
    error = function(err)
    {
      print(err$message)

      shinytoastr::toastr_error(paste("The following error occured:", err$message, sep = "\n"), title = "?!",
                                position = "top-full-width")

    }, finally = {
      return(sim)
    })
  }

  convert_simulation_time_unit <- function(simulation, new_unit){
    sim_design <- simulation$design[[1]]
    sim_outputs <- simulation$outputs[[1]]

    sim_outputs$time <- convert_time(sim_outputs$time, sim_design$time_unit, new_unit)

    simulation$outputs[[1]] <- sim_outputs

    simulation
  }


  # Defines a unique identifier for a widget : "[PARENTOBJECTID]_[MODELCLASSNAME]_[ELEMENTNAME]"
  getUIElementId <- function(parent_id, model, elementName) {
    return(paste(parent_id, model@model, elementName, sep = "_"))
  }
  
  observe({
    
    # shinyjs::toggle(id = "covariates_ui_output", condition = !input$simulate_iiv)
    shinyjs::toggleState(id = "iiv_n", condition = input$simulate_iiv)
    
    shinyjs::toggle(id = "show_iiv", condition = (!is.null(v$history) && (nrow(v$history) > 0) &&
                                                         any(unnest(v$history)$iiv)))
    
  })
  
  # UI events -------------
  observeEvent(input$selected_model, {
    
    study <- studies[[input$selected_model]]
    
    v$selected_study <- study
    clear_admins()
    clear_history()
  })

  # UI events -------------
  observeEvent(input$design_time_unit, {
    if (!is.null(v$x)) {
      v$x <- convert_time(v$x, v$design$time_unit, input$design_time_unit)
    } else {
      v$x <-
        convert_time(input$design_time_range, v$design$time_unit, input$design_time_unit)
    }

    if (!is.null(input$design_time_range)) {
      # Design time range
      max_val <- max(input$design_time_range[2], convert_time(max_last_time, default_time_unit, v$design$time_unit))

      slider_step <- ifelse(input$design_time_unit == "weeks", 0.005,
                            ifelse(input$design_time_unit == "days", 0.05,
                                   ifelse(input$design_time_unit == "hours", 1, 5)))

      updateSliderInput(session, "design_time_range",
                        label = sprintf("Time range (%s)", input$design_time_unit),
                        value = convert_time(input$design_time_range, v$design$time_unit,input$design_time_unit),
                        max = convert_time(max_val, v$design$time_unit, input$design_time_unit),
                        step = slider_step
      )

      # Dosing interval
      updateNumericInput(session, "dosing_interval",
                         label = sprintf("Dosing interval (%s)", input$design_time_unit),
                         value = convert_time(input$dosing_interval, v$design$time_unit, input$design_time_unit))

      # Expo time range
      updateNumericInput(session, "exposition_start",
                         value = convert_time(input$exposition_start, v$design$time_unit, input$design_time_unit))

      updateNumericInput(session, "exposition_end",
                         value = convert_time(input$exposition_end,v$design$time_unit,input$design_time_unit))
    }

    v$design <- convert_design_time_unit(v$design, input$design_time_unit)

    update_run_name()
  })

  # Administrations ----
  observeEvent(input$add_administrations, {
    dose <- shinyreq(input$dose) %>% as.numeric()
    route <- input$route
    n_admins <- shinyreq(input$n_admins)
    dosing_interval <- shinyreq(input$dosing_interval)

    current_admins <- v$design$administrations

    first_admin_time <- ifelse(nrow(current_admins) > 0, max(current_admins$time) + dosing_interval, 0)

    new_regimen <- current_admins %>%
      add_row(time = seq(first_admin_time, by = dosing_interval, length.out = n_admins),
              dose = rep(dose, n_admins),
              route = route)

    v$design$administrations <- new_regimen

    # Extend time range?
    min_value = min(min(v$design$administrations$time), input$design_time_range[1])
    max_value = max(max(v$design$administrations$time) + dosing_interval, input$design_time_range[2])

    updateSliderInput(session, "design_time_range", value = c(min_value, max_value))

    update_run_name()
  })

  observeEvent(input$remove_administrations, {
    selected_administrations <- as.numeric(input$administrations_table_rows_selected)

    if (length(selected_administrations) > 0)
      v$design$administrations <- v$design$administrations %>% slice(-selected_administrations)

    update_run_name()
  })

  clear_admins <- function(){
    v$design$administrations <- v$design$administrations %>% slice(0)
    
    update_run_name()
  }
  
  observeEvent(input$clear_administrations, {
    clear_admins()
  })

  # Simulation ----
  observeEvent(input$simulate_button, {
    sim <- run_simulation()

    v$x <- v$y <- NULL

    if (!is.null(sim)) {
      v$history <- v$history %>% add_row(sim = list(sim))

      hist_df <- v$history %>% unnest(sim) %>% slice(-n())

      runs <- set_names(seq_along(hist_df$simulation), hist_df$simulation)

      if(length(runs) > 0){
        updateCheckboxGroupInput(session, "previous_simulations_selection",
                                 label = "Select previous simulations to display",
                                 choices = runs, selected = NULL)

        freezeReactiveValue(input, "previous_simulations_selection")

      } else {
        v$reset_history_list <- v$reset_history_list + 1
        }

      updateNumericInput(session, "exposition_start", value = sim$design[[1]]$start)
      updateNumericInput(session, "exposition_end", value = sim$design[[1]]$end)
    }
  })

  observeEvent(input$clear_history, {
    clear_history()
  })
  
  clear_history <- function(){
    v$history <- tibble(sim = list())
    v$reset_history_list <- v$reset_history_list + 1
  }

  # Reference lines ----
  observeEvent(input$add_reference_line, {
    ref_line_output <- shinyreq(input$selected_output)
    ref_line_label <- shinyreq(input$reference_line_label)
    ref_line_value <- shinyreq(input$reference_line_value)

    v$reference_lines <- v$reference_lines %>%
      add_row(output = ref_line_output,
              label = ref_line_label,
              value = ref_line_value) %>%
      group_by(output, value) %>%
      slice(n()) %>%
      ungroup()
  })

  observeEvent(input$clear_reference_lines, {
    v$reference_lines <- tibble(output = character(), label = character(), value = numeric())
  })

  # Reactives ----
  get_selected_history <- reactive({
    selected_history <- input$previous_simulations_selection

    if (length(selected_history) == 0)
      return(NULL)

    v$history[as.numeric(selected_history),] %>% unnest(sim)
  })

  app_outputs_table <- reactive({

    shinyreq(v$history)

    sim <- v$history %>% slice(n()) %>% unnest(sim)
    
    shinyreq(nrow(sim) > 0)
    
    sim <- convert_simulation_time_unit(sim, input$design_time_unit)

    selected_history <- get_selected_history()

    sim_outputs <- sim

    if(!is.null(selected_history)){
      for(i in seq_len(nrow(selected_history)))
        selected_history[i,] <-  convert_simulation_time_unit(selected_history[i,], input$design_time_unit)

      sim_outputs <- bind_rows(selected_history, sim)
    }

    exposition_start <- shinyreq(input$exposition_start)
    exposition_end <- shinyreq(input$exposition_end)

    out_table <- sim_outputs %>%
      select(simulation, outputs) %>%
      unnest(outputs) %>%
      filter(between(time, exposition_start, exposition_end))

    out_table
  })

  app_expo_table <- reactive({
    if (!is.null(v$history))
    {
      out_table <- app_outputs_table()

      grouped <- out_table %>%
        gather(output, value, -ID, -simulation, -time) %>%
        group_by(simulation, output, time) %>% 
        summarise(value = median(value, na.rm = TRUE))

      auc <- grouped %>% summarise(AUC = auc_cpp(time, value))

      cc <- grouped %>%
        arrange(value, desc(time)) %>%
        summarise(Cmax = last(value),
                  Tmax = last(time),
                  Cmin = first(value),
                  Tmin = first(time))

      expo_table <- left_join(auc, cc, by = c("simulation", "output"))

    } else {
      expo_table <-
        tibble(output = character(), AUC = numeric(), Tmax = numeric(), Cmax = numeric(), Tmin = numeric(), Cmin = numeric())
    }
  })

  # UI widgets ----
  output$model_problem <- renderText({
    v$selected_study$problem
  })

  output$design_time_unit <- renderUI({
    selectInput("design_time_unit", "Time unit",  choices = available_time_units, selected = default_time_unit)
  })

  output$design_time_range <- renderUI({
    sliderInput(
      "design_time_range",
      label = sprintf("Time range (%s)", default_time_unit),
      min = 0,
      max = max_last_time,
      value = c(0, default_last_time)
    )
  })

  output$precision <- renderUI({
    min_step <- 1 / 12 # 5 minutes
    max_step <- convert_time(1, "days", "hours")
    default_step <-
      convert_time(default_precision, default_time_unit, "hours")

    numericInput("precision",
                 "Time step (hours)",
                 value = default_step,
                 step = 1)
  })

  
  output$route <- renderUI({
    selectInput("route", "Administration route", choices = names(v$selected_study$routes))
  })

  output$covariates_ui_output <- renderUI({
    model_covariates <- v$selected_study$covariates

    if (length(model_covariates) > 0)
    {
      mod_widgets <- app_widgets[[input$selected_model]]
      
      covar_ui <- map(names(model_covariates),
                      function(cov_name) {
                        cov_widget <- mod_widgets[[cov_name]]
                        w_id <- getUIElementId("covariates_ui_output", v$selected_study$model, cov_name)

                        if (!is.null(cov_widget)) {
                          wid_func <- switch(cov_widget$widget,
                                             "slider" = sliderInput,
                                             "numeric" = numericInput,
                                             "select" = selectInput)
                          
                          value_name <- switch(cov_widget$widget,
                                               "select" = "selected",
                                               "value")
                          
                          if(is.null(cov_widget[[value_name]]))
                            cov_widget[[value_name]] <- model_covariates[[cov_name]]

                          do.call(wid_func, c(inputId = w_id, cov_widget[-1]))
                        } else {
                          numericInput(inputId = w_id,
                                       label = cov_name,
                                       value = model_covariates[[cov_name]])
                        }
                      })
    } else {
      covar_ui <- p(em("No covariates associated with the model."))
    }

    div(class = "shiny-flow-layout", covar_ui)
  })

  output$parameters_ui_output <- renderUI({
    input$reset_parameters

    model_parameters <- v$selected_study$parameters
    
    mod_widgets <- app_widgets[[input$selected_model]]
    

    param_ui <- map(names(model_parameters),
                    function(param_name) {
                      param_widget <- mod_widgets[[param_name]]
                      p_id <- getUIElementId("parameters_ui_output", v$selected_study$model, param_name)

                      if (!is.null(param_widget))
                      {
                        target_name <- ifelse(param_widget$widget == "select", "selected", "value")
                        
                        param_widget[[target_name]] <- model_parameters[[param_name]]

                        wid_func <- switch(param_widget$widget,
                                           "slider" = sliderInput,
                                           "numeric" = numericInput,
                                           "select" = selectInput)

                        do.call(wid_func, c(inputId = p_id, param_widget[-1]))
                      } else {
                        numericInput(inputId = p_id,
                                     label = param_name,
                                     value = model_parameters[[param_name]])
                      }
                    })

    div(class = "shiny-flow-layout", param_ui)
  })

  output$selected_output <- renderUI({
    output_choices <- set_names(v$selected_study$outputs$name, sprintf("%s (%s)", v$selected_study$outputs$descr, v$selected_study$outputs$unit))

    selectInput(inputId = "selected_output", label = "Output", choices = output_choices)

  })

  # Design administrations table
  output$administrations_table <- renderDataTable({
    is_empty <- nrow(v$design$administrations) == 0

    if (!is_empty) {
      admins <- v$design$administrations
    } else {
      admins = tibble(
        time = numeric(),
        dose = numeric(),
        route = factor()
      )
    }

    colnames(admins) <- c(
      sprintf("Time (%s)", v$design$time_unit),
      sprintf("Dose (%s)", design_dosing_unit),
      "Route"
    )

    datatable(admins, options = list(dom = "rtp",
                                     scrollY = ifelse(is_empty, NA, "250px"),
                                     scrollCollapse = !is_empty,
                                     paging = F,
                                     language = list(emptyTable = "No administration"), rownames = T))
  })

  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$simulation_plot_dblclick, {
    brush <- input$simulation_plot_brush
    if (!is.null(brush)) {
      v$x <- c(brush$xmin, brush$xmax)
      v$y <- c(brush$ymin, brush$ymax)

    } else {
      v$x <- NULL
      v$y <- NULL
    }
  })

  # Previous simulation runs
  output$previous_simulations_selection <- renderUI({
    v$reset_history_list
    checkboxGroupInput(inputId = "previous_simulations_selection", label = "No simulation history available", choices = NULL)
  })

  # Simulations plot
  output$simulation_plot <- renderPlot({
    shinyreq(nrow(v$history) > 0)

    selected_output <- shinyreq(input$selected_output)
    annotated_output <- v$selected_study$outputs %>% filter(name == selected_output)

    withProgress({

    ref_lines <- v$reference_lines[[selected_output]]

    sim <- v$history %>% slice(n()) %>% unnest(sim)
    sim <- convert_simulation_time_unit(sim, input$design_time_unit)

    selected_history <- NULL

    selected_history <- get_selected_history()

    sim_outputs <- sim

    if(!is.null(selected_history)){
      for(i in seq_len(nrow(selected_history)))
        selected_history[i,] <-  convert_simulation_time_unit(selected_history[i,], input$design_time_unit)

      sim_outputs <- bind_rows(selected_history, sim)
    }
    
    sim_df_noiiv <- sim_outputs %>% 
      filter(!iiv)
    
    sim_df_iiv <- sim_outputs %>% 
      filter(iiv)

    sim_plot <- ggplot()
    
    if(nrow(sim_df_noiiv) > 0){
      sim_df_noiiv <- sim_df_noiiv %>% 
        select(simulation, outputs) %>%
        unnest(outputs) %>%
        select(simulation, ID, time, one_of(selected_output))
      
      sim_plot <- sim_plot+geom_line(data = sim_df_noiiv, 
                                     mapping = aes_string(x = "time", y = selected_output, colour = "simulation"),
                                     size = 0.8)
    }
    
    if(nrow(sim_df_iiv) > 0){
      sim_df_iiv <- sim_df_iiv %>% 
        select(simulation, outputs) %>%
        unnest(outputs) %>%
        select(simulation, ID, time, one_of(selected_output))
      
      sim_df_iiv_summary <- sim_df_iiv %>% 
        group_by(simulation, time) %>% 
        summarise(median = median(!!sym(selected_output), na.rm = TRUE),
                  q05 = quantile(!!sym(selected_output), probs = 0.05, na.rm = TRUE),
                  q95 = quantile(!!sym(selected_output), probs = 0.95, na.rm = TRUE)) %>% 
        ungroup()
      
      sim_plot <- sim_plot+geom_line(data = sim_df_iiv_summary, 
                                     mapping = aes(x = time, y = median, colour = simulation),
                                     size = 0.8)
      
      if(input$show_iiv){
        sim_plot <- sim_plot + geom_ribbon(data = sim_df_iiv_summary, 
                                           mapping = aes(x = time, fill = simulation,
                                                         ymin = q05, ymax = q95), alpha = 0.4, show.legend = FALSE)
        
      }
    }
    
    sim_df <- bind_rows(sim_df_noiiv, sim_df_iiv)

    xlim <- v$x
    ylim <- v$y

    if (input$x_scale == "Linear" & is.null(xlim)) {
      xlim <- range(sim_df$time)
    } else if (input$x_scale == "Logarithmic" & !is.null(v$x)){
      x_values <- sim_df$time

      shiny::validate(need(
        !all(x_values <= 0),
        "No positive value on x-axis: impossible to setup logarithmic scale."
      ))

      if (v$x[1] <= 0)
        xlim[1] <- min(x_values[x_values > 0])
    }

    if (input$y_scale == "Logarithmic" & !is.null(v$y))
    {
      y_values = sim_df[[selected_output]]

      shiny::validate(need(
        !all(y_values <= 0),
        sprintf(
          "%s is always negative along the selected time range: impossible to setup logarithmic scale.",
          selected_output
        )
      ))

      if (v$y[1] <= 0)
        ylim[1] <- min(y_values[y_values > 0])
    }

    output_reference_lines<- v$reference_lines %>% filter(output == selected_output)

    g <- sim_plot

    if(nrow(output_reference_lines) > 0){
      g <- g+
        geom_hline(data = output_reference_lines, aes(yintercept = value), linetype = "dashed")+
        geom_text(data = output_reference_lines, aes(x = Inf, y = value, label = label), inherit.aes = F, hjust = 1, vjust = -0.4)
    }

    g <- g+
      coord_cartesian(xlim = xlim, ylim = ylim)+
      xlab(label = sprintf("Time (%s)", input$design_time_unit))+
      ylab(label = sprintf("%s (%s)", annotated_output$descr, annotated_output$unit))+
      labs(color = "Simulation")+
      theme_bw()

    if (input$x_scale == "Logarithmic")
      g <- g + ggplot2::scale_x_log10()
    
    if (input$y_scale == "Logarithmic"){
      g <- g + ggplot2::scale_y_log10()
    } else if(annotated_output$unit == "%") {
      g <- g + scale_y_continuous(labels = scales::percent)
    }

    if (!is.null(selected_history))
    {
      n_chars <- nchar(c(selected_history$simulation, sim$simulation))
      cols_reduction <- max(n_chars %/% 30) + 1
      g <- g + guides(colour = guide_legend(ncol = max(1, 5 - cols_reduction)))

    }

    g <- g + theme(legend.position = "bottom", 
                   text = element_text(size=18))
    }, message = "Plotting...")

    g
  })

  # Simulated concentrations table
  output$outputs_table <- renderDataTable({
    out_table <- shinyreq(app_outputs_table())
    
    output_names <- colnames(out_table)[-c(1:3)]

    output_headers <- output_names %>%
      map_chr(function(x){
        out <- v$selected_study$outputs %>% filter(name == x)

        ifelse(is.na(out$unit),
               out$descr,
               sprintf("%s (%s)", out$descr, out$unit))

      })

    out_table <- out_table %>%
      rename(!!set_names(c(output_names, "time"), c(output_headers, sprintf("Time (%s)", input$design_time_unit))))

              datatable(
                out_table,
                options = list(dom = "rtip", scrollX = T),
                rownames = F
              )
  })

    # Individual exposure indices table
    output$individual_exposure_table <- renderUI({

      expo <- shinyreq(app_expo_table())
      n_outputs <- nrow(v$selected_study$outputs)
      
      progress <- shiny::Progress$new(session)
      progress$set(message = "Computing individual exposure", value = 0)
      on.exit(progress$close())
      
        ui <-
          purrr::map(seq_len(nrow(v$selected_study$outputs)), function(i) {
            output_name <- v$selected_study$outputs$name[i]
            output_desc <- v$selected_study$outputs$descr[i]
            output_unit <- v$selected_study$outputs$unit[i]
            
            progress$inc(amount = 1 / n_outputs, detail = output_desc)
            
            colnames_df <- tribble(~parameter, ~label,
                                   "AUC", sprintf("AUC (%s.%s)", output_unit, input$design_time_unit),
                                   "Cmax", sprintf("Cmax (%s)", output_unit),
                                   "Tmax", sprintf("Tmax (%s)", input$design_time_unit),
                                   "Cmin", sprintf("Cmin (%s)", output_unit),
                                   "Tmin", sprintf("Tmin (%s)", input$design_time_unit))
            
            output_expo <- expo %>%
              filter(output == output_name) %>%
              select(-output) %>%
              rename(!!set_names(colnames_df$parameter, nm = colnames_df$label))
            
            rh <- rhandsontable::rhandsontable(output_expo)
            list(h5(strong(sprintf("%s (%s)", output_desc, output_unit))),
                 rhandsontable::renderRHandsontable({ rh }))
          })

      ui

    })


  output$download_outputs <- downloadHandler(
    filename = function() {
      middle <-
        ifelse(nrow(v$history) == 1,
               last(v$history)[[1]]$simulation,
               "Simulations")

      sprintf("%s_%s_outputs.csv",
              v$selected_study$model@model,
              middle)
    },
    content = function(file) {
      out_table = app_outputs_table()

      write.csv(out_table, file, row.names = F)
    }
  )

  output$download_exposure_parameters <- downloadHandler(
    filename = function() {
      middle <-
        ifelse(nrow(v$history) == 1,
               last(v$history)[[1]]$simulation,
               "Simulations")

      sprintf("%s_%s_Exposure.csv",
              v$selected_study$model@model,
              middle)

    },
    content = function(file) {
      expo_table <- app_expo_table()

      write.csv(expo_table, file, row.names = F)
    }
  )


  output$computation_message <- renderText({
    sprintf(
      "Last simulation completed in %s seconds: %s computed values per output.",
      round(v$computation_duration, 4),
      v$computation_times
    )
  })

})
