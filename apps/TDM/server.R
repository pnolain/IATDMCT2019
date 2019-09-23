library(tidyverse)
library(lubridate)
library(mrgsolve)
library(rhandsontable)
library(jsonlite)
library(shinyjs)


shiny::shinyServer(function(input, output, session) {
  mod <-
    mread_cache("pkpd", project = "www") %>% update(end = 3600, delta = 48)
  
  theme_set(pmxploit::theme_pmx())
  
  start_date <- ymd(as.character(Sys.Date()))
  
  rv <- reactiveValues(
    admin = tibble(
      time = start_date,
      dose = 75L,
      n = 6L,
      interval = 14L
    ),
    tdm = tibble(
      time = start_date + hours(12 * 6),
      output = "LDLC",
      cc = 70
    )  %>% slice(0),
    cov = tibble(
      AGE = 60,
      SEX = factor(
        0,
        levels = c(0, 1),
        labels = c("Male", "Female")
      ),
      WT = 82.9,
      STATIN = factor(0, levels = c(0, 1), labels = c("No", "Yes")),
      HDSTATIN =  factor(0, levels = c(0, 1), labels = c("No", "Yes")),
      BSLDLC = 134
    ),
    conv = NULL,
    obf_fn_val = 0,
    
    conv_param = tibble(method = "Nelder-Mead",
                        maxit = 1000L)
    
  )
  
  
  hide(
    id = "loading-content",
    anim = TRUE,
    animType = "fade",
    time = 2.5
  )
  
  output$administration_table <- renderRHandsontable({
    rhandsontable(
      as.data.frame(rv$admin),
      useTypes = T,
      rowHeaders = F,
      stretchH = "all",
      overflow = "hidden"
    ) %>%
      rhandsontable::hot_col(col = "time", type = "date")
  })
  
  
  output$covariates_table <- renderRHandsontable({
    rhandsontable(
      as.data.frame(rv$cov),
      rowHeaders = F,
      readOnly = F,
      height = 100,
      stretchH = "all",
      overflow = "hidden"
    ) %>%
      hot_col(col = "AGE", format = 0) %>%
      hot_col(col = "BSLDLC", format = 0)
  })
  
  output$observations_table <- renderRHandsontable({
    output_choice <- c("LDLC", "Alirocumab")
    rhandsontable(
      as.data.frame(rv$tdm),
      useTypes = T,
      rowHeaders = F,
      stretchH = "all",
      readOnly = F,
      overflow = "hidden",
      height = 100
    ) %>%
      rhandsontable::hot_col(col = "time", type = "date") %>%
      hot_col(
        col = "output",
        type = "dropdown",
        source = output_choice,
        strict = TRUE
      )
  })
  
  observeEvent(input$administration_table, {
    df <- rhandsontable::hot_to_r(input$administration_table)
    rv$admin <- as_tibble(df)
  })
  
  # observeEvent(input$observations_table, {
  #   if (!is_empty(input$observations_table$data)) {
  #     df <- rhandsontable::hot_to_r(input$observations_table)
  #     rv$tdm <- as_tibble(df) %>% arrange(time)
  #   }
  #   
  # })
  
  observeEvent(input$covariates_table, {
    df <- rhandsontable::hot_to_r(input$covariates_table)
    rv$cov <- as_tibble(df)
  })
  
  observeEvent(input$add_admin, {
    last_admin_row <- rv$admin %>%
      slice(n())
    
    next_admin_time <-
      (last_admin_row$time + days(last_admin_row$interval * (last_admin_row$n)))
    
    rv$admin <- rv$admin %>%
      add_row(time = next_admin_time) %>%
      fill(everything())
  })
  
  observeEvent(input$clear_admin, {
    rv$admin <- rv$admin %>% slice(1)
    
  })
  
  
  observeEvent(input$reset, {
    rv$cov = tibble(
      AGE = 60,
      SEX = factor(
        0,
        levels = c(0, 1),
        labels = c("Male", "Female")
      ),
      WT = 82.9,
      STATIN = factor(0, levels = c(0, 1), labels = c("No", "Yes")),
      HDSTATIN =  factor(0, levels = c(0, 1), labels = c("No", "Yes")),
      BSLDLC = 134
    )
  })
  
  
  observeEvent(input$add_sample, {
    rv$tdm <- rv$tdm %>%
      add_row(time = Sys.Date(),
              output = "LDLC",
              cc = 70)
  })
  
  observeEvent(input$clear_samples, {
    rv$tdm <- rv$tdm %>%
      slice(0)
    rv$obf_fn_val = 0
    
  })
  
  observeEvent(input$dosing_plot_click, {
    val <- input$dosing_plot_click
    
    obs <- rv$tdm %>% mutate(n_row = row_number())
    obs$time <- as.numeric(obs$time)
    
    np <-
      nearPoints(
        obs,
        coordinfo = val,
        xvar = "time",
        yvar = "cc",
        panelvar1 = "output",
        threshold = 10,
        maxpoints = 1,
        addDist = FALSE,
        allRows = FALSE
      )
    
    if (nrow(np) > 0) {
      rv$tdm <- rv$tdm %>% slice(-np$n_row)
    } else {
      new_sample <-
        tibble(
          time = max(start_date, as.Date(val$x, origin = "1970-01-01")),
          output = val$panelvar1,
          cc = val$y
        )
      
      rv$tdm <- rv$tdm %>%
        bind_rows(new_sample) %>%
        arrange(time)
    }
    
  })
  
  omega_pk <- mod@omega@data$PK
  omega_pd <- mod@omega@data$PD
  
  n_om <- nrow(omega_pd) + nrow(omega_pk)
  omega_full <- matrix(data = 0,
                       nrow = n_om,
                       ncol = n_om)
  omega_full[1:nrow(omega_pd), 1:nrow(omega_pd)] <- omega_pd
  omega_full[(nrow(omega_pd) + 1):(nrow(omega_pd) + nrow(omega_pk)),
             (nrow(omega_pd) + 1):(nrow(omega_pd) + nrow(omega_pk))] <-
    omega_pk
  
  omega_full.inv <- solve(omega_full)
  sigma_ad <-
    matrix(c(0.0465442844971299882 ^ 2, 5.2105220468796993 ^ 2))
  sigma_prop <-
    matrix(c(0.25922786055125713 ^ 2, 0.22391425764936213 ^ 2))
  
  init_theta <- c(
    c(
      ETACLL = 0,
      ETAV2 = 0,
      ETAV3 = 0,
      ETAKM = 0,
      ETAF1 = 0
    ),
    c(
      ETAKOUT = 0,
      ETAEC50 = 0,
      ETAEMAX = 0,
      ETAGAM = 0
    )
  )
  
  
  mapbayes <-
    function(eta,
             y = NULL,
             d = NULL,
             m = NULL,
             pred = FALSE) {
      sig2_ad <- as.numeric(sigma_ad)
      sig2_prop <- as.numeric(sigma_prop)
      eta <- eta %>% as.list
      names(eta) <- names(init_theta)
      eta_m <- eta %>% unlist %>% matrix(nrow = 1)
      
      # m <- mod
      m <- m %>% param(eta)
      
      if (!pred) {
        m <- m %>% obsonly
      } else {
        last_admin_row <- d %>%
          filter(!is.na(amt)) %>%
          slice(n())
        last_admin_end <-
          last_admin_row$time + last_admin_row$ii * (last_admin_row$addl) + last_admin_row$ii
        m <-
          m %>% obsaug %>% update(delta = 12,
                                  end = max(last_admin_end, m@end))
        
      }
      out <- m %>% zero_re() %>% data_set(d) %>%  mrgsim
      
      if (pred)
        return(out)
      # http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3339294/
      
      if (nrow(y) == 0)
        return(0)
      
      sig2j <- out$DV ^ 2 * sig2_prop[out$dv_type] + sig2_ad[out$dv_type]
      sqwres <- log(sig2j) + (1 / sig2j) * (y$cc - out$DV) ^ 2
      
      nOn3 <- diag(eta_m %*% omega_full.inv %*% t(eta_m))
      
      return(sum(sqwres) + nOn3)
    }
  
  
  
  simulation <- reactive({
    input$go
    rv$conv = NULL
    # withBusyIndicator("go",{
    
    isolate({
      obs <- rv$tdm %>%
        mutate(
          ID = 1L,
          cmt = ifelse(output == "LDLC", 4L, 2L),
          evid = 0L,
          dv_type = as.integer(ifelse(cmt == 4L, 2L, 1L))
        ) %>%
        select(-output)
      
      ttt <- rv$admin %>%
        mutate(
          ID = 1L,
          cmt = 1L,
          evid = 1L,
          dose = dose,
          n = n - 1,
          dv_type = 999
        ) %>%
        rename(amt = dose,
               addl = n,
               ii = interval)
      
      obs_ttt <- bind_rows(ttt, obs) %>%
        arrange(time) %>%
        mutate(ii = ii * 24) %>%
        mutate(time = interval(time[1], time) / hours(1))
      
      ds <-
        obs_ttt %>% merge(rv$cov) %>% mutate_if(is.factor, ~ as.integer(.) - 1) %>% as_tibble()
      
      safe_optim <- safely(optim)
      
      fit_bay <-
        safe_optim(
          init_theta,
          mapbayes,
          y = obs,
          d =  ds,
          m = mod,
          method = "BFGS"
        )
      
      if(is.null(fit_bay$error)){
        fit_bay <- fit_bay$result
        
        rv$conv <- fit_bay$convergence
        rv$obf_fn_val <- round(fit_bay$value, 6)
      } else{
        rv$conv <- 9999999
        rv$obf_fn_val <- 9999999
        print("error")
        return(NULL)
      }
      
      list(
        mapbayes = mapbayes(
          eta = fit_bay$par,
          d = ds,
          m = mod,
          pred = TRUE
        ),
        population = mapbayes(
          eta = init_theta,
          d = ds,
          m = mod,
          pred = TRUE
        ),
        eta = fit_bay$par,
        first_admin_time = min(rv$admin$time)
      )
    })
    #  hide(id = "computing", anim = TRUE, animType = "fade", time = 2.5)
  })
  # })
  
  
  output$obj_Fn_value <- renderText({
    ifelse (rv$obf_fn_val != 0,
            paste0("Obj. Fn value: ", rv$obf_fn_val),
            "")
    
  })
  
  
  output$model_problem <- renderUI({
    if (rv$obf_fn_val != 0) {
      if (rv$conv == 0) {
        tags$span(icon("check"),
                  style = paste("color:", "green"),
                  "convergence OK")
      } else{
        tags$span(icon("exclamation-circle"),
                  style = paste("color:", "red"),
                  "convergence error")
      }
      
      
    }
  })
  
  
  
  
  
  forecast <- reactive({
    sim <- shiny::req(simulation())
    
    isolate({
      ttt <- rv$admin %>%
        mutate(
          ID = 1L,
          cmt = 1L,
          evid = 1L,
          dose = dose,
          n = n - 1
        ) %>%
        rename(amt = dose) %>%
        merge(rv$cov) %>% mutate_if(is.factor, ~ as.integer(.) - 1)
      
      last_admin_row <- ttt %>%
        slice(n())
      
      last_admin_end <-
        (last_admin_row$time + days(last_admin_row$interval * (last_admin_row$n + 1)))
      
      
      ttt2 <- ttt %>%
        add_row(time = last_admin_end)
      
      ttt3 <- ttt %>%
        add_row(time = last_admin_end,
                amt = last_admin_row$amt * 0.5)
      
      ttt4 <- ttt %>%
        add_row(time = last_admin_end,
                amt = last_admin_row$amt * 2)
      
      ttt5 <- ttt %>%
        add_row(
          time = last_admin_end,
          amt = last_admin_row$amt * 4,
          interval = last_admin_row$interval * 2
        )
      
      forecast_sims <- tribble(
        ~ name,
        ~ treatment,
        "1: Same dose",
        ttt2,
        "2: Dose / 2",
        ttt3,
        "3: Dose * 2",
        ttt4,
        "4: Dose * 4 & Interval * 2",
        ttt5
      ) %>%
        mutate(dataset = map(treatment, function(x) {
          ds <- x %>%
            fill(everything()) %>%
            arrange(time) %>%
            mutate(time = interval(time[1], time) / hours(1),
                   ii = interval * 24) %>%
            rename(addl = n)
        })) %>%
        mutate(simulation = map(dataset, ~ mapbayes(
          eta = sim$eta,
          d = .,
          m = mod,
          pred = TRUE
        )))
      
    })
    
  })
  
  
  
  output$dosing_plot <- renderPlot({
    res <- shiny::req(simulation())
    
    i_sim <- res$mapbayes
    p_sim <- res$population
    first_admin_time <- res$first_admin_time
    # first_admin_time <- min(rv$admin$time)
    
    
    
    sim_df <- bind_rows(list(
      ipred = as_tibble(i_sim@data),
      pred = as_tibble(p_sim@data)
    ),
    .id = "prediction") %>%
      mutate(time = as.Date(time / 24, origin = first_admin_time)) %>%
      select(time, Alirocumab, LDLC, prediction) %>%
      gather(output, cc, Alirocumab, LDLC)
    
    obs_df <- rv$tdm %>%
      rename(output = output)
    
    ggplot(sim_df, aes(
      x = time,
      y = cc,
      linetype = prediction
    )) +
      geom_line(size = 1.1) +
      geom_point(
        data = obs_df,
        aes(x = time, y = cc),
        inherit.aes = FALSE,
        size = 3,
        color = "red"
      ) +
      scale_x_date(date_breaks = "2 week") +
      facet_wrap( ~ output, scales = "free") +
      scale_linetype_manual(
        breaks = c("ipred", "pred"),
        labels = c("Individual", "Population"),
        values = c("solid", "dashed")
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    
  })
  
  output$parameters_table <- renderRHandsontable({
    res <- shiny::req(simulation())
    
    pop <- res$mapbayes@mod@param
    
    ip <- as_tibble(res$mapbayes@data) %>%
      group_by(ID) %>%
      slice(1) %>%
      ungroup() %>%
      select(BIO:KIN)
    
    
    param <- tribble(
      ~ parameter,
      ~ population,
      ~ individual,
      ~ eta,
      "CLL",
      pop$TVCLL,
      ip$CLL,
      res$eta["ETACLL"],
      "V2",
      pop$TVV2,
      ip$V2,
      res$eta["ETAV2"],
      "V3",
      pop$TVV3,
      ip$V3,
      res$eta["ETAV3"],
      "KM",
      pop$TVKM,
      ip$KM,
      res$eta["ETAKM"],
      "F1",
      pop$TVF1,
      ip$BIO,
      res$eta["ETAF1"],
      "KOUT",
      pop$TVKOUT,
      ip$KOUT,
      res$eta["ETAKOUT"],
      "EC50",
      pop$TVEC50,
      ip$EC50,
      res$eta["ETAEC50"],
      "EMAX",
      pop$TVEMAX,
      ip$EMAX,
      res$eta["ETAEMAX"],
      "GAM",
      pop$TVGAM,
      ip$GAM,
      res$eta["ETAGAM"],
      "KA",
      pop$TVKA,
      ip$KA,
      0,
      "Q",
      pop$TVQ,
      ip$Q,
      0,
      "VM",
      pop$TVVM,
      ip$VM,
      0,
      "LAG",
      pop$TVLAG,
      ip$LAG,
      0
    ) %>%
      gather(" ", value,-parameter) %>%
      spread(parameter, value) %>% select(" ", F1, LAG, KA, V2, CLL, VM, KM, Q, V3, KOUT, EC50, EMAX, GAM)
    
    rhandsontable(
      as.data.frame(param),
      readOnly = F,
      rowHeaders = F,
      stretchH = "all",
      overflow = "hidden"
    ) %>%
      hot_col(
        col = c("KA", "Q", "VM", "LAG"),
        renderer = "
           function (instance, td, row, col, prop, value, cellProperties) {
             Handsontable.renderers.NumericRenderer.apply(this, arguments);
              td.style.color = 'lightgrey';
           }"
      )
  })
  
  output$forecast_plot <- renderPlot({
    res <- shiny::req(forecast())
    
    isolate({
      res_df <-
        res %>% mutate(data = map(simulation, ~ as_tibble(.@data)))
      first_admin <- as.data.frame(res$treatment[1])[1, 1]
      
      last_admin_row <- rv$admin %>%
        slice(n())
      
      last_admin_end <-
        (last_admin_row$time + days(last_admin_row$interval * (last_admin_row$n)))
      sim_df <-
        bind_rows(set_names(res_df$data, res_df$name), .id = "treatment") %>%
        mutate(time = as.Date(time / 24, origin = first_admin)) %>%
        mutate(type = ifelse(time < last_admin_end, "history", "forecast"))
      
      bsldlc <- mod@param$BSLDLC
      
      ggplot(sim_df, aes(
        x = time,
        y = LDLC,
        colour = type
      )) +
        geom_line() +
        geom_hline(yintercept = bsldlc * 0.5, linetype = "dashed") +
        scale_x_date(date_breaks = "2 week") +
        scale_y_continuous(sec.axis = sec_axis(
          ~ (bsldlc - .) / bsldlc,
          name = "LDLC decrease (%)",
          labels = scales::percent
        )) +
        facet_wrap( ~ treatment, nrow = 1, scales = "free_x") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
    })
  })
  
  output$distributions_plot <- renderPlot({
    res <- shiny::req(simulation())
    
    diag_pk <-
      set_names(diag(mod@omega@data$PK),
                c("ETACLL", "ETAV2", "ETAV3", "ETAKM", "ETAF1"))
    diag_pd <-
      set_names(diag(mod@omega@data$PD),
                c("ETAKOUT", "ETAEC50", "ETAEMAX", "ETAGAM"))
    
    full_diag <- c(diag_pk, diag_pd)
    
    ip_df <- enframe(res$eta, name = "param") %>%
      mutate(param_type = ifelse(param %in% names(diag_pk), "PK", "PD"))#%>%
    # filter(param %in% c("ETACLL", "ETAV2", "ETAV3", "ETAKM", "ETAF1"))
    
    pop_df <- tibble(x = seq(-2, 2, length.out = 100)) %>%
      mutate(etas = map(x, function(x) {
        tibble(
          param = names(full_diag),
          param_type = ifelse(param %in% names(diag_pk), "PK", "PD")
        ) %>%
          mutate(value = dnorm(x, mean = 0, sd = sqrt(full_diag[param])))
      })) %>%
      unnest() 
    
    ggplot(pop_df, aes(x)) +
      geom_line(aes(y = value)) +
      geom_vline(
        data = ip_df,
        aes(xintercept = value),
        color = "red",
        linetype = "dashed"
      ) +
      facet_wrap( ~ param, nrow = 2)
  })
})