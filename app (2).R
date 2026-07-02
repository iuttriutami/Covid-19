# ==============================================================================
#  COVID-19 CLUSTER ANALYSIS — SHINY APP — INDONESIAN PROVINCES
#  Interactive version of covid19_cluster_analysis_COMPLETE_EN_v3.R
#
#  Author note: converted from a single-run R script into a Shiny GUI so it
#  can be explored interactively and deployed (shinyapps.io / GitHub + Docker /
#  Posit Connect) instead of only run from the R console.
# ==============================================================================

## ---- 0. PACKAGES -----------------------------------------------------------
required_pkgs <- c(
  "shiny", "shinythemes", "DT", "ggplot2", "dplyr", "tidyr",
  "cluster", "clusterSim", "fpc", "factoextra", "mclust", "e1071",
  "dendextend", "ggrepel", "patchwork", "knitr", "RColorBrewer", "sf", "shinycssloaders"
)
missing_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_pkgs) > 0) {
  message("Missing packages detected: ", paste(missing_pkgs, collapse = ", "),
          "\nInstall them with:\n  install.packages(c(",
          paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))")
}
install.packages(c("shinythemes", "shinycssloaders"))
suppressWarnings(suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(DT)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cluster)
  library(factoextra)
  library(mclust)
  library(e1071)
  library(dendextend)
  library(ggrepel)
  library(patchwork)
  library(RColorBrewer)
  library(shinycssloaders)
}))
# Optional / heavier packages: loaded defensively so the app still runs
# (with reduced features) if they are not installed on the deploy target.
has_clusterSim <- requireNamespace("clusterSim", quietly = TRUE)
has_fpc        <- requireNamespace("fpc",        quietly = TRUE)
has_sf         <- requireNamespace("sf",         quietly = TRUE)
has_knitr      <- requireNamespace("knitr",      quietly = TRUE)

ARI <- mclust::adjustedRandIndex
set.seed(123)

## ---- 1. EMBEDDED DATA (34 PROVINCES) ---------------------------------------
dat <- data.frame(
  province = c("DKI Jakarta","Jawa Barat","Jawa Tengah","Jawa Timur","DI Yogyakarta",
    "Banten","Bali","Sumatera Utara","Sumatera Selatan","Sumatera Barat","Riau",
    "Kepulauan Riau","Jambi","Bengkulu","Lampung","Bangka Belitung","Kalimantan Barat",
    "Kalimantan Tengah","Kalimantan Selatan","Kalimantan Timur","Kalimantan Utara",
    "Sulawesi Utara","Sulawesi Tengah","Sulawesi Selatan","Sulawesi Tenggara",
    "Sulawesi Barat","Gorontalo","Maluku","Maluku Utara","Papua","Papua Barat",
    "Nusa Tenggara Barat","Nusa Tenggara Timur","Aceh"),
  code = c("31","32","33","35","34","36","51","12","16","13","14","21","15","17","18",
    "19","61","62","63","64","65","71","72","73","74","76","75","81","82","94","91",
    "52","53","11"),
  cases = c(14847.2,2127.1,1664.6,1526.5,6570.3,2234.1,3774.8,675.7,850.3,1319.7,1220.0,
    2110.0,789.0,795.8,577.3,1374.0,831.1,1048.7,1054.4,2416.8,853.3,1075.4,828.0,1168.5,
    495.2,598.9,811.2,595.0,764.3,1000.0,986.9,676.7,582.0,549.7),
  cfr = c(0.99,1.34,2.83,3.08,2.32,1.07,2.58,3.20,3.89,2.33,1.67,1.61,2.29,1.94,2.69,1.98,
    2.33,2.71,3.64,1.76,1.69,2.46,1.96,1.98,2.15,2.47,1.89,2.59,2.14,1.33,1.91,2.42,2.71,3.17),
  recovery = c(98.77,98.25,96.86,96.62,97.17,98.50,96.93,96.00,95.83,97.26,97.44,97.73,96.43,
    93.75,96.15,95.00,97.78,96.43,95.45,97.80,92.31,95.83,96.00,97.17,92.31,94.12,94.74,90.91,
    91.84,97.67,93.75,97.22,96.77,96.55),
  dose1 = c(120.5,93.4,88.2,86.7,138.0,87.5,145.3,70.2,65.4,72.1,68.9,85.6,60.3,63.8,58.4,82.3,
    64.2,73.5,67.8,102.4,78.9,78.4,58.6,74.3,59.7,52.4,60.8,57.3,54.1,41.2,44.8,62.7,46.3,50.2),
  dose2 = c(100.2,78.6,73.1,71.4,116.8,72.0,123.0,55.8,51.2,57.3,54.5,69.8,46.7,49.2,44.1,65.4,
    50.1,58.2,53.1,84.7,61.3,63.2,44.3,59.6,45.8,39.1,46.5,42.6,40.2,28.7,31.5,48.3,33.4,36.8),
  stringsAsFactors = FALSE)
rownames(dat) <- dat$province

## ---- 2. CORE FUNCTIONS (same logic as the original script) ----------------
build_features <- function(dose = c("dose2","dose1"),
                            transform = c("log","none"),
                            scaling = c("zscore","minmax","robust")) {
  dose <- match.arg(dose); transform <- match.arg(transform); scaling <- match.arg(scaling)
  cases <- if (transform == "log") log(dat$cases) else dat$cases
  vacc  <- if (dose == "dose2") dat$dose2 else dat$dose1
  M <- cbind(cases = cases, cfr = dat$cfr, recovery = dat$recovery, vacc = vacc)
  Z <- switch(scaling,
    zscore = scale(M),
    minmax = apply(M, 2, function(x) (x - min(x)) / (max(x) - min(x))),
    robust = apply(M, 2, function(x) { s <- IQR(x); (x - median(x)) / ifelse(s == 0, 1, s) }))
  Z <- as.matrix(Z); rownames(Z) <- dat$province; Z
}

cluster_labels <- function(Z, k, method = c("kmeans","ward","gmm","fcm"),
                            dist_method = "euclidean", m = 2) {
  method <- match.arg(method)
  switch(method,
    kmeans = { set.seed(123); kmeans(Z, centers = k, nstart = 50, iter.max = 100)$cluster },
    ward   = cutree(hclust(dist(Z, method = dist_method), method = "ward.D2"), k = k),
    gmm    = suppressMessages(Mclust(Z, G = k, verbose = FALSE)$classification),
    fcm    = { set.seed(123); e1071::cmeans(Z, centers = k, m = m, iter.max = 200)$cluster })
}

validate <- function(Z, cl) {
  cl <- as.integer(cl); dd <- dist(Z)
  sil <- round(mean(silhouette(cl, dd)[, 3]), 3)
  db  <- if (has_clusterSim) round(clusterSim::index.DB(Z, cl, centrotypes = "centroids")$DB, 3) else NA
  if (has_fpc) {
    cs <- fpc::cluster.stats(dd, cl)
    ch <- round(cs$ch, 1); dunn <- round(cs$dunn, 3)
  } else { ch <- NA; dunn <- NA }
  c(Silhouette = sil, Davies_Bouldin = db, Calinski_Harabasz = ch, Dunn = dunn)
}

jaccard_sim <- function(c1, c2) {
  agree <- sum(outer(c1, c1, "==") & outer(c2, c2, "=="))
  total <- sum(outer(c1, c1, "==") | outer(c2, c2, "=="))
  agree / total
}

boot_stability <- function(Z, k, FUN, n_boot = 200) {
  n <- nrow(Z); s <- numeric(n_boot); orig <- FUN(Z, k)
  for (b in seq_len(n_boot)) {
    idx <- sample(n, n, replace = TRUE)
    tryCatch({
      lab_b <- FUN(Z[idx, , drop = FALSE], k)
      s[b]  <- jaccard_sim(orig[idx], lab_b)
    }, error = function(e) s[b] <<- NA)
  }
  s <- s[!is.na(s)]
  list(mean = mean(s), sd = sd(s), ci95 = quantile(s, c(.025, .975)), scores = s)
}
fn_km  <- function(Z, k) { set.seed(sample(1e6, 1)); kmeans(Z, k, nstart = 10, iter.max = 200)$cluster }
fn_hc  <- function(Z, k) cutree(hclust(dist(Z), "ward.D2"), k = k)
fn_gmm <- function(Z, k) suppressMessages(Mclust(Z, G = k, verbose = FALSE)$classification)
fn_fcm <- function(Z, k) { set.seed(sample(1e6, 1)); e1071::cmeans(Z, k, m = 2, iter.max = 200)$cluster }

interp_stability <- function(x) ifelse(x >= .85, "Very stable",
                              ifelse(x >= .75, "Stable",
                              ifelse(x >= .60, "Moderately stable", "Unstable")))

pick_best_method <- function(tbl_validation) {
  PRIORITY <- c("K-Means","Ward","GMM","FCM")
  votes <- tbl_validation$Method[c(
    which.max(tbl_validation$Silhouette),
    which.min(tbl_validation$Davies_Bouldin),
    which.max(tbl_validation$Calinski_Harabasz),
    which.max(tbl_validation$Dunn))]
  vote_tab <- table(factor(votes, levels = PRIORITY))
  names(vote_tab)[order(-as.integer(vote_tab), match(names(vote_tab), PRIORITY))][1]
}

## ==============================================================================
##  UI
## ==============================================================================
ui <- navbarPage(
  title = div(icon("virus-covid"), "COVID-19 Cluster Analysis — Indonesian Provinces"),
  theme = shinytheme("flatly"),
  collapsible = TRUE,

  ## -- TAB: Overview & Data -------------------------------------------------
  tabPanel("Data",
    icon = icon("table"),
    sidebarLayout(
      sidebarPanel(width = 3,
        h4("About"),
        p("Cluster analysis of 34 Indonesian provinces using 4 COVID-19 indicators",
          "(cases/100k, CFR, recovery rate, vaccination dose-2 coverage),",
          "cumulative through Dec 2022."),
        tags$hr(),
        p(tags$b("Method reference (manuscript, IJDNS):")),
        tags$ul(
          tags$li("34 provinces, 4 indicators"),
          tags$li("ln(cases/100k) + z-score standardization"),
          tags$li("Main model: K-Means, k = 3"),
          tags$li("Comparators: Ward, GMM, Fuzzy C-Means")
        ),
        tags$hr(),
        downloadButton("dl_raw", "Download raw data (CSV)")
      ),
      mainPanel(width = 9,
        h4("Raw indicator data (34 provinces)"),
        DTOutput("tbl_raw") %>% withSpinner(color = "#2c3e50")
      )
    )
  ),

  ## -- TAB: Preprocessing & K selection --------------------------------------
  tabPanel("Preprocessing & K",
    icon = icon("sliders"),
    sidebarLayout(
      sidebarPanel(width = 3,
        h4("Feature engineering"),
        selectInput("dose", "Vaccination indicator", choices = c("Dose 2" = "dose2", "Dose 1" = "dose1")),
        selectInput("transform", "Cases transform", choices = c("log(cases)" = "log", "none" = "none")),
        selectInput("scaling", "Standardization", choices = c("z-score" = "zscore", "min-max" = "minmax", "robust (IQR)" = "robust")),
        tags$hr(),
        sliderInput("kmax", "Max k to test", min = 4, max = 10, value = 8, step = 1),
        actionButton("run_kselect", "Run k-selection", icon = icon("play"), class = "btn-primary")
      ),
      mainPanel(width = 9,
        h4("Standardized feature matrix (first rows)"),
        DTOutput("tbl_features") %>% withSpinner(),
        tags$hr(),
        h4("Optimal number of clusters (Elbow / Silhouette / Gap)"),
        plotOutput("plot_kselect", height = "420px") %>% withSpinner()
      )
    )
  ),

  ## -- TAB: Clustering & Validation -------------------------------------------
  tabPanel("Clustering & Validation",
    icon = icon("diagram-project"),
    sidebarLayout(
      sidebarPanel(width = 3,
        h4("Run clustering"),
        sliderInput("k_final", "Number of clusters (k)", min = 2, max = 6, value = 3),
        actionButton("run_cluster", "Run 4 methods", icon = icon("play"), class = "btn-primary"),
        tags$hr(),
        uiOutput("best_method_box")
      ),
      mainPanel(width = 9,
        h4("Cluster assignment per province"),
        DTOutput("tbl_clusters") %>% withSpinner(),
        tags$hr(),
        h4("Internal validation indices"),
        DTOutput("tbl_validation") %>% withSpinner(),
        tags$hr(),
        h4("Cluster profile (original-scale means, K-Means)"),
        DTOutput("tbl_profile") %>% withSpinner()
      )
    )
  ),

  ## -- TAB: Bootstrap Stability -----------------------------------------------
  tabPanel("Bootstrap Stability",
    icon = icon("repeat"),
    sidebarLayout(
      sidebarPanel(width = 3,
        h4("Bootstrap settings"),
        sliderInput("n_boot", "Number of bootstrap resamples", min = 50, max = 1000, value = 200, step = 50),
        helpText("Higher values are slower but more precise. The manuscript uses 1000."),
        actionButton("run_boot", "Run bootstrap stability", icon = icon("play"), class = "btn-primary")
      ),
      mainPanel(width = 9,
        h4("Jaccard stability summary"),
        DTOutput("tbl_boot") %>% withSpinner(),
        tags$hr(),
        h4("Bootstrap Jaccard distribution"),
        plotOutput("plot_boot", height = "400px") %>% withSpinner()
      )
    )
  ),

  ## -- TAB: Sensitivity Analysis -----------------------------------------------
  tabPanel("Sensitivity Analysis",
    icon = icon("magnifying-glass-chart"),
    sidebarLayout(
      sidebarPanel(width = 3,
        h4("Run all sensitivity checks"),
        p("A: number of k, B: preprocessing, C: vaccine dose, D: Ward distance,",
          "E: cross-method concordance, F: FCM fuzzifier m, G: leave-one-out."),
        actionButton("run_sens", "Run sensitivity analysis", icon = icon("play"), class = "btn-primary"),
        helpText("Leave-one-out (G) re-fits K-Means 34 times; may take a few seconds.")
      ),
      mainPanel(width = 9,
        tabsetPanel(
          tabPanel("A. Number of k", DTOutput("tbl_sensA") %>% withSpinner()),
          tabPanel("B. Preprocessing", DTOutput("tbl_sensB") %>% withSpinner()),
          tabPanel("C. Vaccine dose",  DTOutput("tbl_sensC") %>% withSpinner()),
          tabPanel("D. Ward distance", DTOutput("tbl_sensD") %>% withSpinner()),
          tabPanel("E. Cross-method ARI", plotOutput("plot_sensE", height = "380px") %>% withSpinner()),
          tabPanel("F. Fuzzifier m",   DTOutput("tbl_sensF") %>% withSpinner()),
          tabPanel("G. Leave-one-out",
                   plotOutput("plot_sensG", height = "350px") %>% withSpinner(),
                   DTOutput("tbl_sensG") %>% withSpinner())
        )
      )
    )
  ),

  ## -- TAB: Visualizations -----------------------------------------------------
  tabPanel("Visualizations",
    icon = icon("chart-line"),
    tabsetPanel(
      tabPanel("PCA (4 methods)", plotOutput("plot_pca", height = "650px") %>% withSpinner()),
      tabPanel("Silhouette (K-Means)", plotOutput("plot_sil", height = "500px") %>% withSpinner()),
      tabPanel("Ward Dendrogram", plotOutput("plot_dend", height = "500px") %>% withSpinner())
    )
  ),

  ## -- TAB: Cluster Map ----------------------------------------------------------
  tabPanel("Cluster Map",
    icon = icon("map"),
    sidebarLayout(
      sidebarPanel(width = 3,
        h4("Upload shapefile (optional)"),
        p("Upload all shapefile parts together: .shp, .shx, .dbf, .prj",
          "(e.g. indo_by_prov_2023.*). Must contain a PROVNO column matching",
          "the BPS province codes used in this app."),
        fileInput("shp_files", "Select shapefile parts", multiple = TRUE,
                   accept = c(".shp", ".shx", ".dbf", ".prj", ".cpg")),
        actionButton("run_map", "Render map", icon = icon("play"), class = "btn-primary"),
        helpText(if (!has_sf) "Package 'sf' is not installed on this server — map tab is disabled." else "")
      ),
      mainPanel(width = 9,
        plotOutput("plot_map", height = "600px") %>% withSpinner(),
        verbatimTextOutput("map_msg")
      )
    )
  ),

  ## -- TAB: Export -----------------------------------------------------------
  tabPanel("Export",
    icon = icon("download"),
    fluidPage(
      h4("Download results"),
      p("Run the relevant tabs first (Clustering, Bootstrap, Sensitivity) so results exist, then download below."),
      downloadButton("dl_results", "Download complete cluster results (CSV)"),
      br(), br(),
      downloadButton("dl_validation", "Download validation indices (CSV)"),
      br(), br(),
      downloadButton("dl_boot", "Download bootstrap stability table (CSV)")
    )
  )
)

## ==============================================================================
##  SERVER
## ==============================================================================
server <- function(input, output, session) {

  ## ---- Data tab ----
  output$tbl_raw <- renderDT({
    datatable(dat, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })
  output$dl_raw <- downloadHandler(
    filename = function() "covid19_raw_data_34provinces.csv",
    content  = function(file) write.csv(dat, file, row.names = FALSE)
  )

  ## ---- Preprocessing & K selection ----
  Xz <- reactive({
    build_features(input$dose, input$transform, input$scaling)
  })

  output$tbl_features <- renderDT({
    Zr <- round(Xz(), 3)
    df <- data.frame(Province = rownames(Zr), Zr, check.names = FALSE)
    datatable(df, rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE))
  })

  kselect_result <- eventReactive(input$run_kselect, {
    withProgress(message = "Running k-selection...", value = 0.2, {
      Z <- Xz(); set.seed(123)
      p_elbow <- fviz_nbclust(Z, kmeans, method = "wss", k.max = input$kmax) +
        labs(title = "Elbow (WSS)", x = "Number of Clusters (k)", y = "Total WSS") +
        theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"))
      incProgress(0.4)
      p_sil <- fviz_nbclust(Z, kmeans, method = "silhouette", k.max = input$kmax) +
        labs(title = "Silhouette", x = "Number of Clusters (k)", y = "Average Silhouette Width") +
        theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"))
      incProgress(0.6)
      gap_stat <- clusGap(Z, FUN = kmeans, nstart = 25, K.max = input$kmax, B = 50, verbose = FALSE)
      p_gap <- fviz_gap_stat(gap_stat) +
        labs(title = "Gap Statistic", x = "Number of Clusters (k)") +
        theme_minimal(base_size = 11) + theme(plot.title = element_text(face = "bold"))
      incProgress(1)
      p_elbow + p_sil + p_gap
    })
  }, ignoreNULL = FALSE)

  output$plot_kselect <- renderPlot({
    req(kselect_result())
    kselect_result()
  })

  ## ---- Clustering & Validation ----
  cluster_result <- eventReactive(input$run_cluster, {
    withProgress(message = "Running 4 clustering methods...", value = 0.1, {
      Z <- Xz(); k <- input$k_final
      cl_km  <- cluster_labels(Z, k, "kmeans");  incProgress(0.25)
      cl_hc  <- cluster_labels(Z, k, "ward");    incProgress(0.5)
      cl_gmm <- tryCatch(cluster_labels(Z, k, "gmm"), error = function(e) rep(NA, nrow(Z))); incProgress(0.75)
      cl_fcm <- cluster_labels(Z, k, "fcm");     incProgress(1)

      tbl_validation <- data.frame(
        Method = c("K-Means","Ward","GMM","FCM"),
        rbind(validate(Z, cl_km), validate(Z, cl_hc),
              if (all(!is.na(cl_gmm))) validate(Z, cl_gmm) else c(NA,NA,NA,NA),
              validate(Z, cl_fcm)),
        row.names = NULL)

      best_method <- tryCatch(pick_best_method(tbl_validation), error = function(e) "K-Means")

      list(Z = Z, k = k, cl_km = cl_km, cl_hc = cl_hc, cl_gmm = cl_gmm, cl_fcm = cl_fcm,
           tbl_validation = tbl_validation, best_method = best_method)
    })
  }, ignoreNULL = FALSE)

  output$best_method_box <- renderUI({
    req(cluster_result())
    tags$div(class = "well",
      h5("Best method (majority vote of 4 indices):"),
      h3(tags$b(cluster_result()$best_method), style = "color:#2c3e50;")
    )
  })

  output$tbl_clusters <- renderDT({
    req(cluster_result())
    r <- cluster_result()
    results <- data.frame(
      Province = dat$province, Code = dat$code,
      Cluster_KMeans = r$cl_km, Cluster_Ward = r$cl_hc,
      Cluster_GMM = r$cl_gmm, Cluster_FCM = r$cl_fcm,
      Cases_per100k = round(dat$cases, 1), CFR = dat$cfr,
      Recovery = dat$recovery, Vacc_Dose2 = dat$dose2)
    datatable(results, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$tbl_validation <- renderDT({
    req(cluster_result())
    datatable(cluster_result()$tbl_validation, rownames = FALSE,
              options = list(dom = "t", pageLength = 4)) %>%
      formatStyle("Method", target = "row",
                  backgroundColor = styleEqual(cluster_result()$best_method, "#d4edda"))
  })

  output$tbl_profile <- renderDT({
    req(cluster_result())
    r <- cluster_result()
    df_profile <- dat
    df_profile$Cluster_KM <- paste0("K", r$cl_km)
    prof <- df_profile %>%
      group_by(Cluster_KM) %>%
      summarise(n = n(), Cases100k = mean(cases), CFR = mean(cfr),
                Recovery = mean(recovery), Vacc2 = mean(dose2), .groups = "drop") %>%
      arrange(Cluster_KM)
    datatable(prof, rownames = FALSE, options = list(dom = "t")) %>%
      formatRound(c("Cases100k","CFR","Recovery","Vacc2"), 2)
  })

  ## ---- Bootstrap stability ----
  boot_result <- eventReactive(input$run_boot, {
    req(cluster_result())
    withProgress(message = "Bootstrapping cluster stability...", value = 0, {
      Z <- cluster_result()$Z; k <- cluster_result()$k; nb <- input$n_boot
      set.seed(123)
      incProgress(0.1); stab_km  <- boot_stability(Z, k, fn_km,  nb)
      incProgress(0.35); stab_hc  <- boot_stability(Z, k, fn_hc,  nb)
      incProgress(0.6); stab_gmm <- tryCatch(boot_stability(Z, k, fn_gmm, nb),
                                              error = function(e) list(mean=NA,sd=NA,ci95=c(NA,NA),scores=numeric(0)))
      incProgress(0.85); stab_fcm <- boot_stability(Z, k, fn_fcm, nb)
      incProgress(1)

      tbl_boot <- data.frame(
        Method = c("K-Means","Ward","GMM","FCM"),
        Mean_Jaccard = round(c(stab_km$mean, stab_hc$mean, stab_gmm$mean, stab_fcm$mean), 3),
        SD = round(c(stab_km$sd, stab_hc$sd, stab_gmm$sd, stab_fcm$sd), 3),
        CI95_Low = round(c(stab_km$ci95[1], stab_hc$ci95[1], stab_gmm$ci95[1], stab_fcm$ci95[1]), 3),
        CI95_Up  = round(c(stab_km$ci95[2], stab_hc$ci95[2], stab_gmm$ci95[2], stab_fcm$ci95[2]), 3))
      tbl_boot$Interpretation <- interp_stability(tbl_boot$Mean_Jaccard)

      boot_df <- data.frame(
        score = c(stab_km$scores, stab_hc$scores, stab_gmm$scores, stab_fcm$scores),
        method = rep(c("K-Means","Ward","GMM","FCM"),
                     times = c(length(stab_km$scores), length(stab_hc$scores),
                               length(stab_gmm$scores), length(stab_fcm$scores))))
      list(tbl_boot = tbl_boot, boot_df = boot_df)
    })
  })

  output$tbl_boot <- renderDT({
    req(boot_result())
    datatable(boot_result()$tbl_boot, rownames = FALSE, options = list(dom = "t"))
  })

  output$plot_boot <- renderPlot({
    req(boot_result())
    ggplot(boot_result()$boot_df, aes(score, fill = method, color = method)) +
      geom_density(alpha = 0.22, linewidth = 0.7) +
      geom_vline(xintercept = c(0.60, 0.75, 0.85), linetype = "dotted", color = "gray40", linewidth = 0.4) +
      scale_fill_brewer(palette = "Set1") + scale_color_brewer(palette = "Set1") +
      labs(title = paste0("Bootstrap Jaccard Distribution (n = ", input$n_boot, ")"),
           x = "Jaccard similarity", y = "Density", fill = "Method", color = "Method") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
  })

  ## ---- Sensitivity analysis ----
  sens_result <- eventReactive(input$run_sens, {
    req(cluster_result())
    withProgress(message = "Running sensitivity analysis...", value = 0, {
      Xz_ <- cluster_result()$Z; K0 <- cluster_result()$k; base_km <- cluster_result()$cl_km

      incProgress(0.1, detail = "A: number of k")
      sensA <- do.call(rbind, lapply(c("kmeans","ward","gmm","fcm"), function(meth)
        do.call(rbind, lapply(2:6, function(kk) {
          cl <- tryCatch(cluster_labels(Xz_, kk, meth), error = function(e) rep(NA, nrow(Xz_)))
          if (all(is.na(cl))) return(NULL)
          data.frame(Method = meth, k = kk, t(validate(Xz_, cl)))
        }))))

      incProgress(0.3, detail = "B: preprocessing")
      grid <- expand.grid(transform = c("log","none"), scaling = c("zscore","minmax","robust"),
                           stringsAsFactors = FALSE)
      sensB <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
        Z <- build_features(input$dose, grid$transform[i], grid$scaling[i])
        cl <- cluster_labels(Z, K0, "kmeans"); v <- validate(Z, cl)
        data.frame(transform = grid$transform[i], scaling = grid$scaling[i],
                   ARI_vs_baseline = round(ARI(base_km, cl), 3),
                   Silhouette = v["Silhouette"], Davies_Bouldin = v["Davies_Bouldin"], row.names = NULL)
      }))

      incProgress(0.45, detail = "C: vaccine dose")
      Z_d1 <- build_features("dose1", input$transform, input$scaling)
      sensC <- do.call(rbind, lapply(c("kmeans","ward","gmm","fcm"), function(meth) {
        a <- tryCatch(cluster_labels(Xz_, K0, meth), error = function(e) NULL)
        b <- tryCatch(cluster_labels(Z_d1, K0, meth), error = function(e) NULL)
        if (is.null(a) || is.null(b)) return(NULL)
        data.frame(Method = meth, ARI_dose1_vs_dose2 = round(ARI(a, b), 3))
      }))

      incProgress(0.55, detail = "D: Ward distance")
      ward_eucl <- cluster_labels(Xz_, K0, "ward", dist_method = "euclidean")
      sensD <- do.call(rbind, lapply(c("euclidean","manhattan","maximum"), function(dm)
        data.frame(Distance = dm,
                   ARI_vs_euclidean = round(ARI(ward_eucl, cluster_labels(Xz_, K0, "ward", dist_method = dm)), 3))))

      incProgress(0.65, detail = "E: cross-method concordance")
      labs_ <- list(KMeans = cluster_result()$cl_km, Ward = cluster_result()$cl_hc,
                     GMM = cluster_result()$cl_gmm, FCM = cluster_result()$cl_fcm)
      sensE <- outer(seq_along(labs_), seq_along(labs_),
                      Vectorize(function(i, j) round(ARI(labs_[[i]], labs_[[j]]), 3)))
      dimnames(sensE) <- list(names(labs_), names(labs_))

      incProgress(0.75, detail = "F: fuzzifier m")
      sensF <- do.call(rbind, lapply(c(1.5, 2.0, 2.5), function(mm)
        data.frame(m = mm, ARI_vs_KMeans = round(ARI(base_km, cluster_labels(Xz_, K0, "fcm", m = mm)), 3))))

      incProgress(0.85, detail = "G: leave-one-out")
      loo <- sapply(seq_len(nrow(dat)), function(i) {
        M <- cbind(log(dat$cases[-i]), dat$cfr[-i], dat$recovery[-i], dat$dose2[-i])
        Z <- scale(M); set.seed(123)
        cl <- kmeans(Z, centers = K0, nstart = 50, iter.max = 100)$cluster
        ARI(base_km[-i], cl)
      })
      names(loo) <- dat$province
      incProgress(1)

      list(sensA = sensA, sensB = sensB, sensC = sensC, sensD = sensD,
           sensE = sensE, sensF = sensF, loo = loo)
    })
  })

  output$tbl_sensA <- renderDT({ req(sens_result()); datatable(sens_result()$sensA, rownames = FALSE, options = list(pageLength = 10)) })
  output$tbl_sensB <- renderDT({ req(sens_result()); datatable(sens_result()$sensB, rownames = FALSE, options = list(dom = "t")) })
  output$tbl_sensC <- renderDT({ req(sens_result()); datatable(sens_result()$sensC, rownames = FALSE, options = list(dom = "t")) })
  output$tbl_sensD <- renderDT({ req(sens_result()); datatable(sens_result()$sensD, rownames = FALSE, options = list(dom = "t")) })
  output$tbl_sensF <- renderDT({ req(sens_result()); datatable(sens_result()$sensF, rownames = FALSE, options = list(dom = "t")) })

  output$plot_sensE <- renderPlot({
    req(sens_result())
    M <- sens_result()$sensE; nr <- nrow(M); nc <- ncol(M)
    par(mar = c(4, 5, 3.5, 1)); cols <- colorRampPalette(c("#fff5eb", "#fd8d3c", "#7f2704"))(100)
    image(1:nc, 1:nr, t(M[nr:1, , drop = FALSE]), col = cols, zlim = c(min(M, na.rm = TRUE), 1),
          axes = FALSE, xlab = "", ylab = "",
          main = "Adjusted Rand Index between methods", cex.main = 1.1)
    axis(1, at = 1:nc, labels = colnames(M), tick = FALSE, cex.axis = 0.95)
    axis(2, at = 1:nr, labels = rev(rownames(M)), las = 1, tick = FALSE, cex.axis = 0.95)
    for (i in 1:nr) for (j in 1:nc)
      text(j, nr - i + 1, formatC(M[i, j], format = "f", digits = 3),
           col = ifelse(M[i, j] > 0.6, "white", "black"), cex = 1)
    box()
  })

  output$plot_sensG <- renderPlot({
    req(sens_result())
    loo <- sens_result()$loo
    hist(loo, breaks = seq(min(loo) - 0.025, 1.025, by = 0.025),
         col = "#9ecae1", border = "white", main = "Leave-one-province-out stability",
         xlab = "Adjusted Rand Index vs baseline", ylab = "Frequency (provinces)")
    abline(v = mean(loo), col = "#e31a1c", lwd = 2, lty = 2)
    legend("topleft", legend = sprintf("Mean = %.3f", mean(loo)), col = "#e31a1c", lwd = 2, lty = 2, bty = "n")
  })

  output$tbl_sensG <- renderDT({
    req(sens_result())
    loo <- sens_result()$loo
    df <- data.frame(Province = names(loo), ARI_LOO = round(loo, 4)) %>% arrange(ARI_LOO)
    datatable(df, rownames = FALSE, options = list(pageLength = 8))
  })

  ## ---- Visualizations ----
  output$plot_pca <- renderPlot({
    req(cluster_result())
    r <- cluster_result(); Z <- r$Z
    pca <- prcomp(Z, scale. = FALSE)
    pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], province = rownames(Z))
    ve <- round(summary(pca)$importance[2, 1:2] * 100, 1)
    make_pca <- function(labels, title) {
      df <- pca_df; df$cluster <- factor(labels)
      ggplot(df, aes(PC1, PC2, color = cluster, label = province)) +
        geom_point(size = 2.6, alpha = 0.85) +
        geom_text_repel(size = 2.4, max.overlaps = 20, segment.size = 0.3) +
        stat_ellipse(aes(fill = cluster), geom = "polygon", alpha = 0.08, level = 0.90, type = "t") +
        scale_color_brewer(palette = "Dark2", name = "Cluster") +
        scale_fill_brewer(palette = "Dark2", guide = "none") +
        labs(title = title, x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)")) +
        theme_minimal(base_size = 10) +
        theme(plot.title = element_text(face = "bold", size = 11), legend.position = "bottom")
    }
    (make_pca(r$cl_km, "K-Means") | make_pca(r$cl_hc, "Ward")) /
      (make_pca(r$cl_gmm, "GMM")  | make_pca(r$cl_fcm, "FCM")) +
      plot_annotation(title = "PCA Visualization — Four Methods",
                       subtitle = "34 provinces | 90% confidence ellipse",
                       theme = theme(plot.title = element_text(face = "bold", size = 14)))
  })

  output$plot_sil <- renderPlot({
    req(cluster_result())
    r <- cluster_result(); d <- dist(r$Z)
    fviz_silhouette(silhouette(r$cl_km, d), palette = "Dark2", ggtheme = theme_minimal(base_size = 11)) +
      labs(title = "Silhouette — K-Means") + theme(plot.title = element_text(face = "bold"))
  })

  output$plot_dend <- renderPlot({
    req(cluster_result())
    r <- cluster_result()
    hc_obj <- hclust(dist(r$Z, method = "euclidean"), method = "ward.D2")
    dend <- as.dendrogram(hc_obj)
    dend <- color_branches(dend, k = r$k); dend <- set(dend, "labels_cex", 0.75)
    par(mar = c(8, 4, 4, 2))
    plot(dend, main = paste0("Ward Dendrogram (k = ", r$k, ")"), ylab = "Euclidean Distance (Ward)", xlab = "")
  })

  ## ---- Cluster Map ----
  map_result <- eventReactive(input$run_map, {
    req(cluster_result())
    if (!has_sf) return(list(error = "Package 'sf' is not available on this server."))
    files <- input$shp_files
    if (is.null(files)) return(list(error = "Please upload the shapefile parts (.shp/.shx/.dbf/.prj) first."))

    tmpdir <- tempfile(); dir.create(tmpdir)
    shp_path <- NULL
    for (i in seq_len(nrow(files))) {
      dest <- file.path(tmpdir, files$name[i])
      file.copy(files$datapath[i], dest)
      if (grepl("\\.shp$", files$name[i], ignore.case = TRUE)) shp_path <- dest
    }
    if (is.null(shp_path)) return(list(error = "No .shp file found among the uploads."))

    res <- tryCatch({
      map_sf <- sf::st_read(shp_path, quiet = TRUE)
      if (!"PROVNO" %in% names(map_sf))
        return(list(error = paste("Column 'PROVNO' not found. Available columns:",
                                   paste(names(map_sf), collapse = ", "))))
      map_sf$code <- trimws(as.character(map_sf$PROVNO))

      r <- cluster_result()
      label_tbl <- data.frame(
        code = dat$code, province = dat$province,
        `K-Means` = factor(r$cl_km), Ward = factor(r$cl_hc),
        GMM = factor(r$cl_gmm), FCM = factor(r$cl_fcm),
        stringsAsFactors = FALSE, check.names = FALSE)

      best <- r$best_method
      if (!best %in% names(label_tbl)) best <- "K-Means"
      map_join <- dplyr::left_join(map_sf, label_tbl, by = "code")
      map_join$Best <- factor(label_tbl[[best]][match(map_join$code, label_tbl$code)])

      pal3 <- c("1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3",
                "4" = "#e7298a", "5" = "#66a61e", "6" = "#e6ab02")
      p <- ggplot(map_join) +
        geom_sf(aes(fill = Best), color = "white", linewidth = 0.15) +
        scale_fill_manual(values = pal3, name = "Cluster", na.value = "grey90") +
        labs(title = sprintf("COVID-19 Cluster Map of Indonesian Provinces — %s (k = %d)", best, r$k),
             subtitle = "Joined via BPS regional code (PROVNO)") +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(face = "bold"), axis.text = element_blank(), panel.grid = element_blank())
      list(plot = p, error = NULL)
    }, error = function(e) list(error = paste("Failed to read/join shapefile:", conditionMessage(e))))
    res
  })

  output$plot_map <- renderPlot({
    req(map_result())
    if (!is.null(map_result()$error)) return(NULL)
    map_result()$plot
  })
  output$map_msg <- renderText({
    req(map_result())
    if (!is.null(map_result()$error)) map_result()$error else "Map rendered successfully."
  })

  ## ---- Export ----
  output$dl_results <- downloadHandler(
    filename = function() "complete_cluster_results.csv",
    content = function(file) {
      req(cluster_result()); r <- cluster_result()
      results <- data.frame(
        No = seq_len(nrow(dat)), Province = dat$province, Code = dat$code,
        Cluster_KMeans = r$cl_km, Cluster_Ward = r$cl_hc,
        Cluster_GMM = r$cl_gmm, Cluster_FCM = r$cl_fcm,
        Cases_per100k = round(dat$cases, 1), CFR = dat$cfr,
        Recovery = dat$recovery, Vacc_Dose2 = dat$dose2)
      write.csv(results, file, row.names = FALSE)
    }
  )
  output$dl_validation <- downloadHandler(
    filename = function() "internal_validation.csv",
    content = function(file) { req(cluster_result()); write.csv(cluster_result()$tbl_validation, file, row.names = FALSE) }
  )
  output$dl_boot <- downloadHandler(
    filename = function() "bootstrap_stability.csv",
    content = function(file) { req(boot_result()); write.csv(boot_result()$tbl_boot, file, row.names = FALSE) }
  )
}

shinyApp(ui, server)
