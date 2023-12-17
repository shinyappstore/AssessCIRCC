#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
#devtools::install_github("debruine/faux")
library(faux)
library(pwr)
library(psych)
library(ggplot2)
library(cowplot)
library(ggthemes)
library(RColorBrewer)

# Define UI for application that draws a histogram ----
ui <- fluidPage(
    
    # Application title
    titlePanel("Assessing Change in Intervention Research: The Benefits of Composite Outcomes"),
    
    strong("This shiny app accompanies the paper 'Assessing Change in Intervention Research: The Benefits of Composite Outcomes'. Code for the shiny app can be found at"),
    a("https://osf.io/u96em/", href = "https://osf.io/u96em/", target = "_blank"),
    ".",
    p("Intervention research is often time- and resource-intensive, with numerous participants involved over extended periods of time. In order to maximize the value of intervention studies, multiple outcome measures are often included, either to ensure a diverse set of outcomes is being assessed or to refine assessments of specific outcomes. In the paper, we advocate for combining assessments, rather than relying on individual measures assessed separately, to better evaluate the effectiveness of interventions. Specifically, we argue that by pooling information from individual measures into a single outcome, composite scores can provide finer estimates of the underlying theoretical construct of interest, while retaining important properties more sophisticated methods often forego, such as transparency and interpretability.", style = "font-size:12px"),
    p("To illustrate this, here we simulate (N = 100) a simple intervention study that includes multiple measures. The measures are either (i) treated separately in the analysis stage and a significant result in any measure is interpreted as evidence for the effectiveness of the intervention, (ii) one of the measures is specified as the primary outcome, or (iii) the measures are combined into a composite outcome. The statistical test conducted is a one-tailed t-test (only the intervention group is simulated for simplification purposes). The number of participants, number of measures, and the effect size can be determined by the user. The user can also specify the correlation between measures to get at the idea of broad vs narrow constructs (see paper for details). Finally, the user can specify whether the intervention modulates the theoretical construct, or task-specific variance that is not related to the underlying construct (see e.g. Figure 1 in the paper).", style = "font-size:12px"),
    p("Note that we make several assumptions, including that all measures are initially equally unbiased and precise (M = 0, SD = 1), that the initial correlations between measures are the same, and that the effect varies across participants (SD = 0.5) but is correlated highly across all measures (r = .8). If the intervention is set to modulate task-specific variance, only one measure shows the effect, whereas the other measures show no improvements. These assumptions are made for simplification purposes, but can easily be changed by adapting the provided code.", style = "font-size:12px"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            numericInput(
                "npeople", "Number of participants:", value = 30, min = 5, step = 1
            ),
            numericInput(
                "nmeasures", "Number of measures:", value = 3, min = 3, step = 1
            ),
            numericInput(
                "es", "Effect size (d):", min = 0, max = 1, value = 0.1, step = .1
            ),
            numericInput(
                "initial_cor", "Correlation between measures before intervention:", min = 0.1, max = 0.9, value = 0.6, step = .1
            ),
            radioButtons(
                "test_specific_effect", "Is the modulated variance related to:", choices = list("the construct", "a specific task"), selected = "the construct"
            ),
            actionButton("start_sim", "Simulate!")
        ),
        
        mainPanel(
            plotOutput("simPlot")
        )
    )
)


server <- function(input, output) {
    
    run_sim <- function () {
        N_people <- input$npeople
        N_measures <- input$nmeasures
        es <- input$es
        ifelse(input$test_specific_effect == "the construct", test_specific_effect <- FALSE, test_specific_effect <- TRUE)
        
        within <- list(measure = 1:N_measures)
        df_p_values <- data.frame()
        df_p_values_null <- data.frame()
        
        for (i_sim in 1:100) {
            
            # Simulate pretest data ----
            mycor = input$initial_cor #runif(N_measures * (N_measures - 1) / 2, .7, .9)
            df_sample <- sim_design(within, 
                                    n = N_people, mu = 0, sd = 1, r = mycor,
                                    empirical = FALSE, plot = FALSE, long = TRUE)
            
            # Calculate average
            df_sample_average <- data.frame(id = unique(df_sample$id), measure = rep("average", N_people))
            df_sample_average$y <- aggregate(df_sample$y, list(df_sample$id), mean)$x
            
            # Calculate scaled average
            df_sample_wide <- long2wide(df_sample, within = "measure")
            mean_cormat_sample <- sum(cor(df_sample_wide[,2:ncol(df_sample_wide)]))
            grand_mean_sample <- mean(df_sample$y)
            sample_sum <- aggregate(df_sample$y, list(df_sample$id), sum)
            
            df_sample_scaled_average <- data.frame(id = unique(df_sample$id), measure = rep("scaled_average", N_people))
            df_sample_scaled_average$y <- (sample_sum$x - N_measures * grand_mean_sample) / sqrt(mean_cormat_sample) + grand_mean_sample
            
            # Calculate EFA
            fa_result <- fa(df_sample_wide[,2:ncol(df_sample_wide)], nfactors = 1)
            fa_weights <- fa_result$weights
            
            df_sample_fa_coarse <- data.frame(id = unique(df_sample$id), measure = rep("fa_coarse", N_people))
            df_sample_fa_coarse$y <- rowSums(mapply(`*`, df_sample_wide[,2:ncol(df_sample_wide)], fa_result$loadings))
            
            df_sample_fa_refined <- data.frame(id = unique(df_sample$id), measure = rep("fa_refined", N_people))
            df_sample_fa_refined$y <- rowSums(mapply(`*`, df_sample_wide[,2:ncol(df_sample_wide)], fa_result$weights))
            
            
            # Combine dataframes
            df_sample_all <- rbind(df_sample, df_sample_average, df_sample_scaled_average, df_sample_fa_coarse, df_sample_fa_refined)
            
            
            # Simulate posttest data - effect present ----
            ifelse(test_specific_effect == TRUE,
                   es_vector <- rnorm_multi(n = N_people,
                                            vars = N_measures,
                                            mu = c(es, rep(0, N_measures - 1)),
                                            sd = 0.5,
                                            r = 0,
                                            varnames = 1:N_measures,
                                            empirical = FALSE),
                   es_vector <- rnorm_multi(n = N_people,
                                            vars = N_measures,
                                            mu = es,
                                            sd = 0.5,
                                            r = 0.8,
                                            varnames = 1:N_measures,
                                            empirical = FALSE))
            
            #Create correlated copy of pre-test scores to get post-test scores without treatment effect
            for (i in 1:N_measures) {
                df_sample_wide[N_measures + i + 1] <- rnorm_pre(df_sample_wide[,i+1], mu = 0, sd = 1, r = 0.8, empirical = TRUE)
            }
            
            # Add treatment effect if effect present
            for (i in 1:N_measures) {
                df_sample_wide[N_measures * 2 + i + 1] <- df_sample_wide[N_measures + i + 1] + es_vector[,i]
            }
            
            # Add to df_sample
            for (i in 1:N_measures) {
                df_sample[df_sample$measure == i,4] <- df_sample_wide[N_measures * 2 + i + 1]
            }
            names(df_sample) <- c("id", "measure", "y", "y_post")
            
            
            # Calculate average
            df_sample_average$y_post <- aggregate(df_sample$y_post, list(df_sample$id), mean)$x
            
            # Calculate scaled average
            tmp <- df_sample
            tmp$y <- tmp$y_post
            df_sample_post_wide <- long2wide(tmp, within = "measure")
            #mean_cormat_sample_post <- sum(cor(df_sample_post_wide[,2:ncol(df_sample_post_wide)])) # We are taking the correlation structure from pre-test
            grand_mean_sample_post <- mean(df_sample$y_post)
            sample_post_sum <- aggregate(df_sample$y_post, list(df_sample$id), sum)
            
            df_sample_scaled_average$y_post <- (sample_post_sum$x - N_measures * grand_mean_sample_post) / sqrt(mean_cormat_sample) + grand_mean_sample_post
            
            # Calculate EFA
            df_sample_fa_coarse$y_post <- rowSums(mapply(`*`, df_sample_post_wide[,2:ncol(df_sample_post_wide)], fa_result$loadings))
            df_sample_fa_refined$y_post <- rowSums(mapply(`*`, df_sample_post_wide[,2:ncol(df_sample_post_wide)], fa_result$weights))
            
            
            # Combine dataframes
            df_sample_all <- rbind(df_sample, df_sample_average, df_sample_scaled_average, df_sample_fa_coarse, df_sample_fa_refined)
            df_sample_all$y_diff <- df_sample_all$y_post - df_sample_all$y
            
            
            # Simulate posttest data - effect absent ----
            es_vector_null <- rnorm_multi(n = N_people,
                                          vars = N_measures,
                                          mu = 0,
                                          sd = 0.5,
                                          r = 0.8,
                                          varnames = 1:N_measures,
                                          empirical = FALSE)
            
            # Add treatment effect
            df_sample_wide_null <- df_sample_wide
            for (i in 1:N_measures) {
                df_sample_wide_null[N_measures * 2 + i + 1] <- df_sample_wide_null[N_measures + i + 1] + es_vector_null[,i]
            }
            
            # Add to df_sample
            df_sample_null <- df_sample
            for (i in 1:N_measures) {
                df_sample_null[df_sample$measure == i,4] <- df_sample_wide_null[N_measures * 2 + i + 1]
            }
            names(df_sample) <- c("id", "measure", "y", "y_post")
            
            
            # Calculate average
            df_sample_average_null <- df_sample_average
            df_sample_average_null$y_post <- aggregate(df_sample_null$y_post, list(df_sample_null$id), mean)$x
            
            # Calculate scaled average
            tmp <- df_sample_null
            tmp$y <- tmp$y_post
            df_sample_post_wide_null <- long2wide(tmp, within = "measure")
            mean_cormat_sample_post_null <- sum(cor(df_sample_post_wide_null[,2:ncol(df_sample_post_wide_null)]))
            grand_mean_sample_post_null <- mean(df_sample_null$y_post)
            sample_post_sum_null <- aggregate(df_sample_null$y_post, list(df_sample_null$id), sum)
            
            df_sample_scaled_average_null <- df_sample_scaled_average
            df_sample_scaled_average_null$y_post <- (sample_post_sum_null$x - N_measures * grand_mean_sample_post_null) / sqrt(mean_cormat_sample_post_null) + grand_mean_sample_post_null
            
            # Calculate EFA
            df_sample_fa_coarse_null <- df_sample_fa_coarse
            df_sample_fa_coarse_null$y_post <- rowSums(mapply(`*`, df_sample_post_wide_null[,2:ncol(df_sample_post_wide_null)], fa_result$loadings))
            
            df_sample_fa_refined_null <- df_sample_fa_refined
            df_sample_fa_refined_null$y_post <- rowSums(mapply(`*`, df_sample_post_wide_null[,2:ncol(df_sample_post_wide_null)], fa_result$weights))
            
            
            # Combine dataframes
            df_sample_all_null <- rbind(df_sample_null, df_sample_average_null, df_sample_scaled_average_null, df_sample_fa_coarse_null, df_sample_fa_refined_null)
            df_sample_all_null$y_diff <- df_sample_all_null$y_post - df_sample_all_null$y
            
            
            # T-tests ----
            for (i_test in levels(df_sample_all$measure)) {
                i_test <- noquote(i_test)
                test <- t.test(df_sample_all[df_sample_all$measure == i_test,]$y_diff, mu = 0, alternative = "greater")
                test_null <- t.test(df_sample_all_null[df_sample_all_null$measure == i_test,]$y_diff, mu = 0, alternative = "greater")
                
                df_p_values[i_sim, i_test] <- test$p.value
                df_p_values_null[i_sim, i_test] <- test_null$p.value
            }
        }
        return(list(df_p_values = df_p_values, df_p_values_null = df_p_values_null, N_people = N_people, N_measures = N_measures, es = es))
    }
    
    v <- reactiveValues(data = NULL)

    observeEvent(input$test_specific_effect, {
        v$data <- 1
    })

    observeEvent(input$start_sim, {
        v$data <- 0
    })
    
    button_sim <- eventReactive(input$start_sim, {
        run_sim()
    })
    
    output$simPlot <- renderPlot({
        
        if (v$data == 1) return()
        df <- button_sim()
        
        prim_outcome <- 2
        
        pwr <- pwr.t.test(n = df$N_people, d = df$es, type = "one", alternative = "greater")
        pwr_null <- pwr.t.test(n = df$N_people, d = 0, type = "one", alternative = "greater")
        
        
        # Effect present ----
        prop_sig <- vector()
        prop_sig_bonf <- vector()
        for (i in 1:df$N_measures) {
            df$df_p_values[,df$N_measures + 3 + i] <- df$df_p_values[,i] * df$N_measures
            prop_sig[i] <- mean(df$df_p_values[,i] < .05)
            prop_sig_bonf[i] <- mean(df$df_p_values[,df$N_measures + 3 + i] < .05)
        }
        
        positives_rate_across <- mean(ifelse(rowSums(ifelse(df$df_p_values[,1:df$N_measures] < .05, 1, 0)) > 0, 1, 0))
        positives_rate_across_bonf <- mean(ifelse(rowSums(ifelse(df$df_p_values[,(df$N_measures + 4):ncol(df$df_p_values)] < .05, 1, 0)) > 0, 1, 0))
        
        positives_rate_average <- mean(df$df_p_values$average < .05)
        positives_rate_scaled_average <- mean(df$df_p_values$scaled_average < .05)
        positives_rate_fa_coarse <- mean(df$df_p_values$fa_coarse < .05)
        positives_rate_fa_refined <- mean(df$df_p_values$fa_refined < .05)
        
        
        positives_rates <- data.frame(measure = c("Across measures", "Across measures, corr", paste0("Single measure"), "Average", "Scaled average", "EFA coarse", "EFA refined"),
                                      rate = c(positives_rate_across, positives_rate_across_bonf, prop_sig[prim_outcome], positives_rate_average, positives_rate_scaled_average, positives_rate_fa_coarse, positives_rate_fa_refined))
        
        positives_rates$measure <- factor(positives_rates$measure, levels = c("Across measures", "Across measures, corr", paste0("Single measure"), "Average", "Scaled average", "EFA coarse", "EFA refined"))
        
        positives_rates$section <- positives_rates$measure
        levels(positives_rates$section) <- c("Multiple outcomes", "Multiple outcomes", "Primary outcome", "Composite", "Composite", "Composite", "Composite")
        
        
        # Effect not present ----
        prop_sig_null <- vector()
        prop_sig_bonf_null <- vector()
        for (i in 1:df$N_measures) {
            df$df_p_values_null[,df$N_measures + 3 + i] <- df$df_p_values_null[,i] * df$N_measures
            prop_sig_null[i] <- mean(df$df_p_values_null[,i] < .05)
            prop_sig_bonf_null[i] <- mean(df$df_p_values_null[,df$N_measures + 3 + i] < .05)
        }
        
        positives_rate_across_null <- mean(ifelse(rowSums(ifelse(df$df_p_values_null[,1:df$N_measures] < .05, 1, 0)) > 0, 1, 0))
        positives_rate_across_bonf_null <- mean(ifelse(rowSums(ifelse(df$df_p_values_null[,(df$N_measures + 4):ncol(df$df_p_values_null)] < .05, 1, 0)) > 0, 1, 0))
        
        positives_rate_average_null <- mean(df$df_p_values_null$average < .05)
        positives_rate_scaled_average_null <- mean(df$df_p_values_null$scaled_average < .05)
        positives_rate_fa_coarse_null <- mean(df$df_p_values_null$fa_coarse < .05)
        positives_rate_fa_refined_null <- mean(df$df_p_values_null$fa_refined < .05)
        
        
        positives_rates_null <- data.frame(measure = c("Across measures", "Across measures, corr", paste0("Single measure"), "Average", "Scaled average", "EFA coarse", "EFA refined"),
                                           rate = c(positives_rate_across_null, positives_rate_across_bonf_null, prop_sig_null[prim_outcome], positives_rate_average_null, positives_rate_scaled_average_null, positives_rate_fa_coarse_null, positives_rate_fa_refined_null))
        
        positives_rates_null$measure <- factor(positives_rates_null$measure, levels = c("Across measures", "Across measures, corr", paste0("Single measure"), "Average", "Scaled average", "EFA coarse", "EFA refined"))
        
        positives_rates_null$section <- positives_rates_null$measure
        levels(positives_rates_null$section) <- c("Multiple outcomes", "Multiple outcomes", "Primary outcome", "Composite", "Composite", "Composite", "Composite")
        
        
        plot_effect <- ggplot(positives_rates, aes(y = rate, x = measure, fill = section)) +
            geom_bar(stat = "identity") +
            facet_grid(.~section, scales = "free", space = "free") +
            #ylim(0,1) +
            geom_hline(yintercept = pwr$power, linetype = "dashed", colour = "grey", size = 1) +
            theme_fivethirtyeight() +
            scale_fill_brewer(palette = "Pastel2") +
            theme(panel.grid = element_blank(),
                  axis.line = element_line(colour = "grey"),
                  axis.ticks = element_line(colour = "grey"),
                  legend.position = "none") +
            labs(title    = "Statistical power",
                 subtitle = "This plot shows the probability of detecting the specified intervention effect for three different scenarios: \nAcross all outcome measures (green; without correction and with Bonferroni correction), for a single primary outcome measure (orange), \nand for a number of composite scores (blue; average, scaled average, coarse factor scores, refined (Thurstone) factor scores). \nThe gray dashed line indicates the power of a single measure, given the sample size and effect size, as estimated by a power analysis.")
        
        plot_no_effect <- ggplot(positives_rates_null, aes(y = rate, x = measure, fill = section)) +
            geom_bar(stat = "identity") +
            facet_grid(.~section, scales = "free", space = "free") +
            #ylim(0,1) +
            geom_hline(yintercept = pwr_null$power, linetype = "dashed", colour = "grey", size = 1) +
            theme_fivethirtyeight() +
            scale_fill_brewer(palette = "Pastel2") +
            theme(panel.grid = element_blank(),
                  axis.line = element_line(colour = "grey"),
                  axis.ticks = element_line(colour = "grey"),
                  legend.position = "none") +
            labs(title    = "False positives",
                 subtitle = "This plot assumes an effect size of 0 and shows the rate of false positives for three different scenarios: \nAcross all outcome measures (green; without correction and with Bonferroni correction), for a single primary outcome measure (orange), \nand for a number of composite scores (blue; average, scaled average, coarse factor scores, refined (Thurstone) factor scores). \nThe gray dashed line indicates an error rate of .05.")
        
        plot_task_effect <- ggplot(positives_rates, aes(y = rate, x = measure, fill = section)) +
            geom_bar(stat = "identity") +
            facet_grid(.~section, scales = "free", space = "free") +
            #ylim(0,1) +
            #geom_hline(yintercept = pwr$power, linetype = "dashed", colour = "grey", size = 1) +
            theme_fivethirtyeight() +
            scale_fill_brewer(palette = "Pastel2") + 
            theme(panel.grid = element_blank(),
                  axis.line = element_line(colour = "grey"),
                  axis.ticks = element_line(colour = "grey"),
                  legend.position = "none") +
            labs(title    = "False positives",
                 subtitle = "This plot shows the rate of false positives for three different scenarios: \nAcross all outcome measures (green; without correction and with Bonferroni correction), for a single primary outcome measure (orange), \nand for a number of composite scores (blue; average, scaled average, coarse factor scores, refined (Thurstone) factor scores). \nThe gray dashed line indicates an error rate of .05.")

        switch(input$test_specific_effect,
                        "a specific task" = plot_grid(plot_task_effect, NULL, ncol = 1),
                        "the construct" = plot_grid(plot_effect, plot_no_effect, ncol = 1))

    })
}

# Run the application 
shinyApp(ui = ui, server = server)
