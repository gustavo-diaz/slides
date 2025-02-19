library(tidyverse)
library(ggplot2)
library(purrr)
library(estimatr)
library(emmeans)
library(here)

# Load data ----

df_dh <- readRDS("precision_retention/1_df_dh_cleaned.rds")
df_bg <- readRDS("precision_retention/1_df_bg_cleaned.rds")
df_th <- readRDS("precision_retention/1_df_th_cleaned.rds")

# Pre-treatment and post-treatment correlation -----
dh_cor <- lm_robust(speech_approve_scaled ~ dv_quasipre_scaled + Z, data = df_dh, subset = design %in% c(2,3))
bg_cor <- lm_robust(delpref_post ~ delpref_pre + Z, data = df_bg, subset = design %in% c(2, 3))
th_cor <- lm_robust(pension_followup_scaled ~ pension_pre_scaled + Z, data = df_th, subset = design %in% c(2, 3))

cor.test(df_dh$speech_approve_scaled[df_dh$design %in% c(2,3)],
         df_dh$dv_quasipre_scaled[df_dh$design %in% c(2,3)],
         use = "pairwise.complete.obs")
cor.test(df_bg$delpref_post[df_dh$design %in% c(2,3)],
         df_bg$delpref_pre[df_dh$design %in% c(2,3)],
         use = "pairwise.complete.obs")
cor.test(df_th$pension_followup_scaled[df_th$design %in% c(2,3)],
         df_th$pension_pre_scaled[df_th$design %in% c(2,3)],
         use = "pairwise.complete.obs")

# Estimates with full samples (may differ in size due) -----

## Design 1
dh_design1 <- lm_robust(speech_approve_scaled ~ Z, data = df_dh, subset = design == 1)

bg_design1 <- lm_robust(delpref_post ~ Z, data = df_bg, subset = design == 1)

th_design1 <- lm_robust(pension_followup_scaled ~ Z, data = df_th, subset = design == 1)


## Design 2
dh_design2 <- lm_robust(speech_approve_scaled ~ Z + dv_quasipre_scaled, data = df_dh, subset = design == 2)

bg_design2 <- lm_robust(delpref_post ~ Z + delpref_pre, data = df_bg, subset = design == 2)

th_design2 <- lm_robust(pension_followup_scaled ~ Z + pension_pre_scaled, data = df_th, subset = design == 2)


## Design 3
dh_design3 <- lm_robust(speech_approve_scaled ~ Z, fixed_effects = ~ block, data = df_dh, subset = design == 3)

bg_design3 <- lm_robust(delpref_post ~ Z, fixed_effects = ~ block, data = df_bg, subset = design == 3)

th_design3 <- lm_robust(pension_followup_scaled ~ Z + pension_pre_scaled, fixed_effects = ~ block_id, data = df_th, subset = design == 3)

## Extra estimation -- covariate adjustment
# dh_design4 <- lm_robust(speech_approve_scaled ~ Z + age, fixed_effects = ~ block, data = df_dh, subset = design == 3)
# 
# bg_design4 <- lm_robust(delpref_post ~ Z + age, fixed_effects = ~ block, data = df_bg, subset = design == 3)
# 
# th_design4 <- lm_robust(pension_followup_scaled ~ Z + pension_pre_scaled  + age, fixed_effects = ~ block_id, data = df_th, subset = design == 3)


## Pooled treatment effects
dh_pooled <- lm_robust(speech_approve_scaled ~ Z, data = df_dh)

bg_pooled <- lm_robust(delpref_post ~ Z, data = df_bg)

th_pooled <- lm_robust(pension_followup_scaled ~ Z, data = df_th)


## Is treatment effect distinct across designs? (p-values mentioned in SI text)
dh_design_int <- lm_robust(speech_approve_scaled ~ Z*as.factor(design), data = df_dh)
dh_design_int2 <- lm_robust(speech_approve_scaled ~ Z*factor(design, levels = c(2,3,1)), data = df_dh)

bg_design_int <- lm_robust(delpref_post ~ Z*as.factor(design), data = df_bg)
bg_design_int2 <- lm_robust(delpref_post ~ Z*factor(design, levels = c(2,3,1)), data = df_bg)

th_design_int <- lm_robust(pension_followup_scaled ~ Z*as.factor(design), data = df_th)
th_design_int2 <- lm_robust(pension_followup_scaled ~ Z*factor(design, levels = c(2,3,1)), data = df_th)


# Table of treatment effects across designs (SI Table X) -----

dh_mod_list <- list("Pooled" = dh_pooled, "Standard design" = dh_design1, "Pre-post" = dh_design2, "Block randomized" = dh_design3)

bg_mod_list <- list("Pooled" = bg_pooled, "Standard design" = bg_design1, "Pre-post" = bg_design2, "Block randomized" = bg_design3)

th_mod_list <- list("Pooled" = th_pooled, "Standard design" = th_design1, "Pre-post" = th_design2, "Block randomized" = th_design3)

modelsummary::modelsummary(dh_mod_list,
                           coef_map = c("Z" = "Treatment Effect",
                                        "dv_quasipre_scaled" = "Pre-treatment outcome measure",
                                        "(Intercept)" = "Intercept"),
                           output = "latex",
                           stars = c('*' = .05),
                           gof_map = list(list("raw" = "nobs", "clean" = "Observations", "fmt" = 0)),
                           add_rows = data.frame("Block fixed effects", "", "", "", "Yes"),
                           title = "Treatment Effect Estimates By Design, Dietrich and Hayes Replication \\label{tab:dhests}",
                           escape = F)

modelsummary::modelsummary(bg_mod_list,
                           coef_map = c("Z" = "Treatment Effect",
                                        "delpref_pre" = "Pre-treatment outcome measure",
                                        "(Intercept)" = "Intercept"),
                           output = "latex",
                           stars = c('*' = .05),
                           gof_map = list(list("raw" = "nobs", "clean" = "Observations", "fmt" = 0)),
                           add_rows = data.frame("Block fixed effects", "", "", "", "Yes"),
                           title = "Treatment Effect Estimates By Design, Bayram and Graham Replication \\label{tab:bgests}",
                           escape = F)

modelsummary::modelsummary(th_mod_list,
                           coef_map = c("Z" = "Treatment Effect",
                                        "pension_pre_scaled" = "Pre-treatment outcome measure",
                                        "(Intercept)" = "Intercept"),
                           output = "latex",
                           stars = c('*' = .05),
                           gof_map = list(list("raw" = "nobs", "clean" = "Observations", "fmt" = 0)),
                           add_rows = data.frame("Block fixed effects", "", "", "", "Yes"),
                           title = "Treatment Effect Estimates By Design, Tappin and Hewitt Replication \\label{tab:thests}",
                           escape = F)


# Simulation -----

# Paper org:

# To understand the effects of alternative designs on X, Y, Z, A, B, C we did the following. [Describe]

# First, we discuss explicit sample loss. Trivial in single-wave. Substantial in multi-wave
# and more for blocking. Little differential sample inclusion. NO differential post-treatment attrition.

# Second, we address implicit sample loss via simulation. Given the result on
# explicit sample loss pretty much suggest it is not an issue. This simulation is 
# arguably more important for the practical implications of these designs.
# Researchers are absolutely making decisions about sample size that are 
# affected by survey length.... However, in the set up of this exercise,
# we preregistered the decision to not "hard-code" implicit loss in by
# collecting LESS data for the longer surveys because they cost more, etc.
# Instead, we paid the same, aimed to collect the same N, and hcose to SIMULATE
# implicit loss.

# We start by wiping the slate clean of explicit loss.

# Create an indicator for sample loss (from script
# 2_descriptives-explicitsampleloss.R)
df_dh$any_loss <- (df_dh$Finished == 0)
df_bg$any_loss <- (df_bg$Finished == 0)
df_th$any_loss <- ((df_th$finished_w1 != 1) | 
                     (df_th$finished_w2 != 1) | 
                     (df_th$finished_w3 != 1) | 
                     is.na(df_th$start_w2) | 
                     is.na(df_th$start_w3))

# Need to make variable for block id consistent across studies
df_th <- df_th %>%
  mutate(block = block_id)


estimation_fun <- function(pct_loss, des, exp_data, dv, dv_pre, starting_n, exp){
  
  # estimate model under alternative designs
  if(des == 2){
    df_subset <- exp_data %>%
      # starting subset used to estimate treatment effects
      filter(design == 2 & !any_loss & !is.na(Z)) %>%
      # choose random portion to delete to get to starting n
      sample_n(starting_n, replace = F) %>%
      # now chose a random portion according to sim. pct_loss
      sample_n(ceiling(starting_n - starting_n*pct_loss), replace = F)
    
    formula_str <- paste(dv, "~ Z + ", dv_pre)    
    m <- estimatr::lm_robust(as.formula(formula_str),
                              data = df_subset,
                              subset = design == 2)
  }
  if(des == 3){
    df_subset <- exp_data %>%
      # starting subset used to estimate treatment effects
      filter(design == 3 & !any_loss & !is.na(Z)) %>%
      # choose random portion to delete to get to starting n
      sample_n(starting_n, replace = F) %>%
      # now chose a random portion according to sim. pct_loss
      sample_n(ceiling(starting_n - starting_n*pct_loss), replace = F)
    
    if(exp == "th"){
      formula_str <- paste(dv, "~ Z + ", dv_pre)  
    }else{
      formula_str <- paste(dv, "~ Z")  
    }

    m <- estimatr::lm_robust(as.formula(formula_str),
                              fixed_effects = ~ block,
                              data = df_subset, 
                              subset = design == 3)
  }

  # store outputs
  data.frame(exp = exp,
             design = des,
             n = m$n,
             coef = m$coef["Z"],
             se = m$std.error["Z"],
             lower_ci = m$conf.low["Z"],
             upper_ci = m$conf.high["Z"],
             pct_loss = pct_loss)
}

# run simulation
inputs <- expand.grid(pct_loss = seq(0, .5, .05),
                      des = c(2, 3))
n_iters <- 1000
set.seed(123)
sim_iters_dh <- map_dfr(1:nrow(inputs), function(i) {
  # Repeat estimation_fun n_iters times for each combination
  map_dfr(1:n_iters, ~ estimation_fun(pct_loss = inputs$pct_loss[i],
                                   des = inputs$des[i],
                                   exp_data = df_dh,
                                   dv = "speech_approve_scaled",
                                   dv_pre = "dv_quasipre_scaled",
                                   starting_n = 208,
                                   exp = "dh"))
})

sim_iters_bg <- map_dfr(1:nrow(inputs), function(i) {
  # Repeat estimation_fun n_iters times for each combination
  map_dfr(1:n_iters, ~ estimation_fun(pct_loss = inputs$pct_loss[i],
                                  des = inputs$des[i],
                                  exp_data = df_bg,
                                  dv = "delpref_post",
                                  dv_pre = "delpref_pre",
                                  starting_n = 299,
                                  exp = "bg"))
})
sim_iters_th <- map_dfr(1:nrow(inputs), function(i) {
  # Repeat estimation_fun n_iters times for each combination
  map_dfr(1:n_iters, ~ estimation_fun(pct_loss = inputs$pct_loss[i],
                                      des = inputs$des[i],
                                      exp_data = df_th,
                                      dv = "pension_followup_scaled",
                                      dv_pre = "pension_pre_scaled",
                                      starting_n = 608,
                                      exp = "th"))
})


# # check sample size created from 'pct_loss' is comparable
# sim_iters_all %>%
#   group_by(exp, pct_loss, design) %>%
#   summarise(n_mean = mean(n))

## Standard design, simulate many iters at starting sample size
## to estimate average SE
standard_sim_dh <- plyr::adply(1:n_iters, 1, function(i) {
  df_subset <- df_dh %>%
    # starting subset used to estimate treatment effects
    filter(design == 1 & !any_loss & !is.na(Z)) %>%
    # choose random portion to delete to get to starting n
    sample_n(208, replace = F)
  
  m <- lm_robust(speech_approve_scaled ~ Z,
                 data = df_subset,
                 subset = design == 1)
  
  # store outputs
  data.frame(n = m$n,
             coef = m$coef["Z"],
             se = m$std.error["Z"],
             lower_ci = m$conf.low["Z"],
             upper_ci = m$conf.high["Z"],
             exp = "dh")
})



standard_sim_bg <- plyr::adply(1:n_iters, 1, function(i) {
  df_subset <- df_bg %>%
    # starting subset used to estimate treatment effects
    filter(design == 1 & !any_loss & !is.na(Z)) %>%
    # choose random portion to delete to get to starting n
    sample_n(299, replace = F)
  
  m <- lm_robust(delpref_post ~ Z,
                 data = df_subset,
                 subset = design == 1)
  
  # store outputs
  data.frame(n = m$n,
             coef = m$coef["Z"],
             se = m$std.error["Z"],
             lower_ci = m$conf.low["Z"],
             upper_ci = m$conf.high["Z"],
             exp = "bg")
})


standard_sim_th <- plyr::adply(1:n_iters, 1, function(i) {
  df_subset <- df_th %>%
    # starting subset used to estimate treatment effects
    filter(design == 1 & !any_loss & !is.na(Z)) %>%
    # choose random portion to delete to get to starting n
    sample_n(608, replace = F)
  
  m <- lm_robust(pension_followup_scaled ~ Z,
                              data = df_subset,
                              subset = design == 1)
  
  # store outputs
  data.frame(n = m$n,
             coef = m$coef["Z"],
             se = m$std.error["Z"],
             lower_ci = m$conf.low["Z"],
             upper_ci = m$conf.high["Z"],
             exp = "th")
})



# Summaries

sim_iters_all <- bind_rows(sim_iters_dh, sim_iters_bg, sim_iters_th)

standard_sim_all <- bind_rows(standard_sim_dh,  standard_sim_bg, standard_sim_th)


# summary for standard design
standard_summary <- standard_sim_all %>%
  group_by(exp) %>%
  summarise(se_mean = mean(se))


# summary for designs 2 and 3 to plot
# clean up values for plot
sim_summary <- sim_iters_all %>%
  group_by(exp, pct_loss, design) %>%
  summarise(se_mean = mean(se),
            se_sd = sd(se),
            se_max = max(se),
            se_min = min(se),
            se_975 = quantile(se, .975),
            se_025 = quantile(se, .025)) %>%
  # percent change relative to standard design = 
  # [(alt design se - standard se) / standard se] × 100
  mutate(pct_change = case_when(
    exp == "dh" ~ ((se_mean - standard_summary$se_mean[standard_summary$exp == "dh"])/standard_summary$se_mean[standard_summary$exp == "dh"])*100,
    exp == "bg" ~ ((se_mean - standard_summary$se_mean[standard_summary$exp == "bg"])/standard_summary$se_mean[standard_summary$exp == "bg"])*100,
    exp == "th" ~ ((se_mean - standard_summary$se_mean[standard_summary$exp == "th"])/standard_summary$se_mean[standard_summary$exp == "th"])*100,
  )) %>%
  mutate(pct_change_025 = case_when(
    exp == "dh" ~ ((se_025 - standard_summary$se_mean[standard_summary$exp == "dh"])/standard_summary$se_mean[standard_summary$exp == "dh"])*100,
    exp == "bg" ~ ((se_025 - standard_summary$se_mean[standard_summary$exp == "bg"])/standard_summary$se_mean[standard_summary$exp == "bg"])*100,
    exp == "th" ~ ((se_025 - standard_summary$se_mean[standard_summary$exp == "th"])/standard_summary$se_mean[standard_summary$exp == "th"])*100,
  )) %>%
  mutate(pct_change_975 = case_when(
    exp == "dh" ~ ((se_975 - standard_summary$se_mean[standard_summary$exp == "dh"])/standard_summary$se_mean[standard_summary$exp == "dh"])*100,
    exp == "bg" ~ ((se_975 - standard_summary$se_mean[standard_summary$exp == "bg"])/standard_summary$se_mean[standard_summary$exp == "bg"])*100,
    exp == "th" ~ ((se_975 - standard_summary$se_mean[standard_summary$exp == "th"])/standard_summary$se_mean[standard_summary$exp == "th"])*100,
  )) %>%
  mutate(design = factor(design, levels = c(2,3),
                         labels = c("Pre-post", "Block Randomized (and Pre-Post for Tappin & Hewitt)"))) %>%
  mutate(exp = factor(exp, levels = c("dh", "bg", "th"),
                      labels = c("Dietrich & Hayes",
                                 "Bayram & Graham",
                                 "Tappin & Hewitt")))





implicit_loss_pctchange <- ggplot(sim_summary, aes(x = pct_loss,
                                               y = pct_change,
                                               color = factor(design),
                                               shape = factor(design))) +
  facet_wrap(~exp, scales = "fixed", nrow = 1) +
  geom_pointrange(aes(ymin = pct_change_025, ymax = pct_change_975), size = 0.3,
                  position = position_dodge2(width = 0.03)) +
  # Horizontal line at 0
  geom_hline(aes(yintercept = 0),
             linetype = "dashed",
             color = "black") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.minor.x = element_blank(),
    text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  labs(x = "Implicit Sample Loss",
       y = "Percentage Change in Standard Error\nRelative to Standard Design",
       color = "Design") + 
  scale_color_manual(values = c("grey50", "grey20")) +
  scale_shape_manual(values = c(16,17)) + 
  guides(color = guide_legend(title = "Design", override.aes = list(shape = c(16,17)),
                              ncol = 2),
         shape = FALSE) +
  scale_y_continuous(limits = c(-50, 90),
                         breaks = seq(-50, 90, by = 25),
                         labels = paste0(seq(-50, 90, by = 25), "%")) +
  scale_x_continuous(breaks = seq(0, .50, by = .10),
                     labels = paste0(seq(0, 50, by = 10), "%"))


implicit_loss_pctchange
ggsave("3_plot_implicit_loss_pctchange.pdf", implicit_loss_pctchange, width = 10, height = 4)




implicit_loss_rawSE <- ggplot(sim_summary, aes(x = pct_loss,
                                         y = se_mean,
                                         color = factor(design),
                                         shape = factor(design))) +
  facet_wrap(~exp, scales = "fixed", nrow = 1) +
  geom_pointrange(aes(ymin = se_025, ymax = se_975), size = 0.3,
                  position = position_dodge2(width = 0.03)) +
  geom_hline(data = subset(sim_summary, exp == "Dietrich & Hayes"),
             aes(yintercept = standard_summary$se_mean[standard_summary$exp == "dh"]),
             linetype = "dashed",
             color = "black") +
  geom_hline(data = subset(sim_summary, exp == "Bayram & Graham"),
             aes(yintercept = standard_summary$se_mean[standard_summary$exp == "bg"]),
             linetype = "dashed",
             color = "black") +
  geom_hline(data = subset(sim_summary, exp == "Tappin & Hewitt"),
             aes(yintercept = standard_summary$se_mean[standard_summary$exp == "th"]),
             linetype = "dashed",
             color = "black") +
  geom_text(data = sim_summary %>% filter(exp == "Dietrich & Hayes") %>% head(1),
            aes(x = 0, y = standard_summary$se_mean[standard_summary$exp == "dh"]),
            vjust = 1,
            hjust = 0,
            nudge_y = -.002,
            nudge_x = -.02,
            color = "black",
            label = "Standard Design SE") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    panel.grid.minor.x = element_blank(),
    text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 12),
    legend.position = "bottom"
  ) +
  labs(x = "Implicit Sample Loss",
       y = "Standard Error",
       color = "Design") + 
  scale_color_manual(values = c("grey50", "grey20")) +
  scale_shape_manual(values = c(16,17)) + 
  guides(color = guide_legend(title = "Design", override.aes = list(shape = c(16,17)),
                              ncol = 2),
         shape = FALSE)  +
  scale_x_continuous(breaks = seq(0, .50, by = .10),
                     labels = paste0(seq(0, 50, by = 10), "%"))

ggsave("3_plot_implicit_loss_rawSE.pdf", implicit_loss_rawSE, width = 10, height = 4)









