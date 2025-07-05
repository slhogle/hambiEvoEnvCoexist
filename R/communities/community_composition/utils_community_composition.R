# Pair plotting functions -------------------------------------------------

# this function prepares masterplate t0 samples for plotting

format_master <- function(.data){
  .data %>% 
    # make a combined evolution and species identifier and extract the community ID
    dplyr::mutate(sp = paste(str_to_upper(evo_hist), str_extract(strainID, "\\d+"), sep = "_"),
                  community_id = str_extract(sample, "P\\d\\d")) %>%
    # this step is important to ensure that when dfs are pivoted wider the 
    # sp_1 and sp_2 stay consistent 
    dplyr::arrange(community_id, sp) %>%
    # this creates an index for each species present in each community, it is needed
    # for the pivot to be consistent between the master plate and the samples
    dplyr::group_by(community_id) %>% 
    dplyr::mutate(n = 1:n()) %>% 
    dplyr::ungroup() %>% 
    tidyr::pivot_wider(id_cols = community_id, values_from = c(sp, f), names_from = n) %>% 
    dplyr::mutate(transfer = 0) %>% 
    tidyr::expand_grid(strep_conc = c(0, 16, 64, 256)) %>% 
    arrange(sp_1, sp_2, strep_conc) %>% 
    arrange(community_id)
}  

# Binomial sampling and Wilcox test

wc_test <- function(df1, df2){
  # first join the T0 and T8 abundances
  left_join(df1, df2, by = join_by(community_id, strep_conc, sp_1, sp_2)) %>% 
    dplyr::mutate(delta_f_1 = f_1.x - f_1.y) %>% 
    dplyr::select(community_id, strep_conc, sp_1, sp_2, delta_f_1, f_1_8 = f_1.x, f_1_0 = f_1.y, 
                  f_2_8 = f_2.x, f_2_0 = f_2.y) %>% 
    tidyr::nest(data = c(-community_id, -strep_conc)) %>%
    # samples 1000 draws from binomial distribution using f_a median as the probability of success
    dplyr::mutate(f_1_0_rs = purrr::map(data, \(x) map(1:100, \(i) sum(rbinom(1000, 1, x$f_1_0))/1000)),
                  f_1_8_rs = purrr::map(data, \(x) map(1:100, \(i) sum(rbinom(1000, 1, x$f_1_8))/1000))) %>% 
    tidyr::unnest(cols = c(data, f_1_0_rs, f_1_8_rs)) %>% 
    # nest the samples
    tidyr::nest(bs = c(f_1_0_rs, f_1_8_rs)) %>%
    # perform the wilcox test
    dplyr::mutate(wc = purrr::map(bs, \(i) wilcox.test(x = as.numeric(i$f_1_0_rs), y = as.numeric(i$f_1_8_rs)))) %>% 
    # tidy-ify the test output
    dplyr::mutate(tidy_wc = purrr::map(wc, \(x) broom::tidy(x))) %>% 
    tidyr::unnest(cols = c(tidy_wc)) %>% 
    # p-value adjust using bonferroni correction
    dplyr::mutate(p_adjusted = p.adjust(p.value, method = "bonferroni", n = n())) %>% 
    dplyr::arrange(strep_conc, sp_1, sp_2) %>% 
    # define whether change is significantly positive or negative
    dplyr::mutate(change = dplyr::case_when(p.value > 1e-5 ~ 0,
                                            sign(delta_f_1) == -1 & p.value <= 1e-5 ~ -1, 
                                            sign(delta_f_1) == 1 & p.value <= 1e-5 ~ 1)) 
}

# This function plots species pairs in a grid format where with the relative abundance
# of species on the y axis

pair_plot <- function(df){
  outcome_pal <- c("exclusion" = "red", 
                   "towards exclusion" = "red4",
                   "stable coexistence" = "dodgerblue", 
                   "coexistence no mutual invasion" = "dodgerblue4", 
                   "inconclusive" = "green3")
  
  pj <- ggplot2::position_jitterdodge(jitter.width=0.0,
                                      jitter.height = 0.0,
                                      dodge.width = 0.5,
                                      seed=9)
  
  evostatus <- case_when(str_unique(df$evo_group) == "both_anc" ~ "Both Ancestral", 
                         str_unique(df$evo_group) == "both_evo" ~ "Both Evolved",
                         TRUE ~ "Mixed histories")
  
  ggplot2::ggplot(df, aes(x = transfer, y = f_1, group = group)) +
    ggplot2::geom_hline(yintercept = 0, lty = 2, color = "grey70") +
    ggplot2::geom_hline(yintercept = 0.5, lty = 3, color = "grey70") +
    ggplot2::geom_hline(yintercept = 1, lty = 2, color = "grey70") +
    ggplot2::geom_linerange(aes(ymin = fmin_1, ymax = fmax_1, color = outcome), position = pj) +
    ggh4x::geom_pointpath(aes(color = outcome), position = pj, mult = 0.2) +
    ggplot2::scale_color_manual(values = outcome_pal) +
    ggplot2::facet_grid(sp_1 ~ sp_2) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = percent) +
    ggplot2::scale_x_continuous(breaks = c(0, 8)) +
    ggplot2::labs(x = "", y = "", color = "", subtitle = paste0(evostatus, ", ", unique(df$strep_conc), " μg/ml Streptomycin")) +
    ggplot2::theme_bw() + 
    ggplot2::theme(panel.grid = element_blank(),
                   strip.background = element_blank(),
                   legend.position = "noe", 
                   panel.border = element_blank(),
                   axis.text = element_text(size = 8),
                   strip.text = element_text(size = 8))
}




# Network plotting functions ----------------------------------------------

make_pairs <- function(pairs_df, sp){
  pairs_df %>%
    dplyr::filter(transfer == 8) %>%
    dplyr::distinct(evo_group, strep_conc, sp_1, sp_2) %>%
    dplyr::group_by(evo_group, strep_conc) %>%
    dplyr::count( {{ sp }} ) %>% 
    dplyr::rename(sp = {{ sp }} ) %>% 
    dplyr::ungroup()
}

make_nodes <- function(pairs_df, ...){
  sp_1 <- make_pairs(pairs_df, sp_1)
  sp_2 <- make_pairs(pairs_df, sp_2)
  
  total_games <- bind_rows(sp_1, sp_2) %>%
    dplyr::summarize(games = sum(n),
                     .by = c(sp, ...))
  
  win_games <- pairs_df %>%
    dplyr::filter(transfer == 8) %>%
    dplyr::filter(outcome == "exclusion") %>%
    dplyr::group_by(sp_1, sp_2, ...) %>%
    dplyr::mutate(win = case_when(f_1 >= f_2 ~ sp_1, f_2 > f_1 ~ sp_2)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(sp_1, sp_2, win, ...) %>%
    dplyr::group_by(...) %>%
    dplyr::count(win, name = "wins") %>%
    dplyr::rename(sp = win) %>%
    dplyr::ungroup()
  
  lose_games <- pairs_df %>%
    dplyr::filter(transfer == 8) %>%
    dplyr::filter(outcome == "exclusion") %>%
    dplyr::group_by(sp_1, sp_2, ...) %>%
    dplyr::mutate(loss = case_when(f_1 >= f_2 ~ sp_2, f_2 > f_1 ~ sp_1)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(sp_1, sp_2, loss, ...) %>%
    dplyr::group_by(...) %>%
    dplyr::count(loss, name = "losses") %>%
    dplyr::rename(sp = loss) %>%
    dplyr::ungroup()
  
  left_join(total_games, win_games, by = join_by(sp, ...)) %>% 
    dplyr::left_join(lose_games, by = join_by(sp, ...)) %>% 
    dplyr::mutate(across(everything(), ~replace_na(.x, 0))) %>% 
    dplyr::mutate(score = (wins - losses)/games) %>% 
    dplyr::group_by(...) %>% 
    dplyr::arrange(..., desc(score)) %>% 
    dplyr::mutate(rank = dense_rank(desc(score)),
                  plotrank = row_number(desc(score)),
                  id = 1:n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(name = sp, ..., rank, plotrank)
}

make_edges <- function(pairs_df, ...){
  pairs_df %>% 
    dplyr::filter(transfer == 8) %>% 
    dplyr::mutate(from = if_else(f_1 > f_2, sp_1, sp_2),
                  to = if_else(f_1 < f_2, sp_1, sp_2)) %>% 
    dplyr::group_by(from, to, ...) %>% 
    #dplyr::mutate(type = if_else(sum(type == "stable") == 2, "stable", "incomplete")) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(from, to, ..., outcome) %>% 
    dplyr::distinct() %>% 
    dplyr::arrange(from, to, ...)
}

plot_network_hierarchy <- function(net, tune_angle = 1, n_rank = 10, n_break = 10) {
  # code for formatting the positions of the nodes with the ranks was taken from
  # here: https://github.com/Chang-Yu-Chang/emergent-coexistence/blob/v2.0.0/plotting_scripts/Fig3.R
  
  # outcome_pal <- c("exclusion" = "#800020", "coexistence" = "#0659bf", 
  #                  "inconclusive" = "#7c26c7", "bistable" = "#098f07")
  # type_pal <- c("stable" = "solid", "incomplete" = "dashed", 
  #               "potential" = "dotted")
  outcome_pal <- c("exclusion" = "red", 
                   "towards exclusion" = "red4",
                   "stable coexistence" = "dodgerblue", 
                   "coexistence no mutual invasion" = "dodgerblue4", 
                   "inconclusive" = "green3")
  
  node_size <- 3
  edge_width <- 0.8
  
  graph_ranked <- net %>%
    tidygraph::activate(nodes) %>%
    dplyr::select(name, rank, plotrank) %>%
    tidygraph::activate(edges) %>%
    dplyr::mutate(fromRank = .N()$plotrank[match(from, .N()$name)],
                  toRank = .N()$plotrank[match(to, .N()$name)])
  
  graph_ranked <- graph_ranked %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(y = -rank) %>%
    dplyr::group_by(rank) %>%
    dplyr::mutate(x = {seq(0, 1, length.out = n() + 2) %>% `[`(c(-1, -length(.)))}) %>%
    dplyr::ungroup() %>%
    tidygraph::activate(edges) %>%
    dplyr::filter(!str_detect(outcome, "inconclusive")) %>%
    dplyr::arrange(outcome)
  
  ggraph(graph_ranked, layout = "nicely") +
    geom_hline(yintercept = c(-n_rank:-1), color = "grey90") +
    geom_node_text(aes(label = name), repel = TRUE) +
    geom_node_point(size = node_size, shape = 21, fill = "grey", stroke = node_size/5, color = "black") +
    geom_edge_diagonal(aes(color = outcome),
                       arrow = arrow(length = unit(1, "mm"), type = "closed", angle = 30, ends = "last"),
                       start_cap = circle(node_size*.8, "mm"),
                       end_cap = circle(node_size*0.8, "mm")) +
    scale_edge_color_manual(values = outcome_pal) +
    #scale_edge_linetype_manual(values = type_pal) +
    scale_x_continuous(limits = c(0.1, 0.9), expand = c(0,0)) +
    scale_y_continuous(limits = c(-n_break-1, 0), breaks = -n_break:-1, labels = n_break:1) +
    theme_void() +
    theme(
      legend.title = element_blank(),
      axis.title = element_blank(),
      strip.text = element_blank(),
      plot.margin = unit(c(0,0,0,0),"mm")
    )
}

# Trio plotting functions -------------------------------------------------

trio_plot <- function(df) {
  spcols <- c(
    "HAMBI_0403" = "#bd7811",
    "HAMBI_1287" = "#476c9e",
    "HAMBI_1896" = "#31752a",
    "HAMBI_1977" = "#ffc755"
  )
  
  pj <- ggplot2::position_jitterdodge(
    jitter.width = 0.0,
    jitter.height = 0.0,
    dodge.width = 0.5,
    seed = 9
  )
  
  ggplot(df, aes(
    x = transfers,
    y = y,
    group = interaction(strainID, f0, strep_conc, hist)
  )) +
    geom_hline(yintercept = 0.01,
               color = "grey20",
               lty = 2) +
    geom_hline(yintercept = 0.99,
               color = "grey20",
               lty = 2) +
    ggplot2::geom_linerange(aes(
      ymin = ymin,
      ymax = ymax,
      color = strainID
    ), position = pj) +
    ggh4x::geom_pointpath(aes(color = strainID, shape = extinct), position = pj, mult = 0.2) +
    facet_grid(hist ~ strep_conc) +
    ggplot2::labs(x = "Growth cycle", y = "Species frequency", color = "Species") +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = percent
    ) +
    ggplot2::scale_x_continuous(limits = c(-1, 9), breaks = c(0, 8)) +
    scale_shape_manual(values = c(16, 4), guide = "none") +
    scale_color_manual(values = spcols) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      strip.background = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.grid = element_blank(),
      axis.text = element_text(size = 8),
      strip.text = element_text(size = 8)
    )
}

tern_plot_func <- function(df){
  ggtern(df, aes(f_1, f_2, f_3, fill = facetrow, shape = factor(f0))) +
    geom_point(alpha = 0.7, size = 3, position = position_jitter_tern(x = 0.05, y = 0.05, z = 0.05)) +
    ggplot2::scale_fill_manual(values = c("grey20", "#0C59FEFF", "#50D9C1", "#F0D446FF", "#FF8D00FF"), guide="none") +
    scale_shape_manual(values = c(21, 22, 24)) +
    labs(shape="",
         x = "Sp1",
         y = "Sp2",
         z = "Sp3") +
    facet_grid(facetrow~hist) +
    theme(tern.panel.mask.show = FALSE,
          strip.background = element_blank(),
          legend.position = "bottom",
          axis.text = element_text(size=6),
          axis.title = element_text(size = 8),
          tern.panel.grid.minor  = element_blank()) +
    theme_arrowlong()
}


# Assembly rule violation -------------------------------------------------

# Both these functions seem like a very clunky way of doing this... I can think
# of a way in python to do this way faster/simpler, but perhaps I am not
# thinking in the "tidyverse" way. Anyway...

matrixify <- function(df){
  # this function takes a "long" format tibble of species pairs and their
  # coexistence outcome and turns into a symmetric matrix with coexistence == 1
  # and exclusion == 0. It manually adds two species pairs on the matrix
  # diagonal that were missing in the original tibble thus constructing a
  # symmetrical matrix
  df %>% 
    add_row(sp_2 = "ANC_0403", sp_1 = "ANC_0403") %>% 
    add_row(sp_2 = "EVO_1977", sp_1 = "EVO_1977") %>% 
    complete(sp_1, sp_2, fill = list(coexists = NA)) %>% 
    pivot_wider(names_from = sp_1, values_from = coexists) %>% 
    data.frame() %>% 
    column_to_rownames(var = "sp_2") %>% 
    as.matrix()
}

coexistence_rule <- function(coexisting, competition_matrix){
  # This function takes as input a character vector of coexisting species from
  # an experimental trio, quartet or octet and then generates all possible
  # unique pairs of observed coexisting species. It then checks those species
  # pairs in the corresponding pairwise coexistence matrix and records whether
  # those species coexist (1) or not (0) in isolated pairs. 0 values indicate a
  # violation of the assembly rule
  
  # `coexisting` is a character vector of species found at >= 0.01 % <= 0.99 in 
  # a n>2 community 
  # `competition_matrix` is the symmetric binary matrix that records whether
  # species coexist (1) or not (0) in isolated pairs
  
  pairs <- asplit(combn(coexisting, 2), 2)
  tt <- list()
  for (i in 1:length(pairs)){
    sp_1 <- pairs[[i]][1]
    sp_2 <- pairs[[i]][2]
    tt[[i]] <- tibble(community_observed = paste(coexisting, collapse = ", "), 
                      pair_sp1 = sp_1, 
                      pair_sp2 = sp_2, 
                      upper = competition_matrix[sp_1, sp_2], 
                      lower = competition_matrix[sp_2, sp_1]) %>% 
      dplyr::mutate(pair_coexists_alone = dplyr::coalesce(upper, lower)) %>% 
      dplyr::select(-upper, -lower)
  }
  bind_rows(tt)
}

plot_heatmap <- function(df, label){
  ggplot(df) +
    geom_tile(aes(x = strep_conc, y = sp_set_id, fill = perc)) +
    geom_text(aes(x = strep_conc, y = sp_set_id, label = {{label}}, color = perc > 40)) +
    scale_fill_viridis_c(limits = c(0,100)) +
    scale_color_manual(values = c("white", "black"), guide = "none") +
    scale_x_discrete(guide = guide_axis(angle = 0)) +
    labs(x = "Streptomycin concentration (ug/ml)", y = "Species set",
         fill = "Percent\nsamples violating\nassembly rule") +
    coord_fixed() + 
    ggplot2::theme(panel.grid = element_blank(),
                   panel.background = element_blank(), 
                   strip.background = element_blank(),
                   panel.border = element_blank())
}
