tax_locus_copynum <- tibble::tribble(
  ~strainID, ~rRNA16S_cn, ~rRNA16S_locus,             ~genus,        ~species,
  "HAMBI_0006",          7L,  "H0006_04757",      "Pseudomonas",        "putida",
  "HAMBI_0097",          7L,  "H0097_00044",    "Acinetobacter",     "johnsonii",
  "HAMBI_0097",          7L,  "H0097_02759",    "Acinetobacter",     "johnsonii",
  "HAMBI_0097",          7L,  "H0097_01762",    "Acinetobacter",     "johnsonii",
  "HAMBI_0105",          4L,  "H0105_02306",    "Agrobacterium",   "tumefaciens",
  "HAMBI_0262",          3L,  "H0262_00030",    "Brevundimonas",       "bullata",
  "HAMBI_0403",          9L,  "H0403_00517",        "Comamonas",  "testosteroni",
  "HAMBI_0403",          9L,  "H0403_00522",        "Comamonas",  "testosteroni",
  "HAMBI_1279",          7L,  "H1279_03627",           "Hafnia",         "alvei",
  "HAMBI_1279",          7L,  "H1279_00125",           "Hafnia",         "alvei",
  "HAMBI_1279",          7L,  "H1279_03957",           "Hafnia",         "alvei",
  "HAMBI_1287",          7L,  "H1287_03997",      "Citrobacter",        "koseri",
  "HAMBI_1287",          7L,  "H1287_03402",      "Citrobacter",        "koseri",
  "HAMBI_1292",          7L,  "H1292_03239",       "Morganella",      "morganii",
  "HAMBI_1299",          8L,  "H1299_04293",         "Kluyvera",    "intermedia",
  "HAMBI_1299",          8L,  "H1299_01283",         "Kluyvera",    "intermedia",
  "HAMBI_1299",          8L,  "H1279_03957",         "Kluyvera",    "intermedia",
  "HAMBI_1842",          4L,  "H1842_01650",      "Sphingobium",    "yanoikuyae",
  "HAMBI_1896",          4L,  "H1896_00963", "Sphingobacterium",  "spiritivorum",
  "HAMBI_1972",         10L,  "H1972_00343",        "Aeromonas",        "caviae",
  "HAMBI_1972",         10L,  "H1972_03531",        "Aeromonas",        "caviae",
  "HAMBI_1977",          5L,  "H1977_00118",      "Pseudomonas",  "chlororaphis",
  "HAMBI_1988",          5L,  "H1988_05160",     "Chitinophaga",        "sancti",
  "HAMBI_1988",          5L,  "H1988_05152",     "Chitinophaga",        "sancti",
  "HAMBI_1988",          5L,  "H1988_05165",     "Chitinophaga",        "sancti",
  "HAMBI_2159",          4L,  "H2159_01406",        "Trinickia",   "caryophylli",
  "HAMBI_2159",          4L,  "H2159_05851",        "Trinickia",   "caryophylli",
  "HAMBI_2160",          3L,  "H2160_00530",       "Bordetella",         "avium",
  "HAMBI_2164",          5L,  "H2164_03337",      "Cupriavidus",    "oxalaticus",
  "HAMBI_2443",          3L,  "H2443_00128",       "Paracoccus", "denitrificans",
  "HAMBI_2494",          4L,  "H2494_03389", "Paraburkholderia",   "kururiensis",
  "HAMBI_2659",          4L,  "H2659_00367", "Stenotrophomonas",   "maltophilia",
  "HAMBI_2792",          4L,  "H2792_00549",        "Moraxella",         "canis",
  "HAMBI_3031",          2L,  "H3031_00830",         "Niabella",  "yanshanensis",
  "HAMBI_3237",          6L,  "H3237_00875",       "Microvirga",   "lotononidis",
  "HAMBI_1923",          6L,  "H1923_00876",   "Flavobacterium",      "odoratum"
)
            
# Amplicon normalization functions ----------------------------------------

# this function normalizes the abundance of reads for each species by the total number of
# rRNA genes in its genome
normalize_by_copy <- function(.data, tlc = tax_locus_copynum){
  .data %>% 
    # join with the copy number data frame. We join by the locus tag so this will add H1279_03957 to HAMBI_1299
    dplyr::left_join(tlc, by = join_by(rRNA16S_locus), relationship = "many-to-many") %>%
    # get total number of mapping reads per species. This aggregates all the difference ASVs per species
    dplyr::summarize(count = sum(count), .by = c(sample, strainID, rRNA16S_cn)) %>% 
    # group by sample
    dplyr::group_by(sample) %>% 
    # calculate a corrected count which is simply the count divided by copy num for each species
    # dividide by the sum of count divided by copy num for whole sample multiplied by the total
    # number of mapped reads per sample
    dplyr::mutate(count_correct = round(sum(count)*(count/rRNA16S_cn)/sum(count/rRNA16S_cn))) %>%  
    dplyr::ungroup() %>% 
    dplyr::select(sample, strainID, count, count_correct)
}

# this function replaces missing species counts with zero
complete_combos <- function(.data, tlc = tax_locus_copynum, countname = count, remove1923 = TRUE){
  
  # get unique strainIDs
  strainID <- unique(tlc$strainID)
  # table for assigning genus and species names. Doesn't matter if 1923 is there or not
  # because it is filter joined later
  tax <- dplyr::distinct(dplyr::select(tlc, strainID, genus, species))
  if (remove1923) {
    # get unique strainIDs but exclude 1923 if remove1923 is true
    strainID <- strainID[strainID != "HAMBI_1923"]
  }
  
  dplyr::bind_rows(tibble::tibble(strainID = strainID, sample = "dummy"), .data) %>% 
    dplyr::mutate( "{{ countname }}" := dplyr::if_else(sample == "dummy", 1, {{ countname }})) %>% 
    tidyr::complete(sample, strainID) %>% 
    dplyr::filter(sample != "dummy") %>% 
    dplyr::mutate( "{{ countname }}" := dplyr::if_else(is.na({{ countname }}), 0, {{ countname }})) %>% 
    tidyr::replace_na(list(count_correct = 0)) %>% 
    dplyr::left_join(dplyr::distinct(dplyr::select(tlc, strainID, genus, species)),
                     by = join_by(strainID)) %>% 
    dplyr::relocate(genus, species, .after = strainID)
}


# Amplicon plotting functions ---------------------------------------------

# This function plots a barplot with sample on the x-axis and species
# proportions on the y axis. If first takes a dataframe and then filters it to
# invlude only samples with species that shouldn't be in the sample and those
# species should be > 0%, Species fill is based on whether the species should be
# in the sample (orange, yellow, blue, green) or should not be there (grey
# through black).

contaminated_barplot <- function(df, threshold=0, quartet = FALSE, y){
  # colors for the different species
  spcols <- c("HAMBI_0403_anc" = "#faa019",
              "HAMBI_0403_evo" = "#bd7811",
              "HAMBI_1287_anc" = "#75afff",
              "HAMBI_1287_evo" = "#476c9e",
              "HAMBI_1896_anc" = "#59cc4e",
              "HAMBI_1896_evo" = "#31752a",
              "HAMBI_1977_anc" = "#ffd430",
              "HAMBI_1977_evo" = "#ab8e1f",
              "HAMBI_0403_NA"  = "#e6e5e3",
              "HAMBI_1287_NA"  = "#bdbcbb",
              "HAMBI_1896_NA"  = "#8c8c8b",
              "HAMBI_1977_NA"  = "#333333"
  )
  
  if (quartet)
    df1 <- df
  else
    df1 <- semi_join(df, df %>% filter(is.na(evo_hist) & {{ y }} >= threshold), 
            by = join_by(sample)) 
  df1 %>%
    mutate(sp = paste0(strainID, "_", evo_hist)) %>% 
    ggplot() +
    geom_col(aes(x = sample, y = {{ y }}, fill = sp)) +
    labs(y = "Percent abundance", x = "", fill = "") +
    scale_fill_manual(values = spcols) +
    scale_y_continuous(labels = percent) +
    scale_x_discrete(guide = guide_axis(angle=90)) +
    theme_bw() + 
    theme(panel.grid = element_blank(), legend.position = "right")
}

# this function plots a histogram for different x variables
contam_histogram <- function(df, trans, x){
  xtrans <- function(trans = TRUE) {
    if (trans) {
      list(
        scale_x_continuous(
          trans = "log10",
          labels = percent,
          guide = guide_axis(angle = 90)
        ),
        annotation_logticks(sides = "b", color = "grey30")
      )
    }
    else {
      list(scale_x_continuous(labels = percent, guide = guide_axis(angle = 90)))
    }
  }
  df %>% 
    ggplot(aes(x = {{ x }})) +
    geom_histogram(aes(fill = strainID), bins = 20) +
    geom_vline(xintercept = 0.01, linetype = "dashed") +
    scale_fill_manual(values = hambi_colors) +
    xtrans(trans) + 
    facet_grid(~strainID) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          legend.position = "none")
}

