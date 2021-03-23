library(tidyverse)
library(magrittr)

metadata <- read_csv(snakemake@input[[1]])


readStarLog <- function(log_file){

    out = list()
    lines = readLines(log_file)

    out$input_reads = (lines[6] %>% strsplit('\t') %>% unlist)[2] %>% as.integer

    out$uniq_mapped_reads = (lines[9] %>% strsplit('\t') %>% unlist)[2] %>% as.integer

    #out$avg_length = (lines[11] %>% strsplit('\t') %>% unlist)[2] %>% as.numeric
    
    tibble(observation=names(out), value=unlist(unname(out)))
}

read_metrics <- metadata %>%
    select(project_id, sample_id, puck_id, species, sequencing_date) %>%
    mutate(star_log = paste0('/data/rajewsky/projects/slide_seq/projects/', project_id, '/processed_data/', sample_id, '/illumina/complete_data/star_Log.final.out'),
           read_types =paste0('/data/rajewsky/projects/slide_seq/projects/', project_id, '/processed_data/', sample_id, '/illumina/complete_data/split_reads/read_type_num.txt')) %>%

    filter(file.exists(star_log), file.exists(read_types)) %>%
    mutate(star_log = map(star_log,
                          ~ readStarLog(.))) %>%
    unnest(star_log) %>%
    mutate(read_types = map(read_types,
                            ~ read_table2(., col_names=c('rt_obs', 'rt_value')))) %>%
    unnest(read_types) %>%
    mutate(rt_obs = tolower(rt_obs)) %>%
    spread(rt_obs, rt_value) %>%
    spread(observation, value)

read_metrics %>%
    write_delim(snakemake@output[[1]], '\t')
