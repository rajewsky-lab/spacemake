readStarLog <- function(log_file){

		out = list()
		lines = readLines(log_file)
	
		out$input_reads = (lines[6] %>% strsplit('\t') %>% unlist)[2] %>% as.integer

		out$uniq_mapped_reads = (lines[9] %>% strsplit('\t') %>% unlist)[2] %>% as.integer

		#out$avg_length = (lines[11] %>% strsplit('\t') %>% unlist)[2] %>% as.numeric
		
        tibble(observation=names(out), value=unlist(unname(out)))
	}
