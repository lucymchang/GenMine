GenMine <- function(taxa, gene, batch = 10, retmax = 500, rettotal = "all", dir = "", suffix = "genOut.fna", api_key, slow = 0.5){
  # batch = how many taxa to search by at a time
  # retmax = number of records to pull for each taxon at a time
  # rettotal = maximum number of records for each taxon to keep
  # api_key (optional) = API key associated with NCBI account (allows for 10 requests/sec rather than default 3 requests/sec)
  # slow = number of second to wait between requests, so as not to receive 429 (too many requests) error
  
  require(xml2)
  require(ape)
  require(plyr)
  
  # print warnings as they happen
  options(warn = 1)
  
  if(missing(taxa)){
    stop("Needs vector of taxonomic names to search for.")
  }
  if(missing(gene)){
    stop('Needs vector of gene names or other terms in Entrez query format to search for (e.g. "h3[gene]", "28s+rrna[gene]").')
  }
  if(missing(api_key)){
    warning("No api_key. Request frequency is limited to 3 per second (as opposed to 10 per second). API key can be obtained from Settings page of NCBI account (http://www.ncbi.nlm.nih.gov/account/).")
    api_key <- ""
  }else{
    api_key <- paste0("&api_key=", api_key)
  }
  
  dir <- gsub("/$", "", dir)
  if(!dir.exists(dir)){
    dir.create(dir)
  }
  
  sapply(gene, function(gen){
    file <- paste0(dir, "/", paste(strsplit(gen, "[[:punct:]]")[[1]], collapse = "_"), suffix)
    file.create(file)
    taxa <- gsub(" ","+",taxa)
    found <- c()
    base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    
    taxa.batch <- split(taxa, ceiling(seq_along(taxa) / batch))
    
    webEnv <- ""
    key <- c()
    count <- c()
    
    found <- llply(taxa.batch, function(taxa.sub){
      
      cat(paste(gen, which(taxa == taxa.sub[1]), "to", which(taxa == taxa.sub[length(taxa.sub)]), "of", length(taxa), "taxa;", Sys.time(), "\n"))
      
      for(taxon in taxa.sub){
        
        # add a small break to minimize chances that NCBI will return 429 (too many requests) error
        Sys.sleep(slow)
        
        orgn <- paste0(taxon, "[orgn]")
        query <- paste(orgn, gen, sep = "+AND+")
        results <- read_xml(paste0(base, "esearch.fcgi?db=nucleotide&term=", query, "&usehistory=y&webEnv=", webEnv, api_key))
        results <- as_list(results)
        webEnv <- results[["WebEnv"]][[1]]
        key <- c(key, results[["QueryKey"]][[1]])
        count <- c(count, as.numeric(results[["Count"]][[1]]))
      }
      
      hits <- which(count != 0)
      
      if(length(hits) > 0){
        
        l_ply(hits, function(hit){
          
          # start from the beginning
          retstart = 0
          # collect up to either user defined number or all records
          rettotal.fin <- if(rettotal=="all") count[hit] else min(count[hit], rettotal)
          # if total number of records is less than the batch, go up to it instead
          retmax.fin <- min(rettotal.fin, retmax)
          
          # while there are still records to go
          while(retstart < rettotal.fin){
            
            closeAllConnections()
            
            # sequences is either a successful read or NA
            # NA means that batch failed to download
            # but it will move on to the next batch
            sequences <- tryCatch(
              {
                efetch <- url(paste0(base, "/efetch.fcgi?db=nucleotide&WebEnv=", webEnv, "&query_key=", key[hit], "&retstart=", retstart, "&retmax=", min(retmax.fin, rettotal.fin - retstart), "&rettype=fasta&retmode=text", api_key))
                readLines(efetch)
              },
              error = function(cond){
                message(cond)
                message(paste("Records", retstart, "to", retstart + retmax.fin, "for", gen, taxa.sub[hit]))
                return(NA)
              }
            )
            
            close(efetch)
            
            if(!anyNA(sequences)){
              cat(sequences, file = file, sep = "\n", append = T)
            }
            retstart <- retstart + retmax.fin
            
          }
          
        })
        
        return(taxa.sub[hits])
        
      }
      
    })
    
    unlist(found, use.names = F)
    
  })
}





