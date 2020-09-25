# rm(list=ls())
options(stringsAsFactors = F)
# setwd('//iron/prot_proteomics/LabMembers/Karsten/Projects/CPTAC3/PTRC/TNBC/20200327_ptm-sea/data/')

library(pacman)
p_load(cmapR)
p_load(seqinr)
p_load(glue)
p_load(magrittr)

tmp.dir <- tempdir()



# ## input gct
# gct.all.str <- c('TNBC_Phospho_Respective_CR_MedianMADnorm.gct')
# 
# 
# ## database
# ## db.str <- '//iron/proteomics_storage_slow/CPTAC3/RefSeq_20180629/RefSeq.20180629_Human_ucsc_hg38_cpdbnr_mito_264contams.fasta'
# #db.str <- '//cibola/seqdb/RefSeq.20180629_Human_ucsc_hg38_cpdbnr_mito_264contams_553smORFs.fasta'
# #db.str <- '//cibola/seqdb/RefSeq.20160914_Human_ucsc_hg19_customProDBnr_mito_150contams_553smORFS.fasta'
# db.str <- '../db/RefSeq.20171003_Human_ucsc_hg38_cpdb_mito_259contamsnr_553smORFS.fasta'
# 
# flank <- 7
# mod <- '-p'
# 
# ## sm - Spectrum Mill
# ## ph - Philosopher
# pipeline <- c('sm', 'ph', 'pnnl')[3]

#output file name
#flank.out <- "./VMSiteFlank_output.gct"

#Should I write results?
#write.output <- F

## Spectrum Mill
if(pipeline == 'sm'){
  acc.exprs <- '^((..|smORF)_.*\\.{0,1}.*?)_.*' ## extract protein accessions from site ids
  res.expr <- '^.*?_.*?_(.*?)_.*'                   ## extract residue numbers from site ids    
  db.exprs <- '^(.*)(\\|| ).*'                  ## parse protein ids from db names                
}
## Philosopher
if(pipeline == 'ph'){
  acc.exprs <- '^(.._.*\\..*?)_.*'              ## extract protein accessions from site ids
  res.expr <- '^.*_(.*)'                        ## extract residue numbers from site ids
  db.exprs <- '^(.*)(\\|| ).*'                  ## parse protein ids from db names                
}
if(pipeline == 'pnnl'){
  acc.exprs <- '^(.._.*\\..*?)_.*'              ## extract protein accessions from site ids
  res.expr <- '^.*_(.*)@.*'                      ## extract residue numbers from site ids
  db.exprs <- '^(.*)(\\|| ).*'                  ## parse protein ids from db names      
}



## ##########################
##         start
## ##########################

## database
db <- read.fasta(db.str, seqtype = 'AA')

## fix names
names(db) <- names(db) %>% sub(db.exprs, '\\1', .)
db.full <- db

#########################################
## loop over GCTs
for(gct.in.str in gct.all.str){
  
          db <- db.full

          ## data file
          gct.in <- parse.gctx(gct.in.str)
          rdesc <- gct.in@rdesc
          rid <- gct.in@rid %>% gsub(' ', '', .)
          mat <- gct.in@mat
          
          n <- ncol(mat)
          cn <- colnames(mat)
          
          ## clean up
          keep.idx <- apply(mat, 1, function(x) sum(!is.na(x)) )
          keep.idx <- which(keep.idx > 0)
          
          mat <- data.frame(mat[keep.idx, ])
          colnames(mat) <- cn
          
          rdesc <- rdesc[keep.idx, ]
          rid <- rid[keep.idx]  
          
          
          ## database accession numbers
          prot.id <- rid %>% sub(acc.exprs, '\\1', .) %>% unique
          
          if(sum(prot.id %in% names(db)) != length(prot.id)){ 
            rm.idx  <- which(!prot.id %in% names(db))
            warning('Could not find all protein accessions from datafile in database... First 10 accession numbers:\n\n')
            cat(paste(prot.id[rm.idx[1:min(length(rm.idx), 10)]], collapse='\n'), '\n')
            cat('Remove and continue?')
            tmp <- ''
            while(!tmp %in% c('y','n') )
             tmp <- readline(prompt='y/n: ')
            if(tmp =='n' )
               stop('Aborting\n')
            if(tmp == 'y'){
              mat <- mat[-rm.idx, ]
              rdesc <- rdesc[-rm.idx, ]
              rid <- rid[-rm.idx]
            }
            
          }
          
          db <- db[prot.id]
          
          rid.seqwin <- vector('character', length(rid))
          names(rid.seqwin) <- rid
          
          ## loop over sites
          for(acc in rid){
            
              acc.prot <- acc %>% sub(acc.exprs, '\\1', .)
              protseq <- db[[acc.prot]]
              l <- length(protseq)
             
              if(pipeline == 'ph'){
                acc.res.all <- acc %>% sub(res.expr,'\\1', .) %>% gsub('(S|T|Y)', ' ', .) %>% sub('^ ', '', .) %>% strsplit(., ' ') %>% unlist %>% as.numeric
                acc.aa.all <- acc %>% sub(res.expr,'\\1', .) %>%  gsub('[0-9]',' ', .) %>% strsplit(. , ' ') %>% unlist
                acc.aa.all <- acc.aa.all[nchar(acc.aa.all) > 0]
              }
              if(pipeline == 'sm' | pipeline == 'pnnl'){
                acc.res.all <- acc %>% sub(res.expr,'\\1', .) %>% gsub('(S|T|Y)', '', .) %>% strsplit(., 's|t|y') %>% unlist %>% as.numeric
                acc.aa.all <- acc %>% sub(res.expr,'\\1', .) %>% strsplit(., 's|t|y') %>% unlist %>% gsub('[0-9]', '', .)
              }
              
              if(length(acc.aa.all) > 0){
                  
                  seqwin <- c()
                  for(i in 1:length(acc.res.all)){
                    acc.res <- acc.res.all[i]
                    acc.aa <- acc.aa.all[i]
                      
                      ## N-term
                      prefix=''
                      if( (acc.res-flank) < 1 )
                        prefix <- paste(rep('_',abs(acc.res - flank)+1), collapse='')
                      
                      ## C-term
                      sufix <-''
                      if( (acc.res + flank) > l )
                        sufix <- paste(rep('_', (acc.res + flank)-l), collapse='')
                      
                      seqwin[i] <- paste(protseq[max(1, acc.res-flank):min(length(protseq), acc.res + flank)], collapse='')
                      seqwin[i] <- glue(prefix, seqwin[i], sufix) 
                  }
              } else {
                seqwin <- ''
              }
              
              rid.seqwin[acc] <- paste(seqwin, collapse='|')
          }
          
          ## add to rdesc
          rdesc <- data.frame(rdesc, VMsiteFlanks=rid.seqwin )
          
          
          ##############################################
          ## export
          gct.out <- new('GCT')
          gct.out@mat <- data.matrix(mat)
          gct.out@rdesc <- rdesc
          gct.out@cdesc <- gct.in@cdesc
          gct.out@cid <- gct.in@cid
          gct.out@rid <- rid
          
          # fn <- sub('.*/','', gct.in.str)
          # if(grepl('_n[0-9]*x[0-9]', fn)){
          #   fn <-  sub('_n[0-9]*x[0-9]*.*', '_VMsiteFlanks', fn)
          # } else {
          #   fn <-  sub('\\.gct$', '_VMsiteFlanks', fn)
          # }
          if(write.output){
            write.gct(gct.out, ofile = flank.out.dir, appenddim = F)
          }
          
          # wd <- getwd()
          # 
          # setwd(tmp.dir)
          # write.gct(gct.out, ofile = fn, appenddim = F)
          # fn.tmp <- dir('.', pattern = fn)
          # file.copy(fn.tmp, wd, overwrite = T)
          # cat("writting to: ",wd)
          # unlink(fn.tmp)
          # setwd(wd)

}
