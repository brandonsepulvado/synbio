# ==============================================================================
# structural topic models with web of science data 
# ==============================================================================

# preliminaries ================================================================

# first load data in wos_descriptives.R

# load packages
library(stm)
library(stmCorrViz)
library(here)

# set seed
set.seed(1234)

# load data into stm object ====================================================

# ensure that there are no missing abstract
data <- data %>% 
  filter(!is.na(abstract))

# preprocess
processed <- textProcessor(documents = data$abstract,
                           metadata = data,
                           removenumbers = FALSE,
                           removepunctuation = FALSE,
                           ucp = FALSE)
out <- prepDocuments(processed$documents, processed$vocab, processed$meta)
docs <- out$documents
vocab <- out$vocab
meta <-out$meta

# check term frequency 
plotRemoved(processed$documents, lower.thresh = seq(1, 50, by = 1))
out <- prepDocuments(processed$documents, processed$vocab,processed$meta, 
                     lower.thresh = 10)

# topic number diagnostics 
k_diag <- searchK(out$documents, out$vocab, K = seq(5, 20, 1), data = meta, 
                  cores = 5L)

# save result to save time later
# saveRDS(k_diag, 
#         file = here::here('data', 'k_diag_5-20.rds'))


# plot results
plot(k_diag)

# k diagnostics with k = 0
k_diag_0 <- searchK(out$documents, out$vocab, K = 0, data = meta) # 107 topics

# ethics subset ================================================================

# preprocess
ethics_processed <- textProcessor(documents = ethics_data$abstract,
                           metadata = ethics_data,
                           removenumbers = FALSE,
                           removepunctuation = FALSE,
                           ucp = FALSE)
ethics_out <- prepDocuments(ethics_processed$documents, 
                     ethics_processed$vocab, 
                     ethics_processed$meta)
ethics_docs <- ethics_out$documents
ethics_vocab <- ethics_out$vocab
ethics_meta <- ethics_out$meta

# check term frequency 
plotRemoved(ethics_processed$documents, lower.thresh = seq(1, 50, by = 1))
ethics_out <- prepDocuments(ethics_processed$documents, 
                            ethics_processed$vocab,
                            ethics_processed$meta, 
                            lower.thresh = 5)

# topic number diagnostics 
ethics_k_diag <- searchK(out$documents, 
                         out$vocab, 
                         K = seq(5, 20, 1), 
                         data = meta, 
                         cores = 5L)
# save results
saveRDS(ethics_k_diag,
        file = here::here('data', 'ethics_k_diag.rds'))

# estimate stm on ethics data ==================================================

ethics_k10 <- stm(documents = ethics_out$documents, 
                  vocab = ethics_out$vocab,
                  K = 9, 
                  max.em.its = 100, 
                  data = ethics_out$meta,
                  init.type = "Spectral")
