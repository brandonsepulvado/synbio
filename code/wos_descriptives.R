# ==============================================================================
# basic descriptives for wos synbio data
# ==============================================================================

# preliminaries ================================================================

# load packages
library(dplyr)
library(here)
library(stringr)
library(tidytext)
library(ggplot2)
library(forcats)
library(readxl)
library(tidyr)

# load data
data <- read.csv('/Users/brandonsepulvado/Documents/synbio/data/web_of_science/data_subset.csv',
                 header = TRUE)

# descriptives =================================================================

# count publications per year
data %>% 
  filter(!is.na(pub_year), 
         pub_year < 2020) %>% 
  count(pub_year, sort = TRUE) %>% 
  ggplot(aes(x = pub_year, y = n)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(x = 'Publication Year',
       y = 'Number of Publications')

# top keywords
data %>% 
  filter(!is.na(author_keywords)) %>% 
  unnest_tokens(keyword, author_keywords, token = "regex", pattern = ";") %>%
  mutate(keyword = str_trim(keyword, side = 'both')) %>% 
  count(keyword, sort = TRUE) %>% 
  top_n(20) %>% 
  ggplot(aes(x = reorder(keyword, -desc(n)), y = n, fill = keyword)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_minimal() +
  labs(y = 'Number of Occurrences',
       x = NULL)

# open access versions
data %>% 
  mutate(open_access = as.character(open_access),
         open_access = case_when(is.na(open_access) ~ "Not OA",
                                 open_access == "" ~ "Not OA",
                                 TRUE ~ open_access)) %>% 
  count(open_access, sort = TRUE)

# count of wos categories
data %>% 
  filter(!is.na(wos_categs)) %>% 
  unnest_tokens(category, wos_categs, token = "regex", pattern = ";") %>% 
  mutate(category = str_trim(category, side = 'both')) %>% 
  count(category, sort = TRUE) %>% 
  top_n(20) %>% 
  ggplot(aes(x = reorder(category, -desc(n)), y = n, fill = category)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_minimal() +
  labs(y = 'Number of Occurrences',
       x = NULL)

# publication title counts
data %>% 
  filter(!is.na(pub_title)) %>% 
  count(pub_title, sort = TRUE) %>% 
  top_n(20) %>% 
  ggplot(aes(x = reorder(pub_title, -desc(n)), y = n, fill = pub_title)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_minimal() +
  labs(y = 'Number of Occurrences',
       x = NULL)

# ethics-related records =======================================================

# ethics data (filtered in python)
ethics_data <- read_excel('/Users/brandonsepulvado/Documents/synbio/data/web_of_science/data_processed/data_ethics.xlsx')

# # subset data
# ethics_data <- data %>% 
#   filter(str_detect(abstract, "ethic|security|safe|dilemma")) %>% 
#   as_tibble()

ethics_data %>% 
  filter(str_detect(AB, "auto")) %>% 
  View()

# count institutions per year
ethics_data %>% 
  filter(!is.na(PY), 
         PY < 2020) %>% 
  mutate(inst = str_extract_all(C1, "\\][[:alnum:][:blank:]]+,")) %>% 
  unnest(cols = c(inst)) %>% 
  mutate( 
    inst = str_remove_all(inst, "[[:punct:]]"),
    inst = str_trim(inst, "both"),
    inst = str_to_lower(inst)
    ) %>% 
  select(PY, inst) %>% 
  distinct() %>% 
  group_by(PY) %>% 
  summarise(n_inst = n_distinct(inst)) %>% 
  ggplot(aes(x = PY, y = n_inst)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = 2007, color = 'orange2') +
  theme_minimal() +
  labs(x = 'Publication Year',
       y = 'Number of Institutions',
       title = 'Number of Unique Institutions per Year, 1900-2019',
       subtitle = 'Source: Web of Science',
       caption = "Note: Vertical line is at 2007, when the number of publications began exponential growth.")


# number of unique authors per year
ethics_data %>% 
  filter(!is.na(PY), 
         PY < 2020,
         !is.na(AU)) %>% 
  unnest_tokens(author, AU, token = stringr::str_split, pattern = ";") %>% 
  mutate(author = str_trim(author, "both")) %>% 
  select(PY, author) %>% 
  distinct() %>% 
  group_by(PY) %>% 
  summarise(n_auth = n_distinct(author)) %>% 
  ggplot(aes(x = PY, y = n_auth)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = 2007, color = 'orange2') +
  theme_minimal() +
  labs(x = 'Publication Year',
       y = 'Number of Authors',
       title = 'Number of Unique Authors per Year, 1900-2019',
       subtitle = 'Source: Web of Science',
       caption = "Note: Vertical line is at 2007, when the number of publications began exponential growth.")

# count publications per year
ethics_data %>% 
  filter(!is.na(pub_year), 
         pub_year < 2020) %>% 
  count(pub_year, sort = TRUE) %>% 
  ggplot(aes(x = pub_year, y = n)) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(x = 'Publication Year',
       y = 'Number of Publications')

# top keywords around ethics
ethics_data %>% 
  filter(!is.na(author_keywords)) %>% 
  unnest_tokens(keyword, author_keywords, token = "regex", pattern = ";") %>%
  mutate(keyword = str_trim(keyword, side = 'both')) %>% 
  count(keyword, sort = TRUE) %>% 
  top_n(20) %>% 
  ggplot(aes(x = reorder(keyword, -desc(n)), y = n, fill = keyword)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_minimal() +
  labs(y = 'Number of Occurrences',
       x = NULL) 

# count of wos categories
ethics_data %>% 
  filter(!is.na(wos_categs)) %>% 
  unnest_tokens(category, wos_categs, token = "regex", pattern = ";") %>% 
  mutate(category = str_trim(category, side = 'both')) %>% 
  count(category, sort = TRUE) %>% 
  top_n(20) %>% 
  ggplot(aes(x = reorder(category, -desc(n)), y = n, fill = category)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_minimal() +
  labs(y = 'Number of Occurrences',
       x = NULL)
