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

# subset data
ethics_data <- data %>% 
  filter(str_detect(abstract, "ethic|security|safe|dilemma"))

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
