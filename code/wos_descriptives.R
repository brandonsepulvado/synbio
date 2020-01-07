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
data <- read.csv(here('data', 'web_of_science', 'data_subset.csv'), header = TRUE)

# descriptives =================================================================

# count publications per year
data %>% 
  filter(!is.na(pub_year), 
         pub_year < 2020) %>% 
  count(pub_year, sort = TRUE) %>% 
  ggplot(aes(x = pub_year, y = n)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 1990, color = 'orange2') +
  theme_minimal() +
  labs(x = 'Publication Year',
       y = 'Number of Publications',
       title = 'Number of Publications Per Year, 1900-2019',
       subtitle = 'Source: Web of Science',
       caption = 'Note: Vertical line denotes 1990; years without publications have no point.')

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
       x = NULL,
       title = 'Top 20 Author Keywords, 1900-2019',
       subtitle = 'Source: Web of Science'
       )

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
  labs(y = 'Number of Publications',
       x = NULL,
       title = 'Most Common Fields, 1900-2019',
       subtitle = 'Source: Web of Science')

# publication title counts
data %>% 
  filter(!is.na(pub_title)) %>% 
  count(pub_title, sort = TRUE) %>% 
  top_n(20) %>% 
  ggplot(aes(x = reorder(pub_title, -desc(n)), y = n, fill = pub_title)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_minimal() +
  labs(y = 'Number of Publications',
       x = NULL,
       title = 'Most Common Publication Outlets, 1900-2019',
       subtitle = 'Source: Web of Science')

# ethics-related records =======================================================

# subset data
ethics_data <- data %>% 
  filter(str_detect(abstract, "ethic|security|safe|dilemma"))

# count publications per year
ethics_data %>% 
  filter(!is.na(pub_year), 
         pub_year < 2020) %>% 
  count(pub_year, sort = TRUE) %>% 
  ggplot(aes(x = pub_year, y = n)) +
  geom_point() +
  geom_line() +
  geom_vline(xintercept = 2007, color = 'orange2') +
  theme_minimal() +
  labs(x = 'Publication Year',
       y = 'Number of Publications',
       title = 'Number of Ethics-Related Publications Per Year, 1900-2019',
       subtitle = 'Source: Web of Science',
       caption = 'Note: Vertical line denotes 2007; years without publications have no point; These publications are a subset of the previous records \n whose abstract contains ethic*, security, safe*, or dilemma..')


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
  labs(y = 'Number of Publications',
       x = NULL,
       title = 'Top Keywords in Ethics-Related Publications',
       subtitle = 'Source: Web of Science',
       caption = 'Note: These publications are a subset of the previous records whose abstract contains ethic*, security, safe*, or dilemma.') 

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
       x = NULL,
       title = 'Most Common Fields in Ethics-Related Publications',
       subtitle = 'Source: Web of Science',
       caption = 'Note: These publications are a subset of the previous records whose abstract contains ethic*, security, safe*, or dilemma.')

# publication title counts
ethics_data %>% 
  filter(!is.na(pub_title)) %>% 
  count(pub_title, sort = TRUE) %>% 
  top_n(20) %>% 
  ggplot(aes(x = reorder(pub_title, -desc(n)), y = n, fill = pub_title)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_minimal() +
  labs(y = 'Number of Publications',
       x = NULL,
       title = 'Most Common Outlets in Ethics-Related Publications',
       subtitle = 'Source: Web of Science',
       caption = 'Note: These publications are a subset of the previous records whose abstract contains ethic*, security, safe*, or dilemma.')

# distribiution of title counts
# TODO: figure out how to limit number of titles with same counts
ethics_data %>% 
  filter(!is.na(pub_title), !is.na(pub_year)) %>% 
  mutate(period_indic = case_when(pub_year <= 2007 ~ '<= 2007',
                                  TRUE ~ '> 2007')) %>% 
  group_by(period_indic) %>% 
  count(pub_title, sort = TRUE) %>% 
  top_n(5, n) %>% 
  ggplot(aes(x = reorder(pub_title, -desc(n)), y = n, fill = pub_title)) +
  facet_wrap(~period_indic) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  theme_minimal() +
  labs(y = 'Number of Publications',
       x = NULL,
       title = 'Most Common Outlets in Ethics-Related Publications',
       subtitle = 'Source: Web of Science',
       caption = 'Note: These publications are a subset of the previous records whose abstract contains ethic*, security, safe*, or dilemma.')

