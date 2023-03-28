library(tidyverse)

aggregate_network.dendrogram = function(m, method = 'binary') {
  if (method == 'binary') {
    distance_matrix = dist(m, method = method)
    dendrogram = hclust(distance_matrix, method = "average")
  } else if (method == 'rege') {
    dendrogram = m %>%
      as.matrix() %>%
      blockmodeling::REGE.for() %>%
      pluck('E') %>%
      as.dist() %>%
      hclust(method = "average")
    dendrogram$height = round(dendrogram$height, 6)
  }
  
  return(dendrogram)
}

aggregate_network = function(m,
                             interaction_cutoff,
                             clustering_cutoff = NULL,
                             nr_of_clusters = NULL,
                             method = 'binary') {
  nr_of_original_species = length(m)
  
  ### --- CREATE CLUSTERS --- ###
  
  if (method == 'binary') {
    distance_matrix = dist(m, method = method)
    dendrogram = hclust(distance_matrix, method = "average")
  } else if (method == 'rege') {
    dendrogram = m %>%
      as.matrix() %>%
      blockmodeling::REGE.for() %>%
      pluck('E') %>%
      as.dist() %>%
      hclust(method = "average")
    dendrogram$height = round(dendrogram$height, 6)
  }
  
  species_list = data.frame(cutree(dendrogram, k = nr_of_clusters, h = clustering_cutoff))
  species_list = species_list %>%
    rename(membership = 1) %>%
    rownames_to_column('species')
  
  cluster_list = species_list %>%
    count(membership, name = "size") %>%
    rename(cluster = membership)
  
  ### --- CONNECT CLUSTERS  --- ###
  
  nr_of_clusters = max(species_list$membership)
  
  
  possible_connections = find_possible_connections(cluster_list, nr_of_clusters)
  
  realised_connections = find_realised_connections(nr_of_original_species,
                                                   nr_of_clusters,
                                                   m,
                                                   species_list)
  ratio_possible_realised = realised_connections / possible_connections
  final_web = decide_which_clusters_interact(ratio_possible_realised,
                                             nr_of_clusters,
                                             interaction_cutoff)
  
  graph_and_species_list = list(
    graph = as.matrix(final_web),
    species_list = species_list %>% arrange(membership) %>% relocate(membership)
  )
  
  return(graph_and_species_list)
}
