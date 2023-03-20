#' This package provides functions to apply geographical filters to biodiversity survey data.
#'
#' @docType package
#' @name geotag

library(rgbif)
library(logger)
library(taxonbridge)

#' Search for GBIF occurrences in a geographical area.
#'
#' @param tax.list A list or vector of taxon keys from the GBIF backbone.
#' @param geometry (character) Searches for occurrences inside a polygon in Well Known Text (WKT) format.
#' @param country (character) The 2-letter country code (ISO-3166-1) in which the occurrence was recorded.
#' @param min.occ (integer) The minimum number of occurrences for a taxon to be considered present in the geometry/country.
#' @example \dontrun{has.gbif.occ(c(7191157), country="fr")}
#' @import rgbif
has.gbif.occ <- function(tax.list, geometry = NULL, country = NULL, min.occ = 1) {
  occ <- occ_data(
    taxonKey = tax.list,
    geometry = geometry,
    country = country,
    hasCoordinate = TRUE,
    hasGeospatialIssue=FALSE,
    limit = 1
  )
  lapply(occ, function(x) { x$meta$count >= min.occ })
}

# TODO : if rank information available, use to resolve ambiguous names
#' Prepare a taxonomic reference by merging the NCBI taxonomy and the GBIF backbone taxonomy.
#'
#' @param data.frame
#' @param col.name
#' @param geometry
#' @param country
#' @param min.occ
#' @param taxid
#' @param ref.taxo
#' @param taxonomy
#' @example \dontrun{gbif.filter()}
#' @export gbif.filter
#' @import logger
gbif.filter <- function(data.frame, col.name, geometry=NULL, country=NULL, min.occ = 1, taxid=FALSE, ref.taxo="gbif", taxonomy=NULL) {
  ids = data.frame[,col.name]
  origin.ids <- unique(ids[!is.na(ids)])
  log_info('Found {length(origin.ids)} unique taxa')
  log_info("Map taxonomic identifiers to GBIF Backbone")
  if(taxid) {
    if (ref.taxo == "ncbi") {
      tax.list <- lapply(origin.ids, function(x) {
        matches <- taxonomy[which(taxonomy$ncbi_id == x),]
        if(nrow(matches) > 0) taxonomy[which(taxonomy$ncbi_id == x),]$taxonID else NA
      })
    }
    else {
      tax.list <- origin.ids
    }
  }
  else {
    tax.list <- lapply(origin.ids, 
      function(x) { 
        matches = taxonomy[which(taxonomy$canonicalName == x),]$taxonID
        if(length(matches) == 1) return(matches) else NA
      })
  }
  names(tax.list) <- origin.ids
  search.tax.list = tax.list[!is.na(tax.list)]
  log_info("Mapped {length(search.tax.list)} taxonomic identifiers")
  not_mapped = origin.ids[is.na(tax.list)]
  if(length(not_mapped) > 0)
    log_info("Could not map the following taxa : {list(not_mapped)}")
  log_info("Search for GBIF occurrences")
  has_occ <- has.gbif.occ(search.tax.list, geometry, country, min.occ)
  names(has_occ) <- search.tax.list #origin.ids
  log_info("Add gbif_filter flag")
  data.frame$gbif_filter <- lapply(data.frame[,col.name], 
    function(x) {
       if(x %in% names(search.tax.list) & length(search.tax.list[[as.character(x)]]) == 1) 
         has_occ[[search.tax.list[[as.character(x)]]]] 
       else NA
    }
  )
  return(data.frame)
}

#' Prepare a taxonomic reference by merging the NCBI taxonomy and the GBIF backbone taxonomy.
#'
#' @param path.to.gbif Path to a pre-downloaded version of the the GBIF backbone taxonomy.
#' @param path.to.ncbi Path to a pre-downloaded version of the the NCBI taxonomy.
#' @example \dontrun{prepare.taxonomy()}
#' @export prepare.taxonomy
#' @import taxonbridge
#' @import logger
prepare.taxonomy <- function(path.to.gbif=NULL, path.to.ncbi=NULL) {
  path.to.gbif <- if (!is.null(path.to.gbif)) path.to.gbif else download_gbif()
  path.to.ncbi <- if (!is.null(path.to.ncbi)) path.to.ncbi else download_ncbi(taxonkitpath = "/usr/local/bin/taxonkit")
  taxonomy <- load_taxonomies(path.to.gbif, path.to.ncbi)
  log_info('Saved GBIF Backbone to {path.to.gbif}')
  log_info('Saved NCBI Taxonomy to {path.to.ncbi}')
  lineages <- get_lineages(taxonomy)
  invalid.kingdom <- get_validity(lineages, rank = "kingdom", valid = FALSE)
  valid.taxonomy <- taxonomy[!taxonomy$taxonID %in% invalid.kingdom$taxonID,]
  log_info('Save taxonomy to file taxonomy.rds')
  saveRDS(valid.taxonomy, file = "taxonomy.rds")
}