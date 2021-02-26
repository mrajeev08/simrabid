# get month from start date
get_timestep <- function(date, origin_date = "01-01-2002",
                         date_fun = lubridate::dmy,
                         units = "weeks") {

  as.numeric(lubridate::as.duration(date_fun(date) - date_fun(origin_date)), units)
}


# getting cellid of raster based on utm coords
# so as not have to match names!
get_ids <- function(x_coord, y_coord, rast, id_col) {

  cell_id <- cellFromXY(rast, cbind(x_coord, y_coord))

  admin_id <- rast[][cell_id]
  admin_code <- id_col[admin_id]

  # Use match nearest here to match to closest non-NA vals!
  if(sum(is.na(admin_id)) > 0) {
    out_rast <- match_nearest(cell_ids = cell_id[is.na(admin_id)],
                              to_match = rast,
                              max_adjacent = 10)
    admin_id <- out_rast[][cell_id]
    admin_code <- id_col[admin_id]

    if(sum(is.na(admin_id)) > 0) {
      warning(paste0(sum(is.na(admin_id)), " points went unmatched to an admin unit!"))
    }
  }

  return(admin_code)

}


#' Match to closest non-NA cell
#'
#' Matches cells in pop raster which do not have a match to a shapefile feature
#' to nearest cell that does have a match. Mainly issue with coastal
#' populations or shapefiles with holes.
#'
#' @param cell_ids cells with NA values for shapefile feature
#' @param to_match raster with associated values from shapefile
#' @param max_adjacent max cell window to look for nearest non-NA cells
#'
#' @import raster
#' @keywords internal
#'
match_nearest <- function(cell_ids, to_match, max_adjacent = 10) {

  # For each cell_id to match: get the adjacent cell_ids
  find <- data.table(adjacent(to_match, cell_ids))

  # Get the admin values at the adjacent cell_ids
  find$match <- to_match[find$to]
  suppressWarnings({
    matched <- find[, .(match = min(match, na.rm = TRUE)), by = "from"] # match is the friction grid index
  })
  matched <- matched[!is.infinite(match)]
  find <- find[!(from %in% matched$from)]


  for (i in 1:max_adjacent) {
    if (nrow(find > 0)) {
      adj_next <- adjacent(to_match, unique(find$to))
      adj_next <- data.table(to = adj_next[, "from"],
                             to_i = adj_next[, "to"]) # next to_i
      adj_next$match <- to_match[adj_next$to_i]
      # this throws warnings: throw out now!
      suppressWarnings({
        matches <- adj_next[, .(match = min(match, na.rm = TRUE)), by = "to"]
      })
      find <- find[, c("from", "to")][matches, on = "to"]
      matches <- find[, .(match = min(match, na.rm = TRUE)), by = "from"]
      matched <- rbind(matched, matches[!is.infinite(match)])

      find <- find[, c("from", "to")][adj_next, on = "to", allow.cartesian = TRUE]
      find <- unique(find[, c("to", "to_i") := .(to_i, NULL)][!(from %in% matched$from)]) # search out
    } else {
      break
    }
  }

  if (nrow(find > 0)) {
    print("Warning: still unmatched pixels!")
  }

  to_match[cell_ids] <- matched$match[match(cell_ids, matched$from)]

  return(to_match)
}

# get latest fun
get_latest <- function(path, pattern) {
  list.files(path, full.names = TRUE)[grep(pattern, list.files(path))][1]
}

