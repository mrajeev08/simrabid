
# simrabid

`simrabid` is an implementation of a spatially explicit, individual-based model of canine rabies.

It is a work in progress, and eventually will be an R package with vignettes and documentation! But for now to play with the 
development version: clone the repo and from within the directory, run `devtools::load_all()` in R or Rstudio.

## Dependencies

- Just [`data.table`](https://rdatatable.gitlab.io/data.table/) so far!
- And devtools in the short term for dev.

## Roadmap

- A relative probability based movement model (i.e. movement probability is relative 
  to index location, for instance if you want to account for landscape barriers (roads, rivers, etc.)); 
  see old commit code for how this might work with a probability list;

- Use non utm coordinate space, instead use cellFromXY in raster and haversine distances 
 to get location of cell id movements so that you can simulate across larger 
 spatial scales; one issue is then your scale of aggregation gets distorted, aggregate grid cells somehow
 so it's approximately 1 km?

- Implement carrying capacity on pop growth, and replacement of dogs removed due to infection.

- Implement colonization of uninhabited patches (with some limits so households can't pop up in space
  that is uninhabitable, i.e. rivers/roads/etc.)

- Build a constructor class that gets passed to simrabid function so that only **valid** combinations
  of model arguments can be passed, and you only have to test this once for N simulations. 

- Profile and speed up! 
  - Vaccination function
  - Filtering data.table
  - Only track currently infectious + exposed linelist
  - Easy fixes (i.e. storing things in the appropriate type, keys onf filters, etc.)

- Issue with intermediate scales where small admin units do not get matched to 
  any grid cell

- Village metapopulation (separate function?)

- Example output environment for customizing summary functions

- Applying mortality to exposed class?

- Documentation on how to use & customize

- Construct synthetic populations or use high res pop data to get estimates of spatial dog pops (popcompr?)

- Benchmarks across scales, etc.
## Known limitations / potential extensions
- Doesn't simulate expansion of occupied cells (i.e. colonization of patches by doggos)
- Carrying capacity for growth

