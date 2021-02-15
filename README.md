
# simrabid

`simrabid` is an implementation of a spatially explicit, individual-based model of canine rabies.

It is a work in progress, and eventually will be an R package with vignettes and documentation! But for now to play with the 
development version: clone the repo and from within the directory, run `devtools::load_all()` in R or Rstudio.

## Dependencies

- Just [`data.table`](https://rdatatable.gitlab.io/data.table/) so far!
- And devtools in the short term for dev.

## To Do

- Keys for faster data.table filtering (i.e. for boolean values especially! store the 0/1 as logical!)
- Do tests to make sure it makes sense with SD data (unit tests?!)
- Document and build! (managing namespace & dependencies)
- Test summary stats (for the longer ones write in Rcpp)
- Try them out on the cluster (and ask HPC folks if it's smarter to do it as array or mpi cluster)
- Make a branch for probability based movements (i.e. cell based movements)
- Keep both options (continuous space vs. discretized -- does it amount to
the same @ fine scales)(vignette)
- Separate branch for Rcpp functions for movement (continuous time & cell based), transmission;
- Last optimization should be moving to indexed vector world (rather than lists) (accumulators should get spit out at beginning and doubled when necessary)(exposed gets inserted as needed)(bind together as data.table at end)
- Benchmark vignette for decisions on filtering and syntax within functions (worthwhile to run again and on multiple systems maybe?) (MEMORY & TIME & ACROSS OBJECT SIZES)
- Pass everything that can be passed as vector/mat/list as such instead of data.table/frame (i.e. try vacc_dt as vector filter & use matching on loc_id to join)--vector of location ids, steps, and number vaccinated)
- It should be faster to pass args as vectors, because every $ call costs something
- For drastic optimizing then store everything as vectors with matching indices that get updated & make data.table @ end (but this shouldn't matter too much, maybe check how many calls there are?) And probably gets over written by how long it takes to do split apply vs. data.table...
- So only for ones where you don't have to do any aggregation!
- Make init more flexible in how it initializes I & E seeds
- Or seed location & tstep explicitly of starting cases?
- components to add = puppy vaccination (in between campaigns)
- and also immigration, colonization of empty patches?
- Figure out how to flexibly pass additional arguments to functions passed through (i.e. incursion functions, etc.); if you want to pass to fit, pass in param list!
- Else args
- Double check using appropriate & vs. &&, | vs. || throughout
- Limit ifelse statements
- Store things as integers where possible

## Known limitations / potential extensions
- Doesn't simulate expansion of occupied cells (i.e. colonization of patches by doggos)
- Carrying capacity for growth

