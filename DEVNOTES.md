# Developer Notes

## Adding functions to the package

When adding new functions to the package, roxygen must also be updated, so once a new function has been added (or a function has been modified), make sure the roxygen skeleton is correct and run the following lines:
```
devtools::document()
devtools::load_all()
```

## Adding Stan models or issues with existing ones

When adding a new Stan model and sometimes when modifying an existing one, it seems to be necessary to run the following three lines to recompile the models:
```
pkgbuild::compile_dll()
roxygen2::roxygenize()
devtools::install()
```

Not running these commands can lead to strange behaviour, for example `rstan::sampling()` appeared to use an old version of the model while `rstan::vb()` used a newer version of the model.
