# Convert duration data to binary

Converts data on the time dyads spent interacting into approximate
binary data, for use with a binary edge model. Both the time spent
interacting and the time observed are divided by `max_event_length` and
rounded, so each row records the number of periods of that length with
an event, out of the number of periods observed. This is a rough
approximation: use it with caution, and make sure all durations are in
the same units.

## Usage

``` r
convert_duration_to_binary(df, event_name, divisor_name, max_event_length)
```

## Arguments

- df:

  Original data frame, ready for use with
  [`bison_model()`](https://bisonr.org/reference/bison_model.md).

- event_name:

  Name of the column with the time spent interacting.

- divisor_name:

  Name of the column with the time observed.

- max_event_length:

  The maximum length of time a single social event could last, in the
  same units as the durations.

## Value

The input data frame with the `event_name` and `divisor_name` columns
replaced by counts suitable for a binary model, with the converted
`divisor_name` column as the effort, e.g.
`(time_together | time_observed) ~ dyad(node_1, node_2)`.

## See also

[`bison_model()`](https://bisonr.org/reference/bison_model.md)

## Examples

``` r
df <- data.frame(
  node_1 = c("A", "A", "B"), node_2 = c("B", "C", "C"),
  time_together = c(30, 0, 10), time_observed = c(60, 60, 60)
)
convert_duration_to_binary(df, "time_together", "time_observed", max_event_length = 10)
#>   node_1 node_2 time_together time_observed
#> 1      A      B             3             6
#> 2      A      C             0             6
#> 3      B      C             1             6
```
