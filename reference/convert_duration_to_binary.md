# Convert duration data to binary

Uses a naive method for converting duration data to binary data. Use
with caution. Ensure that all time periods are in the same units.

## Usage

``` r
convert_duration_to_binary(df, event_name, divisor_name, max_event_length)
```

## Arguments

- df:

  Original dataframe, ready for use with `bison_model`.

- event_name:

  The variable name in the dataframe corresponding to event durations.

- divisor_name:

  The variable name in the dataframe corresponding to divisor durations
  (e.g. observation time).

- max_event_length:

  The maximum length of time a single social event could take place for.

## Value

A dataframe the same as the input but with `event_name` and
`divisor_name` columns replaced with values compatible with binary
models.
