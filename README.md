

# Some context

This here is a demo repository. It's a small script how to
create upper air surface wind barb maps for lectures.

The script is not very generic as the input grib files
have to be __univariate__ (i am just reading field number one,
and expect to get what I expect).

If you `source('georg.R')` or `Rscript georg.R` the first
part will be executed where some random points will be picked
to draw wind barbs.

After the `stop(...)` command there is a second part which
mainly does the same, but you can pick (`n=10`) locations
manually where to draw the map.

Adding errors is not yet implemented, but you mainly have
to modify `barbs` before creating `spbarbs` and before you
plot these guys.

                        Pre-Alpha :)

# Extensions

Remove wind barbs (if needed). Package that stuff. A more
generic read-grib version. I would suggest the `getgrib`
package I've written, which would, indeed, create a maybe
unnecessary dependency.
