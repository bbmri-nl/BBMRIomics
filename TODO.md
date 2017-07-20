Here are a few things that were not high on the priority list but
could be beneficial for the *BBMRIomics*-package:

### Adding unit-test ###

R has several ways to incorporate unit-tests to packages which can
even be added to GitHub CI hook on pushing. 

We could think of checking existence of data files and the
accessibility of the database urls.

### Using couchapp for interaction with the metadatabase ###

The current use of couchDBs' functionality is very limited this could
be extended to leverage the stability and ease of use of the
metadatabase. For example, using the
[couchapp python module](https://github.com/couchapp/couchapp). One
useful extentions would be using so-called `list`-functions that can
transform JSON to csv and no R-parsing is required. See the couchapp
_design view already available in the metadatabase for beta-testing.


There are currently too many views and not all information is directly
useful. Reducing the number of views and reevaluating the information
returned could be useful.

### Improve/Extend Molgenis R-API ###

The molgenis R-API currently lacks proper conversion of data-types and
does not allow extraction of complete tables. These limitations are
currently filled-in by *BBMRIomics*-functionality but ideally these
are ported to the molgenis R-API. 

Furthermore, currently, the R-API is available as a script using
`global`-variables to keep information for the current database
session, i.e. the user-token. This is not good pratice. A better
implementation in a separate molgenis R package should be easy to
create.
