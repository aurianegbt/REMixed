## Comments regarding first releases

For warning : "You have examples for unexported functions. Please either omit these
examples or export these functions.
Examples for unexported function
   model.pasin() in:
      dynFUN_demo.Rd
      model.clairon.Rd
      model.pasin.Rd
      model.pk.Rd"
      
model.pasin is an exported function, i am not sure to understand which function is concerned. 

For warning : "\dontrun{} should only be used if the example really cannot be executed
(e.g. because of missing additional software, missing API keys, ...) by
the user. That's why wrapping examples in \dontrun{} adds the comment
("# Not run:") as a warning for the user. Does not seem necessary.
Please replace \dontrun with \donttest.
Please unwrap the examples if they are executable in < 5 sec, or replace
dontrun{} with \donttest{}.
For more details:
<https://contributor.r-project.org/cran-cookbook/general_issues.html#structuring-of-examples> "

It is the case, as this package rely heavily on Monolix software, linked to R through connectors but that does need license key. Function that could be tested/run whithout Monolix do not use "\dontrun{}". 

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
