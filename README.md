## Filters

The filtering steps follow the general unmatched variant calling approach as employed in previous papers. In brief, they involve recounting all variants found in one patient across all samples from that patient, after which the germline variants are filtered out using an exact binomial test (germline_exact_binom.R). After that, further artefacts are filtered out using a filter based on a beta-binomial overdispersion estimation (beta_binom_filter.R). More information on these filters can be found at https://github.com/TimCoorens/Unmatched_NormSeq. The script filtering_placenta.R has the general code and workflow for the filtering of SNVs and indels.


