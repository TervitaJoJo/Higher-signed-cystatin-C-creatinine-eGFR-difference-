suppressPackageStartupMessages(library(dplyr))
root <- commandArgs(TRUE)[1]
d <- read.csv(file.path(root,"tables/clinical_enhancement_model_results.csv"))
d <- d %>% filter(analysis=="primary_linear",term=="delta10") %>%
  mutate(transition=recode(outcome,recovery_event="independent",persistent_event="difficulty",death_event="dead"),
    outcome=recode(outcome,recovery_event="Recovery",persistent_event="Persistent difficulty",death_event="Death")) %>%
  select(cohort,transition,outcome,odds_ratio,ci_lower,ci_upper,p,n_people,n_intervals,events)
write.csv(d,file.path(root,"tables/final_competing_transition_results.csv"),row.names=FALSE)
