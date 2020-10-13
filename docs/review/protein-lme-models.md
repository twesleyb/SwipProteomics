# MSstatsTMT Protein-wise linear mixed modeling



```R
# Code for repeated measure design from MSstats

} else { 
    if (singleSubject) {
	fit.full <- lm(ABUNDANCE ~ GROUP, data=data2) # [1]

    } else { 
	if (!TechReplicate) {
	      fit.full <- lmer(ABUNDANCE ~ GROUP + (1|SUBJECT), data=data2) # [2]
	      df.full <- lm(ABUNDANCE ~ GROUP + SUBJECT, data=data2)$df.residual

	} else {
	      fit.full <- lmer(ABUNDANCE ~ GROUP + (1|SUBJECT) + (1|GROUP:SUBJECT), data=data2) # [3]
	      df.full <- lm(ABUNDANCE ~ GROUP + SUBJECT + GROUP:SUBJECT,  data=data2)$df.residual
	}
    }	
}
```
There are three linear models:
1. `lm(ABUNDANCE ~ GROUP, data)`
2. `lmer(ABUNDANCE ~ GROUP + (1|SUBJECT), data)`
3. `lmer(ABUNDANCE ~ GROUP + (1|SUBJECT) + (1|GROUP:SUBJECT), data)`

> Model [1] is the most basic design. No mixed effects. Model [2] includes
> a mixed effect of Subject for repeated measures.
> Model [3] is for Technical replication?

Model __[2]__ is for repeated measures:
`fit.full <- lmer(ABUNDANCE ~ GROUP + (1|SUBJECT), data=data2) # [2]`
> mixed effect of Subject/BioReplicate

## lmerTest::lmer Models fit by MSstatsTMT:
> NOTE: Group == Condition; Subject == BioReplicate;
1. `lmer(Abundance ~ 1 + (1|Mixture) + (1|Mixture:TechRepMixture) + Group, data)`

2. `lmerAbundance ~ 1 + (1|Run) + Group, data)`

3. `lmer(Abundance ~ 1 + Group, data)`

4. `lmer(Abundance ~ 1 + (1|Mixture) + (1|Mixture:TechRepMix) + Group + (1|Subject:Group:Mixture), data)`

5. `lmer(Abundance ~ 1 + (1|Run) + Group + (1|Subject:Group), data)`

#### lmer Variables
> Mixture: Concatenation of TMT samples (a TMT experiment)
> TechRepMixture: technical replicates of a mixture
> Run: From annotation_dt, should match Spectrum.File in raw.pd. Therefore
>    this cooresponds to the n Spectrum.Files in raw.pd. A MS analysis.
> Group: Condition
> Subject: BioReplicate

## MSstatsTMT lmer wrappers:
1. `fit <- fit_full_model_spikedin(sub_data)`
2. `fit <- fit_reduced_model_mulrun(sub_data)`
3. `fit <- fit_reduced_model_onerun(sub_data)`
4. `fit <- fit_full_model(sub_data)`
5. `fit <- fit_reduced_model_techrep(sub_data)`

```R
if (sub_singleSubject) {
  # if no biological variation within each condition and mixture

  if (sub_TechReplicate & sub_bioMixture) {
    # and if multiple mixtures and technical replicates

    # fit the full model with mixture and techrep effects:
    # Abundance ~ 1 + (1|Mixture) + (1|Mixture:TechRepMixture) + Group
    fit <- fit_full_model_spikedin(sub_data) # [1]

    # if full model is not applicable
    if (is.null(fit)) {

      # then fit the reduced model with only run effect:
      # Abundance ~ 1 + (1|Run) + Group
      fit <- fit_reduced_model_mulrun(sub_data) # [2]
    }

    # if the second model is not applicable
    if (is.null(fit)) {

      # then fit one-way ANOVA model:
      # Abundance ~ 1 + Group
      fit <- fit_reduced_model_onerun(sub_data) # [3]
    }

  } else {

    # if there are multiple technical replicates or multiple mixtures
    if (sub_TechReplicate | sub_bioMixture) {

      # then fit the reduced model with only run effect:
      # Abundance ~ 1 + (1|Run) + Group
      fit <- fit_reduced_model_mulrun(sub_data) # [2]

      # if the second model is not applicable
      if (is.null(fit)) {

        # then fit one-way ANOVA model:
        # Abundance ~ 1 + Group
        fit <- fit_reduced_model_onerun(sub_data) # [3]
      }

      # else single run case
    } else {

      # fit one-way ANOVA model:
      # Abundance ~ 1 + Group
      fit <- fit_reduced_model_onerun(sub_data) # [3]
    }
  }

  # else biological variation exists within each condition and mixture
} else {

  # if there are multiple biological mixtures
  if (sub_bioMixture) {

    # and if there are multiple technical replicate MS runs
    if (sub_TechReplicate) {

      # then fit the full model with mixture, techrep, subject effects:
      # Abundance ~ 1 + (1|Mixture) + (1|Mixture:TechRepMix) + 
      #   Group + (1|Subject:Group:Mixture), data)
      fit <- fit_full_model(sub_data) # [4]

      # if the full model is not applicable
      if (is.null(fit)) {

        # then fit the reduced model with run and subject effects
        # Abundance ~ 1 + (1|Run) + Group + (1|Subjec:Group)
        fit <- fit_reduced_model_techrep(sub_data) # [5]
      }

      # if the second model is not applicable
      if (is.null(fit)) {
        # then fit one-way ANOVA model:
        # Abundance ~ 1 + Group
        fit <- fit_reduced_model_onerun(sub_data) # [3]
      }

      # else single technical replicate MS run
    } else {

      # fit the reduced model with only run effect:
      fit <- fit_reduced_model_mulrun(sub_data) # [2]

      # if the second model is not applicable
      if (is.null(fit)) {

        # then fit one-way ANOVA model:
        # Abundance ~ 1 + Group
        fit <- fit_reduced_model_onerun(sub_data) # [3]
      }
    }

    # else single biological mixture
  } else {

    # if multiple technical replicate MS runs
    if (sub_TechReplicate) {

      # then fit the reduced model with run and subject effects:
      # Abundance ~ 1 + (1|Run) + Group + (1|Subjec:Group)
      fit <- fit_reduced_model_techrep(sub_data) # [5]

      # if second model is not applicable
      if (is.null(fit)) {

        # then fit one-way ANOVA model:
        fit <- fit_reduced_model_onerun(sub_data) # [3]
      }

      # else single run
    } else {

      # fit one-way ANOVA model:
      fit <- fit_reduced_model_onerun(sub_data) # [3]

    } # single technical replicate MS run
  } # single biological mixture
} # biological variation

## estimate variance and df from linear models
# ... code goes on ...
```
