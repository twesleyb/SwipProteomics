# MSstatsTMT code to perform protein-wise linear mixed modeling

## lmerTest::lmer Models fit by MSstatsTMT:
1. `lmer(Abundance ~ 1 + (1|Mixture) + (1|Mixture:TechRepMixture) + Group, data)`

2. `lmerAbundance ~ 1 + (1|Run) + Group, data)`

3. `lmer(Abundance ~ 1 + Group, data)`

4. `lmer(Abundance ~ 1 + (1|Mixture) + (1|Mixture:TechRepMix) + Group + (1|Subject:Group:Mixture), data)`

5. `lmer(Abundance ~ 1 + (1|Run) + Group + (1|Subjec:Group), data)`

#### lmer Variables
> Mixture: Concatenation of TMT samples (a TMT experiment)
> TechRepMixture: technical replicates of a mixture
> Run: From annotation_dt, should match Spectrum.File in raw.pd. Therefore
>    this cooresponds to the n Spectrum.Files in raw.pd. A MS analysis.
> Group: ???
> Subjec: ???

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
