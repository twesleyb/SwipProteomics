# MSstatsTMT code to perform protein-wise linear mixed modeling

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
  fit <- fit_reduced_model_onerun(sub_data)  # [3]
}

} else {

if (sub_TechReplicate | sub_bioMixture) { 
  # if there are multiple technical replicates or multiple mixtures
  # then fit the reduced model with only run effect:
  # Abundance ~ 1 + (1|Run) + Group
  fit <- fit_reduced_model_mulrun(sub_data) # [2]

  # the second model is not applicable
  if (is.null(fit)){ 

    # then fit one-way ANOVA model: 
    # Abundance ~ 1 + Group
    fit <- fit_reduced_model_onerun(sub_data) # [3]

  }

} else { 

	# else single run case
	# then fit one-way ANOVA model:
	# Abundance ~ 1 + Group
	fit <- fit_reduced_model_onerun(sub_data) # [3]
}
}

} else { 
    # else biological variation exists within each condition and mixture

    # if there are multiple biological mixtures
    if (sub_bioMixture) {  

	    # and if there are multiple technical replicate MS runs
	    if (sub_TechReplicate) {           

	      # then fit the full model with mixture, techrep, subject effects:
	      # Abundance ~ 1 + (1|Mixture) + (1|Mixture:TechRepMix)
	      # + Group + (1|Subject:Group:Mixture), data)
	      fit <- fit_full_model(sub_data) # [4]

  # if the full model is not applicable
  if (is.null(fit)) { 

    # then fit the reduced model with run and subject effects
    # Abundance ~ 1 + (1|Run) + Group + (1|Subjec:Group)
    fit <- fit_reduced_model_techrep(sub_data)  # [5]
  }

  # if the second model is not applicable
  if (is.null(fit)) { 
    # then fit one-way ANOVA model:
    # Abundance ~ 1 + Group
    fit <- fit_reduced_model_onerun(sub_data) # [3]
  }

} else { 
	# else single technical replicate MS run

	# fit the reduced model with only run effect:
	fit <- fit_reduced_model_mulrun(sub_data) # [2]

  # if the second model is not applicable
  if(is.null(fit)) { 

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
	      if(is.null(fit)){

		      # then fit one-way ANOVA model:
		      fit <- fit_reduced_model_onerun(sub_data)  # [3]
  }

# else single run
} else { 

  # fit one-way ANOVA model:
  fit <- fit_reduced_model_onerun(sub_data) # [3]

} # single technical replicate MS run
} # single biological mixture
} # biological variation

## estimate variance and df from linear models
# ...
