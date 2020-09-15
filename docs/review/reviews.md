# SwipProteomics eLife Reviewer comments

Please provide more information about the experimental design and data analysis
to address the following reviewer comments and concerns:

> Why was WASH1 used as the fusion target for BioID? Is BioID2 alone a good
> control for this experiment? Had the authors considered using WASH1 with its
> N-terminal WASH assembly domain removed as a better control?

I think the suggestion to use a mutant WASH with its N-terminal domain removed
is really interesting. Jia et al., 2010 found that WASH1's N-terminus (aa 1-51)
are required for interactions with Fam21, SWIP, Strumpellin, CCDC53 and CAPZa.
So, you'd predict that WASH1-N51-BiRA would localize elsewhere and may provide 
a control for WASH1 non-specific binding. However, the unintendend consequences
of overexpressing this mutant protein are undesirable as a control.
BioID is a suitable control as it (more uniformly) controls for BirA expression 
throughtout the cell.

> Why was SPS-MS3 not used as the mass spectrometry approach? More accurate
> data may have been achieved even with the application of FAIMS.
Pretty sure it is SPS-MS3. Clarify with core.
```Greg

From Greg:

We did not use MS3 since we now have FAIMS to add an additional level of
specificity.  The reporter ions can be quantitated via MS2. (The additional MS
level in MS3 v. MS2 is to further isolate the peptide from interference before
reporter ion quant).

Here is an explanation from the PD 2.4 manual on reporter ion quant:  We used
the default S/N method.

"Reporter Abundance Based On: For reporter ion quantification, the PSM abundance
can be based on either S/N or peak intensity. The Automatic setting bases the
abundance on S/N for quantification spectra acquired using Orbitrap MS and bases
the abundance on intensity for quantification spectra acquired suing ion trap
MS."
```

> There are currently no plots of the data and the authors use a bespoke data
> analysis pipeline. It is not clear whether the experiment was successful and
> the results crucially depend on successful separation of the organelles. The
> paper is currently missing western blots to confirm separation of the
> organelles, upon which the analysis pipeline depends.
* To add: Boxplot showing normalization of samples
* We can show good seperation of the Fractions by PCA.

Previously, we omitted the protein PCA that was in the supplement because our fear that it 
may point towards overclustering. This plot would also go towards showing
separation of the 'organelles' == clusters. 
Looking back at this plot, yes I think we overclustered things.

> To the reader the power of using a spatial proteomics approach such as the
> LOPIT-DC method of Geladaki et al seems to be a bit lost here. The authors
> need to be clear what added information this experimental design gives over
> simply just looking at changes in the total proteome. It is hard to reconcile
> the data with organelle re-localization without firstly showing total
> abundance changes between the MUT and wildtype, and secondly giving some
> indication that the necessary organelle resolution has been achieved (as
> above). This should be clarified in the text. For example, if a protein
> resides in more than one location, and there is a change in the abundance of
> this protein in one of these locations, this may manifest itself and a
> perceived change is localization.
RE clarify the text: we can work on that.


> Furthermore, if endosomal properties are changed as a result of the mutation
> it might result in the endosomal protein pelting at different speeds. This
> will complicate the statistical analysis because fractions cannot be compared
> like for like.
If the fundamental properties of the properties of endosomes was altered this
would change things. We assume that things are globally unaffected.

> Was the final cytosolic (supernatant) fraction was discarded? If the
> abundance of the cytosolic pool of these proteins is changing then the
> observed results here could simply be due to that, rather than any
> interpreted changes at the endosome.
Cytosolic pool up/down --> OBS endosome changes

> The Authors claim (line 451) 'In addition to highlighting the neuronal roles
> of WASH in CCC- and Retriever mediated endosomal sorting, our proteomics
> approach also identified protein modules 21 with increased abundance in SWIP
> P1019R mutant brain.' This is very confusing as the protein abundance change
> are not shown, what is known is that the more of the proteins is likely to be
> associated with a complex, but not the overall abundance of the complex
> components in the cell. These arguments persist throughout the discussion
> section. If the authors take only a proportion of each fraction and this is
> not consistent across fractions or WT and MUT, then they do not know what the
> total abundance changes are.
I keep reading this thinking they are onto a good comment, but then am confused
by the time I finish the sentence.

> The authors should confirm that the endosomal enriched fraction is the
>  same in both WT and mutation experiments. 
Which fraction is 'endosome enriched'?

> In their proteome data, the authors argue that 37 out of 255 modules exhibit
> significant differences in WT and MUT brains. These data indicates that in
> addition to endo-lysosomal modules, many other pathways are also affected in
> the MUT brains. These include endoplasmic reticulum (ER) module (M83),
> synaptic modules (M35 and M248) and many others that the authors did not
> specify...did the authors observe any defects in other organelles or cellular
> compartments, such as ER, mitochondria, synapse...etc?
Point to discussion of modules that are not changing in supplement.
Review other significant modules. The reason we did not go over all significant
modules is because not all of them looked biologically tractable--that is
some looked better than others. most of the modules that exhibited change
but are not discussed have no strong enrichment for known biological processes.

> Finally, in their TEM images in figure 5, the authors argue that the
> electrical-dense inclusions in the cell bodies of MUT neurons are "visually"
> consistent with lipofusine accumulation. The authors need to use biochemical
> or histological methods to prove their point. This will significantly
> strengthen their arguments.  

I thought the TEM images of lipofuscin. Jamie wrote 'visually consistent' to be
conservative. These look like lipofuscin and this is the fields standard.
>
> It would be beneficial if the
> authors could do some IPs with the WT and mutant SWIP vectors to validate the
> proteomic data
IP what?
IP WASH to show complex is dysrupted = done. We don't have strong evidence of other
protein complexes that are changing.

Maybe they mean fractionation + western blot. Jamie will add this. Tyler will
add sample PCA. These could go in the same supplemental figure.

> The analysis method used was chosen as previous approaches to deal with
> spatial proteomics data in the literature make use of well curated organelle
> markers. The authors claim that they did not have access to a robust set of
> marker lists, but other studies have used mouse neuronal cell lines (Itzhak
> DN et al Cell Reports 2017) and also mouse ES cells (Christoforou. A et al,
> Nature Commun. 2016). These lists could easily have been adapted and used to
> visualize organelle separation using straightforward approaches such as PCA.
 This is easy: take their lists and label PCA with it--I expect it wont look
 very good. I can do this. But then I have to do it with out data. Our plot
 shows that we overclustered.

Indeed, this is what I first did. Truthfully we just didn't get the pRoloc
workflow to work. 

Need to think back to when we tried to make modules by semi-manual curration.
Ultimately, I favor a data driven approach.

Include discussion of two PATHS

> The analysis of the (spatial) proteomics data is currently not clear and
> there is some confusion. Firstly, edgeR was originally developed to handle
> RNA-sequencing data, not scRNA-sequencing data. Furthermore, RNA-sequencing
> data are indeed interpreted as counts and a negative binomial distribution is
> appropriate. This is not the case for proteomics data, as an integral under
> the isotopic envelope is involved in computing the intensity. Thus the
> analysis is not appropriate for the task. LIMMA, DEP, MSqRob, DeqMS,
> MSstatsTMT would all be appropriate methods.
Boils down to: ARE THE DATA COUNTS?

I thought it was counts:
integral under isotopic envelop = count of number of reporter ions in a time
window?
I think we are SPS3 and data is count of number of reporter ions

ADD pvalue plot:
A recommended way to evaluate
statistical model appropriateness in differential expression
omics experiments (where most genes/proteins do not have
any differential expression) is to plot the distribution of raw p
values (not corrected for multiple testing) (40). The raw p value
distributions (supplemental Fig. S5) had the expected appearance
of a flat p value distribution (for unchanged protein
expression) with modest “spikes” at small p values for the
differentially expressed proteins. We identified significantly
changing proteins in all experimental groups:

I have used DEP and limma and maybe MSstatsTMT for MRM.

I will work on making it more clear.

> The GLM framework for differential protein abundance between modules is not
> quite clear and the analysis is not quite correct. Instead of summarizing a
> module as the sum of the proteins, linear models should be fit on the data
> directly with a global module term and a factor for each protein. The protein
> factor will probably need to be encoded as a random effect. Lme4 and gam
> packages in R should be able to do this analysis. This section would gain a
> lot of clarity from some more precise descriptions.

The reviewer is correct. linear models could be fit to the data with a module
term and a factor for protein. This would be more appropriate. I can try this.
Practicallity was a motivating factor here.

> From the figure it looks like the spatial proteomics data was normalized so
> that the max intensity in the most intense fraction is 1 - is that the case?
I need to look back, but I think this is correct.

> Usually spatial proteomics data are normalized so that protein intensity sums
> to 1 across the fraction. I also find all the normalization and filtering for
> the TMT analysis quite confusing - a table might help with the desired effect
> in a column. Why did the authors not summarize peptides to proteins via the
> median or sum and then normalize so that proteins sum to 1 across the
> fractions?
Yes, this is different from other spatial proteomics. 
I can add a table which describes what we did.
What the reviewer describes is the proloc norm appraoch often used in spatial 
proteomics. because we did not do proloc i did not do it this way.

> There are no plots of the finally normalized spatial proteomics data to see
> whether the experiment was successful or not.
add sample level boxplots and pca

> The connection with human findings is a bit overstated. The findings suggest
> that the clinical phenotypes between humans and mice are similar. However,
> the mechanistic insights are only shown in mouse models. The text is slightly
> overstated and the mechanistic insights in humans should be toned down.
agree

> In their proteome data, the authors argue that 37 out of 255 modules exhibit
> significant differences in WT and MUT brains. These data indicates that in
> addition to endo-lysosomal modules, many other pathways are also affected in
> the MUT brains. These include endoplasmic reticulum (ER) module (M83),
> synaptic modules (M35 and M248) and many others that the authors did not
> specify. A major concern is that whether the endo-lysosomal dysfunction is
> the only factor that contributes to the behavioral defects? A rescue
> experiment can solve most of this concern. It has been shown recently the
> R33, a retromer chaperone, can strengthen retromer function and improves
> memory in a mouse model of AD (PMID: 31964406). The authors can consider
> testing this drug in their model.  
>
> Reviewer 3 • The authors suggest that the
> WASH complex may not interact as closely with retromer as it does in other
> cells. This is a bold statement to make based on BioID and given the existing
> literature associated with the retromer-WASH axis. For example, the
> VPS35-D620N disruption of binding to FAM21. Could the authors expand on this?
> Is it possible that the retromer complex is not present in the
> proximity-based proteomics due to the use of WASHC1?
It is a bold statement. It is in the discussion right? can reword to be less
strong, but I think its a good hypothesis.

> Could the authors expand more on their result showing that many of the lysosomal
> protein interactors are enriched in the SWIP mutant condition compared to the WT
> when many of these proteins have been shown to be lost in neurodegenerative
> disease? Do the authors think that if they looked at longer aged animal they
> would see a drop as the lysosomes become impaired and that their model is
> looking at how the cells try to compensate for the endosomal dysfunction (ie
> early stages of neurodegeneration)?

> It is interesting that the SWIP(P1019R) mutant mice exhibit such significant
> progressive motor deficits. The authors found no difference in the cleaved
> caspase 3 staining in the striatum but did
> they look at whether there was a loss of dopaminergic neurons in the substantia
> nigra pars compacta (or a loss of dopaminergic innervation or dopamine levels in
> the striatum) to account for these motor deficits? I would expect there to be a
> drastic loss of dopamine due to the significant motor deficits shown.
> Interestingly, SNCA is also present as an interactor of the WASHC1. Could the
> authors expand on whether they think alpha synuclein could therefore, also be
> playing a role (particularly as the authors also suggest an elevation of ER
> stress modulators in the SWIP mutant mice proteomics?)
good thoughts, need to look at the data again
re domapine: not sure, need to learn more about this

> [Introduction Lines 42-44] I do not quite understand the first sentence - is
> "throughout their elaborate processes" needed?  The first paragraph of the
> introduction lacks citations.  

need  to reference text

> 86-103, it would good to mention highlight that other spatial proteomics
> approaches have successfully explore trafficking pathways:

Agreed, need to elaborate on how our approach is different, and why, and 
reference previous spatial proteomics:

* AP-4 vesicles contribute to spatial control of autophagy via RUSC-dependent
  peripheral delivery of ATG9A (Davies, AK et al.Nature Comm. 9 (1) 3958)

* Role of the AP-5 adaptor protein complex in late endosome-to-Golgi retrieval
  (Hirst et al., PLOS Biology 16(1) e2004411)

* Determining the content of vesicles captured by golgin tethers using LOPIT-DC
  (Shin et al. bioRxiv 841965; doi: https://doi.org/10.1101/841965)

> Results Line 132, Benjamini-Hochberg correction is missing a citation

Need to add the following citation:

* Benjamini, Y. and Hochberg, Y., 1995. Controlling the false discovery rate: a
   practical and powerful approach to multiple testing. Journal of the Royal
   statistical society: series B (Methodological), 57(1), pp.289-300.  

> Line 143-144, the authors might find it useful that SNX1, SNX3 are regulators of
> golgin-97-vesicles destined to the trans-Golgi from the endosomes
> Determining the content of vesicles captured by golgin tethers using LOPIT-DC
> (Shin et al. bioRxiv 841965; doi: https://doi.org/10.1101/841965 ) 
Interesting

> Lines 211, GLMs are missing a citation 
add reference

> Figure 2 is quite difficult to follow, some of
> these figures could be moved to the supplement to focus on the key points.
> Since the authors have included graphs, as in nodes with edges, and graphs as
> in barcharts it would be good to be more precise in the latter case.

clarify nomenclature: 
* graph: matrix representation of a collection of nodes/edges
* figure: a bar chart or any plot
* graph: the visaul represeantion of a network as done e.g. Cytoscaope


> Figure
> 4 F it looks like a t-test was used but this is count data so should be
> corrected with the appropriate test.  Figure 4 J and K IT is not clear which
> test was applied here. The t-test is inappropriate for frequencies and the
> Wilcoxon-Mann-Whitney or chi-squared test is more appropriate.  Figure 5 F It
>  is not clear what the Kruskal-Wallis test is referring to.
need to reference text

> In the section TMT Proteomics Quantitative Analysis, it would be good to
> describe a typical spatial proteomics workflow so that readers can easily
> identify the difference. The current exposition assumes that people are familiar
> with spatial proteomics.  The authors then might consider a standard workflow

Agreed

* Breckels, L.M., Mulvey, C.M., Lilley, K.S. and Gatto, L., 2016. A
  Bioconductor workflow for processing and analysing spatial proteomics data.
  F1000Research, 5.  Then some more bespoke methods and reference might be

* Itzhak, D.N., Tyanova, S., Cox, J. and Borner, G.H., 2016. Global,
  quantitative and dynamic mapping of protein subcellular localization. Elife,
  5, p.e16950.

* Beltran, P.M.J., Mathias, R.A. and Cristea, I.M., 2016. A portrait of the
   human organelle proteome in space and time during cytomegalovirus infection.
   Cell systems, 3(4), pp.361-373.

* Itzhak, D.N., Davies, C., Tyanova, S., Mishra, A., Williamson, J., Antrobus,
  R., Cox, J., Weekes, M.P. and Borner, G.H., 2017. A mass spectrometry-based
  approach for mapping protein subcellular localization reveals the spatial
  proteome of mouse primary neurons. Cell reports, 20(11), pp.2706-2718.

* Crook, O.M., Mulvey, C.M., Kirk, P.D., Lilley, K.S. and Gatto, L., 2018. A
  Bayesian mixture modelling approach for spatial proteomics. PLoS
  computational biology, 14(11), p.e1006516.

* Crook, O.M., Breckels, L.M., Lilley, K.S., Kirk, P.D. and Gatto, L., 2019. A
  Bioconductor workflow for the Bayesian analysis of spatial proteomics.
  F1000Research, 8.

> In Figure 2B is the data normalised to the protein loading control. 
> I could not tell from the description.
Jamie


> It would be beneficial to clarify the different subcellular fractions 
yes

> Discussion has a typo
fix
