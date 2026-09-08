# Brain dump about my CATE simulations poster

## Context

I submitted an abstract for a poster at the International Clinical Trials Methodology Conference 2026, which was accepted. The abstract detailed simulations to explore the sample size and missing covariate handling requirements when using causal machine learning estimators to estimate conditional average treatment effects (CATEs) in RCTs. This abstract also discussed confidence intervals for CATE and HTE testing procedures. I also want to include some results I have exploring some techniques for evaluating and selecting CATE models, and some exploratory work I've done on subgroup identification and validation within a single RCT.

The poster I will present will be A0 and will use a template from Imperial College London; templates are available in both PowerPoint and LaTeX. I am leaning towards using PowerPoint due to the ease of moving elements of the poster around.

I would like to create 6 figures to illustrate headline results from each component of my simulations: sample size, missing data, HTE testing, confidence intervals, model selection/evaluation, and CATE/subgroup validation. I am open to having fewer figures, but I am struggling to figure out how best to communicate results from each element of these simulations.

The following sections have both sample text that could be included in the poster, and notes on the ideas that I would like to be clear from the poster, or at least intriguing enough to encourage conversations with poster viewers.

## Poster sections/structure

### Motivation/Background

HTEs are fundamental to personalised medicine. Recent CML advances have enabled estimation of causal quantities, including the conditional average treatment effect (CATE), though applications typically rely on large observational datasets. As RCTs inherently satisfy overlap and no unmeasured confounding assumptions, they are well-suited candidates for CML application.

However, there are several concerns that need to be addressed to facilitate wider application of CML in RCTs:
- Sample size requirements
- Uncertainty estimation
- Statistical testing
- Missing covariate handling
- Model selection
- Internal validation

This research explores requirements in relation to the above via extensive simulation, using a selection of popular CML methods.

### Methods/Simulation Details

*This bit is quite boring, so it would be nice to keep it brief and refer people to a more extensive document, which I'll probably make in the next week with the full ADEMP framework.*

#### Data Generating Mechanisms

2-arm RCTs with simple randomisation. Sample size varies across components between $n=100$ and $n=1000$. The scenarios varied across the different components. I want to focus on 4 key scenarios:
- Null scenario: no heterogeneity, only an average treatment effect
- Simple HTE: single covariate drives heterogeneity
- Multivariable HTE: multiple effect modifiers with interaction terms
- Non-linear HTE: cosine of an effect modifier

Both binary and continuous outcomes were explored.

Missing data mechanisms: MAR, where missingness in covariates is dependent on the observable data; MNAR, where missingness is dependent on unobserved covariates; and MNAR-Y, where missingness is dependent on an unobserved covariate which is also an effect modifier.

#### Methods

CML models considered:
- Causal forest, a popular CML method with a user-friendly package and a wealth of documentation. Also, it's a special case of the R learner.
- DR learner, the doubly-robust learner, which is a two-stage estimator for CATEs. The first stage requires estimation of a propensity score and a conditional marginal outcome function, which are used to construct pseudo-outcomes, noisy estimates for the CATE, which are then refined in a second-stage regression of the pseudo-outcomes $\phi$ on covariates.
  - Random forest learners
  - SuperLearner algorithm, an ensembling method for combining many models to achieve optimal prediction
  - Semi-oracle: use a fixed propensity score of 0.5 instead of estimating it
  - Oracle: benchmark where the conditional marginal outcomes are known *don't think I should actually include this in the poster because it's not realistic*

Use a selection of different crossfitting procedures depending on the model; this was explored before to figure out the best approaches, and explaining this on the poster is probably too much detail.

Three ways to test heterogeneity were considered:
- Best linear predictor: fits a linear model of the treatment effect and tests for a non-zero coefficient of the HTE component of the model.
- Permutation test: testing permutation independence between the covariates and either the pseudo-outcomes or the final CATEs. (I think the permutation test on the final CATEs could be removed from the poster because it's always overconservative by construction and not actually correct to use, I think)

Several ways of handling missing data were considered:
- Complete case analysis
- IPW
- Single imputation:
  - Mean
  - missForest
  - Parametric model (GLM)
- Missing indicators (MPA): impute missing values with the mean and create a binary flag for each variable to indicate whether it was missing for each observation
- Missingness incorporated in attributes (MIA): specific to random forests; within each tree, decide whether to split on missingness or send all missing observations to the left/right of a split
- Multiple imputation: 50 imputations using a forest-based model and Rubin's rules to recombine

Confidence intervals were either generated using the causal forest's inbuilt approach, or, for RF models, using the half-sample bootstrap procedure.

Model selection procedures were compared using the DR score (pseudo-outcome) as a surrogate for the true CATE in a PEHE metric, or using a bias-corrected version of the PEHE to account for the use of a surrogate truth (influence functions method). Also compared using XGBoost or AutoML (fitting many ML methods and choosing the best predictor, or ensembling). Also compared different crossfitting approaches (fitting surrogate models on the same data or the hold-out data) or fitting surrogates on the entire dataset. Nine candidate models were used for evaluation, all DR learners with various base learners (random forests with different hyperparameter settings, various penalised regressions, and SuperLearners with different combinations of base learners).

Model validation exploration tested the use of an interim analysis to explore HTE, varying the interim proportion between 0.25 and 0.75 ($n=1000$). At interim analysis, fit a classifier on the estimated CATEs to identify the top 10% of responders. Then in the remainder of the trial, use the trained classifier to generate a subgroup label and perform a subgroup interaction analysis.

#### Estimand

We are interested in the conditional average treatment effect:
$$ \tau(x) = E[Y(1) - Y(0) | X = x] $$
where $Y(W)$ is the potential outcome under treatment $W$. For binary outcomes, this is the expected risk difference; for continuous outcomes, this is the expected mean difference.

#### Performance Metrics

I computed soooo many different performance metrics across the different components of the simulation study. I've been finding it difficult to summarise them, and the metrics differ slightly from using simulation studies in RCTs because the parameter of interest is a vector (individual-level) rather than a single value.

The precision in estimation of heterogeneous effects (PEHE) is basically an MSE metric for CATEs specifically:
$$ PEHE = \frac{1}{n} \sum_{i=1}^n  (\hat{\tau}(x_i) - \tau(x_i))^2$$

I'd been calling it MSE as I went through the simulations, but it can be called either MSE or PEHE. Representing the square root of this is sometimes better when the MSEs are large, to improve the scaling.

Representing bias has been a bit challenging. So far, I've been summarising bias in each iteration as:
$$ bias = \frac{1}{n} \sum_{i=1}^n (\hat{\tau}(x_i) - \tau(x_i))$$

or relative bias:
$$ bias = \frac{1}{n} \sum_{i=1}^n \frac{(\hat{\tau}(x_i) - \tau(x_i))}{\tau(x_i)}$$

But this can blow up when the $\tau(x_i)$ are close to 0.

In the missing data simulations, there are additional relative PEHE and relative bias metrics that compare the PEHE and bias to the complete data setting. I'm not sure how helpful these metrics really are.

I computed a Pearson correlation metric between $\tau(x_i)$ and $\hat{\tau}(x_i)$ and a Spearman rank correlation. Additionally, I computed the sign accuracy of CATE estimates, estimating the proportion of observations where $\hat{\tau}(x_i)$ has the correct sign.

For the HTE tests, I extracted the p-values from each of the tests and so can look at the distribution of p-values or the power of the tests using a 10% or 5% threshold.

Confidence intervals were assessed using marginal coverage (the proportion of confidence intervals that contained the true CATE) and simultaneous coverage (a binary indicator of whether the confidence intervals contained the true CATE). This was evaluated using the initial trial data (used for model fitting) and a test grid of covariates.

Model evaluation assessed the ranking agreement of the true and estimated model rankings, using Spearman and Kendall correlations. I also looked at the selection accuracy of the top model selected, but I probably won't display that because it was shit across all metrics. I also have the mean true rank of the model selected by each metric. The interesting metric I have is regret, which is the excess PEHE incurred through the metric selecting a model over the actual best one.

For the validation stuff, I have the power of the subgroup interaction test in the stage 2 data. I also had some stuff about the stability of variable importance metrics, but this wasn't very interesting because the scenario also has a single effect modifier. But these simulations run super quickly, so I could explore other scenarios ASAP and have some more interesting results.

### Results

I am not sure whether I should present results from both binary and continuous outcome simulations, or whether it would be better to just focus on one outcome to avoid making the poster too busy.

*update: going to restrict to continuous outcome setting only*

Across all plots, I want to use `theme_light` and set `strip.background = "white"` and `strip.text = "black"`. The colour palette is a bit of a struggle. The Imperial poster has a white background and uses Imperial blue for text (title, authors, logo, etc.), which is #0000cd. There is also a core black, #161a1d, and white, #ffffff.

There is an expanded palette of an additional 27 colours:
- #232333
- #000080
- #8b4513
- #008080
- #c71585
- #4b0082
- #dc143c
- #ff4500
- #006400
- #708090
- #0000cd
- #ffff00
- #40e0d0
- #ee82ee
- #7b68ee
- #ff0000
- #ff8c00
- #00ff7f
- #f5f5f5
- #00bfff
- #f0e68c
- #afeeee
- #ffb6c1
- #e6e6fa
- #fa8072
- #ffa500
- #98fb98

Different plots will have different numbers of colours. Since the background is white, I need to choose colours that stand out against it, and would love it if they were colour-blind friendly, but I don't think this palette is. I don't have to use this palette; I could use any.

#### Sample size requirements

The key takeaway I have is that the forest-based CML methods are unbiased and have pretty low MSE until $n=250$, and at $n=100$ the estimates become crap. When the SuperLearner is used, the $n=250$ estimates are also pretty rubbish.

I'm not sure whether to present the bias/relative bias/MSE/RMSE/Spearman correlation. If presenting bias/MSE, then I think that a lollipop plot would be most useful (the rsimsum package uses lollipop plots and then adds little brackets to show the MCSE confidence intervals). The correlation could also be presented as a lollipop, but I think that 1 would need to be the reference point, as lower correlation means that the model is less good.

*Update: lollipop plot is chaotic as hell. Switch back to a point a line for MCSE confidence intervals (alpha = 0.5)
- x axis: sample size
- y axis:
    - top row: Mean absolute error (MAE)
    - bottom row: Mean PEHE (MSE)
- colour: CATE model
- faceting: scenario (columns)

#### Missing data handling

My main message from this work was that many imputation methods had comparable performance. The computation resource for multiple imputation does not reward improved performance. Missing indicators only work for RF-based methods, and it's not really very principled in the handling of covariates. MPA presents the best trade-off of handling missing data and considering implications of missingness in future predictions. Choice is down to context, and you shouldn't just use inbuilt methods.

Sample size was fixed at $n=500$ for this comparison. I have results from each missing data method used and the CATE model it was used in, but I think that, to make this an understandable plot for a poster, I might need to average across all the CATE models and keep the comparison to the missing data handling method used. Some of the methods aren't used across all the learners (e.g. MIA, only applicable for RF-type models), so the numbers of iterations vary slightly.

I think the lollipop plots would be best again, for bias/relative bias/relative complete bias/etc.

*update: no lollipop. go back to point and confidence intervals, average over all models and use colour to show missingness handling approach.*
- x axis: missing data mechanism, MAR, MNAR, MNAR-Y
- y axis:
    - top row: MAE
    - bottom row: MSE
- colour: missing data handling method
- faceting: scenario (cols)

#### HTE testing methods

Key message: extremely variable performance means that these can't be the only thing used to make inference on HTE. It needs to be a multidisciplinary discussion, taking into account the context and the MCID. It can be a useful piece of evidence when exploring HTEs in more detail, e.g. along covariates.

Again, I think it would be better to average across all CATE models and focus on comparisons between testing procedures. I also tested the performance of the tests when they had the true CATEs. These are the true BLP results and the DR-oracle PO permutation test. I am not sure whether I could include the performance of the tests with both estimates and true values in the same plot, or if I should just show the performance of the tests when they have access to the true data.

I think that power (or the false positive rate in the null scenario) would be the best summary metric. I can then use lollipop plots showing deviation from desired power (0.05 in the null setting and 0.95 in the other settings, or 0.1/0.9).
- x axis: sample size
- y axis: power/type 1 error
- colour: test type
- faceting: scenario

#### Confidence interval coverage/uncertainty estimation

Sidenote: I'm not sure if I should be computing a bias-corrected coverage estimate, because some of the methods did display bias, and their poor coverage is more likely due to the bias in the point estimates than the width of the confidence intervals.

Key message: half-sample bootstrap is computationally intensive and only available for RF-based methods, and you need to carefully consider hyperparameter settings to make sure that you get the desired coverage. CF-based methods also have dramatic undercoverage when HTE is more complex, but, more often than not, confidence intervals are overly conservative.

It might be best to show the simultaneous coverage. But when the goal is to explore HTE within a single dataset, the marginal coverage might be more interesting... Use a lollipop plot again to show deviation from desired coverage.

- x axis: subsampling ratio
- y axis: simultaneous/marginal coverage (not sure whether this should be in the covariate grid or the observations)
- colour: method - CF variance estimates or half-sample bootstrap
- faceting: scenario

#### Model evaluation

The key takeaway - influence-corrected metrics generally perform less well at model selection - bias correction can make the metric negative, and then it's not clear how to proceed... But even the DR score is not great at selecting the absolute best CATE model. You're never gonna know the actual best model, and will probably have a set of candidate CATE models with comparable performance. The DR score could therefore be a good way of removing the bad candidates and enable a researcher to choose their final CATE model with other considerations in mind - interpretability etc. Heavy ML is required to get surrogate CATEs, but AutoML can work in small-ish data (n=250), just not when you only estimate surrogates in the hold-out data - you need the full data, or to use the same data that you did to fit the candidates, which is kind of bad practice. This would be most useful in scenarios where you have multiple trials or something like that, to explore external validity better.

I think the regret is the most interesting thing to show here, but this simulation study was way more factorial, so I could show an axis plot, which compares each factor and marginalises out the rest. I also only want to show the metrics where the propensity was fixed at 0.5 and not estimated. A lollipop plot would be good for showing the regret, but the true rank plot of the various different modelling configurations is the most visually appealing (using a viridis colour scale).

lollipop:
- x axis: sample size
- y axis: regret
- colour: surrogate model (AutoML / XGBoost) x metrics (DR score / IF score)
- facet: scenario
- shape: crossfitting method

true rankings:
Would it be bad to show the distribution of true rankings in the axis way? As in, comparing all XGBoost vs AutoML, all DR vs IF, etc.?

#### Subgroup validation

Key takeaway: ML-assisted subgroup exploration could enhance subgroup analysis in RCTs by exploring HTE early on (which we know is possible now that we've established sample size requirements). We could define subgroups in a complex manner, or use variable importance to refine subgroup analyses...

This is what I also want to focus on in the future work section.

I want a plot which shows the proportion of successful stage 2 interaction tests when the subgroup was defined using a classifier trained on the stage 1 data. I already have this somewhere, and it's a super simple plot at the minute, because there's only one model and one subgroup being defined...

### Discussion

We've established some parameters for where and when we can use ML to explore HTEs in RCTs... Other general frameworks have been developed, but they have not considered requirements for sample size and missing data handling... This work has shown that medium-sized trials can leverage CML to gain a better picture of HTEs...

Uncertainty estimates at the individual level are not great, and HTE tests are also lacking, and so further development of these is required if we want to explore HTE more formally and use individual-level estimates in practice...

Internal validation work shows potential utility of CML within trials, and potentially within adaptive trials... We would like to talk about enrichment designs and sample size re-estimation to estimate subgroup effects, and potentially even covariate response adaptive randomisation stuff...