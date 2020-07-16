#' Regression models for localised mutations: Evaluating differential mutation rates across classes of sites
#' 
#' RM2() uses negative binomial regression to evaluate local mutation rates and processes between stacked sites of the same class to flanking control regions
#' @param maf Data frame of mutations prepared by get_mut_trinuc_strand
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{start}{the start position of the mutation in base 1 coordinates}
#'     \item{end}{the end position of the mutation in base 1 coordinates}
#'     \item{ref}{the reference allele as a string containing the bases A, T, C or G}
#'     \item{alt}{the alternate allele as a string containing the bases A, T, C or G}
#'     \item{mut_trinuc}{trinucleotide context - where middle is C or T - with alternate allele}
#'     \item{mut_strand}{character indicating Watson (w) or Crick (c)}
#'     \item{ref_alt}{character indicating single-base substitution}
#' }
#' @param sites Data frame
#' \describe{
#'     \item{chr}{autosomal chromosomes as chr1 to chr22 and sex chromosomes as chrX and chrY}
#'     \item{start}{the start position of the mutation in base 0 coordinates}
#'     \item{end}{the end position of the mutation in base 0 coordinates}
#' }
#' @param mut_class_column Character corresponding to column of mutation classes for grouped analysis
#' @param cofactor_column Character corresponding to column of cofactors
#' @param window_size Integer indicating the half-width of sites and flanking regions (added to left and right for full width). Default 100bp
#' @param n_min_mut Integer indicating the minimum number of mutations required to perform analysis
#'
#' @return Data frame containing the regression estimates and likelihood ratio test output with the following columns: mut_type,
#' pp, this_coef, obs_mut, exp_mut, exp_mut_hi, exp_mut_lo, fc, n_sites_tested
#' \describe{
#'     \item{mut_type}{A string identifying the mutation class}
#'     \item{pp}{The p-value from the likelihood ratio test}
#'     \item{this_coef}{The coefficient from is_site}
#'     \item{obs_mut}{The total number of observed mutations of that class}
#'     \item{exp_mut}{The expected number of mutations determined by the model}
#'     \item{exp_mut_lo}{Lower bound of 95% confidence interval}
#'     \item{exp_mut_hi}{Upper bound of 95% confidence interval}
#'     \item{fc}{Observed mutations divided by expected mutations}
#'     \item{pp_cofac}{The p-value from the likelihood ratio test of site:cofactor interaction}
#'     \item{this_coef_cofac}{The coefficient from the site:cofactor interaction term}
#'     \item{n_sites_tested}{The number of sites that were tested - all sites if no downsampling}
#' }
#' @export
RM2 = function(maf, sites, mut_class_columns = NA, cofactor_column = NA,
			window_size = 100, n_min_mut = 100) {

	prepared_sites = NULL
	if (is.na(cofactor_column)) {
		prepared_sites =
			.prepare_bins_of_sites(sites, window_size, maf,
						n_bin = 10, n_min_mut = n_min_mut)
	} else {
		split_maf = split(maf, maf[,cofactor_column])
		prepared_sites = lapply(names(split_maf), function(cof)
				.prepare_bins_of_sites(sites, window_size, split_maf[[cof]],
						n_bin = 10, n_min_mut = n_min_mut))
		names(prepared_sites) = names(split_maf)
	}

	maf_list = do.call(c, lapply(mut_class_columns, function(mc) {
		if (is.na(mc)) {
			list("total_muts__total_muts" = maf)
		} else {
			res = split(maf, maf[,mc])
			names(res) = paste(mc, names(res), sep="__")
			res
		}
	}))

	# maf1_title = "total_muts__total_muts"
	res = do.call(rbind, lapply(names(maf_list), function(maf1_title) {
		maf1 = maf_list[[maf1_title]]

		dfr = NULL
		if (is.na(cofactor_column)) {
			dfr = .maf_to_dfr(maf1, prepared_sites)
		} else {
			split_maf1 = split(maf1, maf1[,cofactor_column])
			dfr = do.call(rbind, lapply(names(split_maf1), function(cof)
					data.frame(.maf_to_dfr(split_maf1[[cof]], prepared_sites), cof, stringsAsFactors=F)))
		}

		.test_NB_model(dfr, maf1_title, n_min_mut = n_min_mut, test_cofactor = !is.na(cofactor_column))
	}))

	if (!is.null(res[[1]])) {
		res$n_sites_tested = nrow(sites)
	}

	res
}



#' Runs negative binomial regression and likelihood ratio test
#'
#' @param dfr Data frame generated from maf_to_dfr that contains mutation and trinucleotide counts for each quadnucleotide context + indel
#' @param mut_type Character corresponding to total mutations or name of subclass
#' @param n_min_mut Integer (100) representing the minimum number of mutations required to run test
#' @param test_cofactor Boolean indicating whether to analyze cofactor
#'
#' @return Data frame containing the results of the regression with the following columns: mut_type,
#' pp, this_coef, obs_mut, exp_mut, exp_mut_lo, exp_mut_hi, fc, pp_cofac, this_coef_cofac
#' \describe{
#'     \item{mut_type}{A string identifying the mutation class}
#'     \item{pp}{The p-value from the likelihood ratio test}
#'     \item{this_coef}{The coefficient from is_site}
#'     \item{obs_mut}{The total number of observed mutations of that class}
#'     \item{exp_mut}{The expected number of mutations determined by the model}
#'     \item{exp_mut_lo}{Lower bound of 95% confidence interval}
#'     \item{exp_mut_hi}{Upper bound of 95% confidence interval}
#'     \item{fc}{Observed mutations divided by expected mutations}
#'     \item{pp_cofac}{The p-value from the likelihood ratio test of site:cofactor interaction}
#'     \item{this_coef_cofac}{The coefficient from the site:cofactor interaction term}
#' }
.test_NB_model = function(dfr, mut_type, n_min_mut = 100, test_cofactor = F) {

	empty_res = data.frame(mut_type, pp=NA, this_coef=NA, obs_mut=NA,
		exp_mut=NA, exp_mut_lo=NA, exp_mut_hi=NA, fc=NA, pp_cofac=NA, this_coef_cofac=NA,
		stringsAsFactors=F)

	if (test_cofactor & length(unique(dfr$cof)) != 2) {
		print(paste0(mut_type,  ": Incorrect number of cofactor values (",
				length(unique(dfr$cof)), "), skipping."))
		return(empty_res)
	}

	if (sum(dfr$muts_count) < n_min_mut) {
		print(paste0(mut_type,  ": Too few mutations in sites+flanks (", sum(dfr$muts_count), "), skipping."))
		return(empty_res)
	}

	silent_quadnucs = names(which(c(by(dfr$muts_count, dfr$quadnuc, sum))==0))
	dfr = dfr[!dfr$quadnuc %in% silent_quadnucs, ]

	this_formula = "muts_count ~ log1p(mut_rate) + is_site + offset(log(posi_count))"

	# adjust formula if only one quad-nuc context of mutations is available; most commonly indel models
	if (length(unique(dfr$quadnuc)) > 1)	{
		this_formula = paste(this_formula, "+ quadnuc")
	}

	#
	if (test_cofactor)	{
		this_formula = paste(this_formula, "+ cof")
	}

	h1 = try(MASS::glm.nb(stats::as.formula(this_formula), data=dfr), silent=T)

	if (any(class(h1)=="try-error")) {
		print(paste0(mut_type, ": glm.nb failed at h1, n_mut=", sum(dfr$muts_count), "."))
		return(empty_res)
	}

	h0 = try(stats::update(h1, . ~ . - is_site, data=dfr), silent=T)

	if (any(class(h0)=="try-error")) {
		print(paste0(mut_type, ": glm.nb failed at h0, n_mut=", sum(dfr$muts_count), "."))
		return(empty_res)
	}

	# "manual" calculation is more precise, thanks to lower.tail.
	# if using 1-p(), again limited to 10-16.
	# anova(h1) provides wald test that is less conservative
	pp = stats::pchisq(-(h0$twologlik - h1$twologlik), 1, lower.tail=F)
	this_coef = stats::coef(h1)[['is_siteTRUE']]

	# compute predicted mutation counts with ranges
	h0_mut = .predict_muts_ranges_probs(h0, which(dfr$is_site))
	h1_mut = .predict_muts_ranges_probs(h1, which(dfr$is_site))
	real_obs_mut = sum(dfr[dfr$is_site,"muts_count"])
	exp_mut = h0_mut[['mid']]
	exp_mut_lo = h0_mut[['lo']]
	exp_mut_hi = h0_mut[['hi']]
	fc = h1_mut[['mid']] / h0_mut[['mid']]

	pp_cofac = this_coef_cofac = NA
	if (test_cofactor) {

		h2 = try(stats::update(h1, . ~ . + cof:is_site, data=dfr), silent=T)

		pp_cofac = stats::pchisq(-(h1$twologlik - h2$twologlik), 1, lower.tail=F)
		# coefficient is computed with first level of factor as reference
		this_coef_cofac_name = paste0("is_siteTRUE:cof", sort(unique(dfr$cof))[2])
		this_coef_cofac = stats::coef(h2)[[this_coef_cofac_name]]
	}

	data.frame(mut_type, pp, this_coef, obs_mut=real_obs_mut,
		exp_mut, exp_mut_lo, exp_mut_hi, fc, pp_cofac, this_coef_cofac, stringsAsFactors=F)
}


#' Extract expected mutations count estimates by probability sampling
#'
#' @param hyp List containing glm.nb regression outputs
#' @param select_positions Numeric vector containing indices of site positions
#'
#' @return Vector of mid, lo and hi mutation count predictions
.predict_muts_ranges_probs = function(hyp, select_positions) {
  hyp_predict_probs = replicate(1000, sum(MASS::rnegbin(stats::fitted(hyp), theta = (summary(hyp))$theta)[select_positions]))
	mid = stats::median(hyp_predict_probs)
	lo = stats::quantile(hyp_predict_probs, 0.025, na.rm=T)[[1]]
	hi = stats::quantile(hyp_predict_probs, 0.975, na.rm=T)[[1]]

	c(mid=mid, lo=lo, hi=hi)
}
